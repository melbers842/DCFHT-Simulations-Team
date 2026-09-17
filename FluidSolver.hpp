#pragma once
//
// FluidSolver.hpp
//
// Multi-dimensional generalisation of the 1-D Hall thruster discharge model,
// Eq. 13.22 - 13.26 of Lecture 14 & 15.
//
// WHAT CHANGED FROM THE LECTURE
// -----------------------------
// The lecture solves a steady 1-D system by reformulating in ion-flux space
// (Eq. 13.27 - 13.32) to avoid the sonic singularities, then integrating from
// the anode. That reformulation does not survive the move to 2-D: ion flux
// stops being usable as an independent variable and the "constants of the
// flow" stop being constant along a single path.
//
// So this solver goes back to the underlying conservation laws in vector form
// and marches them in time until they stop changing. Steady state is then the
// converged solution rather than the result of an integration.
//
//   continuity   d n_s/dt + div(n_s v_s) = +/- n_e nu_i        (from 13.22)
//   ion mom.     see ionMomentumRHS()                          (from 13.24)
//   neutral mom. see neutralMomentumRHS()                      (from 13.25)
//   electron mom. inertialess drift-diffusion, see electronVelocity()  (13.23)
//   electron en. see electronEnergyRHS()                       (from 13.26)
//
// The electron momentum equation keeps magnetisation, which the 1-D model does
// not need. That anisotropy is the physics the DCFHT study is actually about.
//
// STILL OPEN -- these are yours to decide, see the SRD
//   * potential closure: Boltzmann (implemented) vs full Poisson solve
//   * anomalous transport coefficient alpha
//   * wall boundary conditions
//
// DCFHT Simulation Team
//

#include "Collisions.hpp"
#include "Constants.hpp"
#include "Mesh.hpp"
#include "PlasmaState.hpp"
#include "Vec.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

template <int D>
class FluidSolver {
public:
    // =====================================================================
    // How the electric potential is obtained. This is a real modelling
    // decision, not a tuning knob, and the two options are not interchangeable.
    //
    // Imported  -- phi is read from file and held fixed. This matches the plan
    //              of importing field data from ANSYS. The electron drift-
    //              diffusion equation is then well-posed: a known potential
    //              drives a computed velocity.
    //
    // Boltzmann -- phi = phi_ref + T_e ln(n_e/n_ref), from Eq. 13.23 with the
    //              collisional drag dropped.
    //
    //              WARNING, AND THIS IS THE IMPORTANT PART: Boltzmann must NOT
    //              be combined with computing v_e from the drift-diffusion
    //              residual. The Boltzmann relation is derived by ASSUMING the
    //              pressure gradient balances the electric field, so by
    //              construction e n_e grad(phi) is very nearly equal to
    //              grad(n_e k T_e). Computing v_e from the difference of those
    //              two terms means computing a small difference of two large
    //              nearly-equal numbers, then multiplying it by a mobility of
    //              order 10^3. The result is amplified round-off, which feeds
    //              the ohmic term of the energy equation, which raises T_e,
    //              which raises phi, which raises grad(phi). That loop diverges
    //              in single-digit timesteps regardless of how small dt is.
    //
    //              If you want a self-consistent potential, the physically
    //              correct closure is current continuity, div(j) = 0 with
    //              j = e(n_i v_i - n_e v_e). That determines v_e from charge
    //              conservation rather than from a cancelling residual, and it
    //              turns phi into an elliptic solve -- which is what Eigen's
    //              sparse solvers are in the dependency list for.
    // =====================================================================
    enum class PotentialClosure {
        Imported,
        Boltzmann
    };

    // =====================================================================
    // How the electron velocity is obtained.
    //
    // THIS IS THE CENTRAL CONSTRAINT OF A QUASINEUTRAL FLUID HALL CODE, and it
    // is worth stating plainly because getting it wrong produces electron
    // velocities faster than light rather than an obvious error message.
    //
    // Electron velocity and potential cannot BOTH be free. In a quasineutral
    // plasma they are tied together by current continuity:
    //
    //     j = e (n_i v_i - n_e v_e)        div(j) = 0
    //
    // So exactly one of the two is prescribed and the other follows:
    //
    //   CurrentContinuity -- the discharge current is prescribed (which is what
    //       Eq. 13.68 does at the anode) and electron flux follows from
    //       n_e v_e = n_i v_i - j/e. Well-posed with an imported potential.
    //       This is the default and what the demo uses.
    //
    //   DriftDiffusion -- v_e from the mobility form of Eq. 13.23. Only valid
    //       when phi is itself solved from div(j) = 0. Combined with a
    //       PRESCRIBED phi it is unconstrained: nothing stops the mobility,
    //       which reaches 10^5 m^2/Vs where B is weak, from multiplying a
    //       5 x 10^4 V/m field into a 10^9 m/s velocity. Physically, real
    //       electrons are held back by the space charge that would build up;
    //       with phi fixed, that feedback is switched off.
    // =====================================================================
    enum class ElectronVelocityClosure {
        CurrentContinuity,
        DriftDiffusion
    };

    struct Parameters {
        PotentialClosure closure = PotentialClosure::Imported;
        ElectronVelocityClosure electronClosure = ElectronVelocityClosure::CurrentContinuity;
        // Discharge current density, A/m^2, for the CurrentContinuity closure.
        // Direction is the axial (last) axis. Order 10^3 for a small thruster;
        // cf. the ion flux panel of Fig. 13.4, where e * Gamma_i ~ 10^3 A/m^2.
        double dischargeCurrentDensity = 1040.0;
        double anomalousAlpha   = 1.0 / 16.0;  // Bohm coefficient, sweep this
        double referenceDensity = 1.0e17;      // m^-3, Boltzmann closure only
        double referencePhi     = 0.0;         // V, anode potential (Eq. 13.67: zero)
        double floorDensity     = 1.0e10;      // m^-3, keeps divisions finite
        double floorTemperature = 0.1;         // eV, same reason
        double maxTemperature   = 100.0;       // eV, numerical guard
        double cflNumber        = 0.05;        // advection limit
    };

    FluidSolver(const Mesh<D>& mesh, Parameters params)
        : mesh_(mesh), p_(params) {}

    // ---------------------------------------------------------------------
    // Advance the whole state by dt using forward Euler.
    //
    // Forward Euler is chosen deliberately as a starting point: it is easy to
    // read against the equations above and easy to debug. It is also only
    // first-order accurate and conditionally stable, so once the physics is
    // verified this is the obvious thing to replace with RK2/RK4 or an
    // implicit treatment of the stiff electron energy term.
    // ---------------------------------------------------------------------
    void step(PlasmaState<D>& s, double dt) const {
        const auto n = static_cast<std::size_t>(mesh_.numNodes());

        // Work on copies so every update sees the same old state, rather than
        // a half-updated one. Cheap at these mesh sizes and removes a whole
        // category of ordering bug.
        std::vector<double> dn_i(n, 0.0), dn_n(n, 0.0), dEnergy(n, 0.0);
        std::vector<Vec<D>> dv_i(n, Vec<D>{}), dv_n(n, Vec<D>{});

        for (std::size_t k = 0; k < n; ++k) {
            const int i = static_cast<int>(k);
            dn_i[k]    = ionContinuityRHS(s, i);
            dn_n[k]    = neutralContinuityRHS(s, i);
            dv_i[k]    = ionMomentumRHS(s, i);
            dv_n[k]    = neutralMomentumRHS(s, i);
            dEnergy[k] = electronEnergyRHS(s, i);
        }

        for (std::size_t k = 0; k < n; ++k) {
            // --- densities, Eq. 13.22 ---
            s.n_i[k] = std::max(s.n_i[k] + dt * dn_i[k], p_.floorDensity);
            s.n_n[k] = std::max(s.n_n[k] + dt * dn_n[k], p_.floorDensity);

            // --- heavy-species velocities, Eq. 13.24 and 13.25 ---
            s.v_i[k] += dv_i[k] * dt;
            s.v_n[k] += dv_n[k] * dt;

            // --- electron energy, Eq. 13.26 ---
            // Stored as the volumetric energy density 3/2 n_e k T_e, advanced,
            // then converted back to a temperature.
            const double oldEnergy = 1.5 * s.n_e[k] * constants::kTe_J(s.T_e[k]);
            const double newEnergy = std::max(oldEnergy + dt * dEnergy[k], 0.0);
            const double density   = std::max(s.n_e[k], p_.floorDensity);
            // Bounded on both sides. The upper bound is a guard, not physics:
            // electron temperature in a Hall thruster channel is tens of eV,
            // so hitting this ceiling means something diverged and you want to
            // notice, not silently integrate a runaway.
            s.T_e[k] = std::clamp(newEnergy / (1.5 * density * constants::kElementaryCharge),
                                  p_.floorTemperature, p_.maxTemperature);
        }

        // --- closures, applied after the time-advanced quantities -----------
        applyQuasineutrality(s);
        updatePotential(s);
        updateElectronVelocity(s);
    }

    // ---------------------------------------------------------------------
    // Eq. 13.22 -- continuity.
    //
    //   d n_i/dt + div(n_i v_i) = + n_e nu_i
    //   d n_n/dt + div(n_n v_n) = - n_e nu_i
    //
    // Electron continuity is not integrated separately; quasineutrality ties
    // n_e to n_i. If you later move to a Poisson solve, electron continuity
    // becomes a real equation and this is where it goes.
    // ---------------------------------------------------------------------
    double ionContinuityRHS(const PlasmaState<D>& s, int i) const {
        const auto k = static_cast<std::size_t>(i);
        const double nu_i = collisions::ionizationFrequency(s.n_n[k], s.T_e[k]);
        return -mesh_.divergenceOfFlux(s.n_i, s.v_i, i) + s.n_e[k] * nu_i;
    }

    double neutralContinuityRHS(const PlasmaState<D>& s, int i) const {
        const auto k = static_cast<std::size_t>(i);
        const double nu_i = collisions::ionizationFrequency(s.n_n[k], s.T_e[k]);
        return -mesh_.divergenceOfFlux(s.n_n, s.v_n, i) - s.n_e[k] * nu_i;
    }

    // ---------------------------------------------------------------------
    // Eq. 13.24 -- ion momentum, divided through by m_i n_i to give dv_i/dt.
    //
    //   dv_i/dt = -(v_i . grad) v_i
    //             - (e/m_i) grad(phi)
    //             - (nu_i + nu_in + nu_cex)(v_i - v_n)
    //
    // Ions are treated as unmagnetised, following the lecture's note that they
    // interact only weakly with the thruster's magnetic field because of their
    // large gyroradius.
    // ---------------------------------------------------------------------
    Vec<D> ionMomentumRHS(const PlasmaState<D>& s, int i) const {
        const auto k = static_cast<std::size_t>(i);
        using namespace constants;

        Vec<D> rhs = advection(s.v_i, s.v_i[k], i) * -1.0;

        const Vec<D> gradPhi = mesh_.gradient(s.phi, i);
        rhs += gradPhi * (-kElementaryCharge / kXenonMass);

        const Vec<D> slip = s.v_i[k] - s.v_n[k];
        const double relSpeed = slip.norm();
        const double nu_i   = collisions::ionizationFrequency(s.n_n[k], s.T_e[k]);
        const double nu_in  = collisions::ionNeutralFrequency(s.n_n[k], relSpeed);
        const double nu_cex = collisions::chargeExchangeFrequency(s.n_n[k], relSpeed);
        rhs += slip * -(nu_i + nu_in + nu_cex);

        return rhs;
    }

    // ---------------------------------------------------------------------
    // Eq. 13.25 -- neutral momentum, divided through by m_i n_n.
    //
    //   dv_n/dt = -(v_n . grad) v_n
    //             + (nu_in + nu_cex)(v_i - v_n)
    //             - (k T_n / (m_i n_n)) grad(n_n)
    //
    // T_n is held at 300 K, so the pressure gradient collapses from
    // grad(n_n k T_n) to k T_n grad(n_n).
    // ---------------------------------------------------------------------
    Vec<D> neutralMomentumRHS(const PlasmaState<D>& s, int i) const {
        const auto k = static_cast<std::size_t>(i);
        using namespace constants;

        Vec<D> rhs = advection(s.v_n, s.v_n[k], i) * -1.0;

        const Vec<D> slip = s.v_i[k] - s.v_n[k];
        const double relSpeed = slip.norm();
        const double nu_in  = collisions::ionNeutralFrequency(s.n_n[k], relSpeed);
        const double nu_cex = collisions::chargeExchangeFrequency(s.n_n[k], relSpeed);
        rhs += slip * (nu_in + nu_cex);

        const double density = std::max(s.n_n[k], p_.floorDensity);
        const Vec<D> gradN   = mesh_.gradient(s.n_n, i);
        rhs += gradN * (-kTn_J() / (kXenonMass * density));

        return rhs;
    }

    // ---------------------------------------------------------------------
    // Eq. 13.26 -- electron energy, as a volumetric energy density rate.
    //
    //   d/dt(3/2 n_e k T_e) = -div(5/2 k T_e n_e v_e)
    //                         + e n_e v_e . grad(phi)
    //                         - n_e nu_i E_i
    //
    // Returns J/m^3/s.
    // ---------------------------------------------------------------------
    double electronEnergyRHS(const PlasmaState<D>& s, int i) const {
        const auto k = static_cast<std::size_t>(i);
        using namespace constants;

        // Convective transport of enthalpy, div(5/2 k T_e n_e v_e).
        // Built as a scalar field times the electron velocity so the same
        // conservative flux-divergence routine can be reused.
        std::vector<double> enthalpy(s.n_e.size());
        for (std::size_t j = 0; j < s.n_e.size(); ++j) {
            enthalpy[j] = 2.5 * kTe_J(s.T_e[j]) * s.n_e[j];
        }
        const double transport = -mesh_.divergenceOfFlux(enthalpy, s.v_e, i);

        // Ohmic heating, e n_e v_e . grad(phi).
        const Vec<D> gradPhi = mesh_.gradient(s.phi, i);
        const double ohmic   = kElementaryCharge * s.n_e[k] * s.v_e[k].dot(gradPhi);

        // Ionization cost, n_e nu_i E_i.
        const double nu_i = collisions::ionizationFrequency(s.n_n[k], s.T_e[k]);
        const double ionizationLoss =
            s.n_e[k] * nu_i * kXenonIonizationEnergy * kElementaryCharge;

        return transport + ohmic - ionizationLoss;
    }

    // ---------------------------------------------------------------------
    // Eq. 13.23 -- electron momentum with electron inertia neglected.
    //
    // Rearranged for velocity:
    //
    //   v_e = mu_e grad(phi) - (mu_e/(e n_e)) grad(n_e k T_e)
    //
    // with mu_e = e/(m_e nu_e) the unmagnetised mobility. Magnetisation is
    // then applied by splitting that drive into components parallel and
    // perpendicular to B: the parallel component is unaffected, the
    // perpendicular one is reduced by 1/(1 + Omega^2) where Omega = w_c/nu_e
    // is the Hall parameter.
    //
    // That anisotropy is the whole point of the DCFHT study. Where B is strong
    // and nearly parallel to a wall, Omega is large, cross-field transport is
    // choked, and electrons are confined -- the mirroring behaviour you are
    // trying to characterise. It also works unchanged in 3-D, which is why
    // this is done with a projection rather than an explicit cross product.
    //
    // UPGRADE PATH: this is a scalar reduction, not the full mobility tensor.
    // It omits the Hall (E x B) component of the drift, which is real and is
    // what carries the azimuthal Hall current. Add the full tensor once the
    // scalar version is verified.
    // ---------------------------------------------------------------------
    Vec<D> electronVelocity(const PlasmaState<D>& s, int i) const {
        const auto k = static_cast<std::size_t>(i);
        using namespace constants;

        const double density = std::max(s.n_e[k], p_.floorDensity);

        // --- current continuity closure ------------------------------------
        // n_e v_e = n_i v_i - j/e, i.e. electrons carry whatever current the
        // ions do not. A uniform axial j is divergence-free by construction,
        // so div(j) = 0 is satisfied exactly. When you move to a non-uniform
        // current this becomes an elliptic solve for phi -- that is the point
        // at which Eigen's sparse solvers earn their place in the build.
        if (p_.electronClosure == ElectronVelocityClosure::CurrentContinuity) {
            Vec<D> flux = s.v_i[k] * s.n_i[k];
            flux[D - 1] -= p_.dischargeCurrentDensity / kElementaryCharge;
            return flux * (1.0 / density);
        }

        // --- drift-diffusion closure ---------------------------------------
        const double B       = s.B[k].norm();
        const double nu_e    = std::max(
            collisions::electronMomentumFrequency(s.n_n[k], s.T_e[k], B, p_.anomalousAlpha),
            1.0);
        const double mobility = kElementaryCharge / (kElectronMass * nu_e);

        // Electron pressure gradient, grad(n_e k T_e).
        std::vector<double> pressure(s.n_e.size());
        for (std::size_t j = 0; j < s.n_e.size(); ++j) {
            pressure[j] = s.n_e[j] * kTe_J(s.T_e[j]);
        }
        const Vec<D> gradP   = mesh_.gradient(pressure, i);
        const Vec<D> gradPhi = mesh_.gradient(s.phi, i);

        Vec<D> drive = gradPhi * mobility
                     + gradP * (-mobility / (kElementaryCharge * density));

        // No field, no magnetisation.
        if (B <= 0.0) return drive;

        const double hall = collisions::cyclotronFrequency(B) / nu_e;
        const double perpFactor = 1.0 / (1.0 + hall * hall);

        Vec<D> bHat = s.B[k] * (1.0 / B);
        const double along = drive.dot(bHat);
        const Vec<D> parallel = bHat * along;
        const Vec<D> perpendicular = drive - parallel;

        return parallel + perpendicular * perpFactor;
    }

    // ---------------------------------------------------------------------
    // Timestep limits.
    //
    // Returns the largest stable dt for the heavy species, from a CFL
    // condition on advection. Call this every step -- as the discharge builds
    // up, ion velocities rise and the allowable step shrinks.
    //
    // NOTE: this does NOT bound the electron energy equation, which is much
    // stiffer and will usually be the real constraint. If runs blow up while
    // satisfying this limit, the energy equation is why, and the fix is to
    // sub-cycle it or make it implicit rather than to shrink dt globally.
    // ---------------------------------------------------------------------
    double maxStableTimestep(const PlasmaState<D>& s) const {
        double fastest = 1.0;
        for (std::size_t k = 0; k < s.n_e.size(); ++k) {
            fastest = std::max(fastest, s.v_i[k].norm());
            fastest = std::max(fastest, s.v_n[k].norm());
            // Electrons are by far the fastest species and appear in the
            // enthalpy transport term of the energy equation, so they set the
            // real limit -- typically two to three orders of magnitude below
            // the heavy-species limit. Leaving them out of this calculation
            // produces a dt that looks reasonable and blows up in a few steps.
            fastest = std::max(fastest, s.v_e[k].norm());
        }
        return p_.cflNumber * mesh_.minSpacing() / fastest;
    }

    // Electron plasma period, for reference. If you later replace the
    // Boltzmann closure with a Poisson solve, the timestep must resolve this,
    // and it is typically several orders of magnitude smaller than the CFL
    // limit above -- which is the main practical cost of going self-consistent.
    static double plasmaPeriod(double n_e) {
        using namespace constants;
        if (n_e <= 0.0) return 1.0;
        const double omega_pe = std::sqrt(n_e * kElementaryCharge * kElementaryCharge
                                          / (kEpsilon0 * kElectronMass));
        return 2.0 * kPi / omega_pe;
    }

private:
    // (v . grad) v, the advective term shared by both momentum equations.
    Vec<D> advection(const std::vector<Vec<D>>& field, const Vec<D>& velocity, int i) const {
        Vec<D> result;
        for (int comp = 0; comp < D; ++comp) {
            std::vector<double> component(field.size());
            for (std::size_t j = 0; j < field.size(); ++j) component[j] = field[j][comp];
            result[comp] = velocity.dot(mesh_.gradient(component, i));
        }
        return result;
    }

    // n_e = n_i. Replace this with a Poisson solve if you go self-consistent.
    void applyQuasineutrality(PlasmaState<D>& s) const {
        for (std::size_t k = 0; k < s.n_e.size(); ++k) s.n_e[k] = s.n_i[k];
    }

    // Boltzmann potential: Eq. 13.23 with the collisional drag dropped and
    // integrated, giving phi = phi_ref + T_e ln(n_e/n_ref).
    //
    // T_e is in eV, so this comes out in volts with no extra conversion --
    // one of the few places the eV convention pays for itself directly.
    void updatePotential(PlasmaState<D>& s) const {
        // Imported potential is held fixed -- nothing to recompute.
        if (p_.closure == PotentialClosure::Imported) return;

        // The log is clamped because the Boltzmann relation is only valid
        // where there is actually plasma. Left unclamped, a node sitting at
        // the density floor gives ln(1e10/3e17) ~ -17, a 50 V swing across one
        // cell. The clamp bounds the excursion to +/- kMaxLogRatio * T_e volts.
        //
        // A numerical guard, not physics. If results depend on its value, the
        // mesh is reaching into vacuum and the domain needs attention instead.
        constexpr double kMaxLogRatio = 4.0;
        for (std::size_t k = 0; k < s.n_e.size(); ++k) {
            const double ratio = std::max(s.n_e[k], p_.floorDensity) / p_.referenceDensity;
            const double logRatio = std::clamp(std::log(ratio), -kMaxLogRatio, kMaxLogRatio);
            s.phi[k] = p_.referencePhi + s.T_e[k] * logRatio;
        }
    }

    void updateElectronVelocity(PlasmaState<D>& s) const {
        std::vector<Vec<D>> updated(s.v_e.size());
        for (std::size_t k = 0; k < s.v_e.size(); ++k) {
            updated[k] = electronVelocity(s, static_cast<int>(k));
        }
        s.v_e = std::move(updated);
    }

    const Mesh<D>& mesh_;
    Parameters     p_;
};
