//
// fluid_demo_main.cpp
//
// Drives the 2-D fluid solver on a simplified DCFHT channel and writes a .vti
// time series for ParaView.
//
// WHAT THIS IS
//   A wiring harness. Every piece is connected and running: mesh, state,
//   equations 13.22-13.26, output. The geometry, the magnetic field and the
//   boundary conditions are all simplified stand-ins.
//
// WHAT THIS IS NOT
//   A validated Hall thruster simulation. The magnetic field is an analytic
//   mirror-like guess, not the ANSYS import. The boundary conditions are
//   crude. The cross sections in Collisions.hpp are placeholders. Do not read
//   physics off this yet.
//
// The value is that it runs end to end, so from here every change is a
// physics change against a working pipeline.
//
// Change kDim to 3 and the whole thing compiles and runs in three dimensions.
//
// DCFHT Simulation Team
//

#include "Constants.hpp"
#include "FluidSolver.hpp"
#include "Mesh.hpp"
#include "PlasmaState.hpp"
#include "Vec.hpp"
#include "vti_writer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

// ===========================================================================
// The one line that switches between 2-D and 3-D.
// ===========================================================================
constexpr int kDim = 2;

// --- geometry, loosely BHT-200 scale (see Fig. 13.4) -----------------------
constexpr double kChannelLength      = 0.025;   // m, axial extent
constexpr double kChannelWidth       = 0.010;   // m, radial extent
constexpr int    kNodesRadial        = 24;
constexpr int    kNodesAxial         = 60;

// --- operating conditions --------------------------------------------------
constexpr double kAnodeNeutralDensity = 1.0e20;  // m^-3
constexpr double kSeedPlasmaDensity   = 3.0e17;  // m^-3, cf. Eq. 13.67 range
constexpr double kSeedElectronTemp    = 3.0;     // eV, cf. Eq. 13.67 range
constexpr double kNeutralInjectSpeed  = 300.0;   // m/s, roughly sonic at 300 K
constexpr double kPeakB               = 0.02;    // T, cf. SPT-100 example
constexpr double kDischargeVoltage    = 250.0;   // V, anode-to-exit drop

// --- run control ------------------------------------------------------------
constexpr int kNumSteps    = 20000;
constexpr int kOutputEvery = 400;    // 50 output frames

// Analytic stand-in for the imported field: axial Gaussian peaked near the
// exit plane, with a radial component that grows toward the walls so the field
// lines bow outward. Crude, but it produces the field-strength variation along
// a wall that magnetic mirroring depends on, which is enough to exercise the
// anisotropic mobility path in the solver.
//
// REPLACE THIS with the ANSYS .fld import. It is a placeholder.
Vec<kDim> placeholderMagneticField(const Vec<kDim>& pos) {
    const double r = pos[0];
    const double z = pos[1];

    const double zPeak = 0.85 * kChannelLength;
    const double width = 0.20 * kChannelLength;
    const double envelope = std::exp(-std::pow((z - zPeak) / width, 2.0));

    const double rMid = 0.5 * kChannelWidth;
    const double radialFraction = (r - rMid) / rMid;

    Vec<kDim> B{};
    B[0] = kPeakB * envelope;                        // radial, the confining component
    B[1] = kPeakB * envelope * 0.35 * radialFraction; // axial bowing
    return B;
}

// Analytic stand-in for the imported potential, shaped after the "Potential
// Structure" panel of Fig. 13.4: flat and near zero through the upstream
// channel, then dropping steeply across the acceleration region toward the
// exit. Anode is at zero potential, per Eq. 13.67.
//
// REPLACE THIS with the ANSYS import. It is a placeholder.
double placeholderPotential(const Vec<kDim>& pos) {
    const double z = pos[1];
    const double zDrop = 0.60 * kChannelLength;
    const double width = 0.18 * kChannelLength;
    // Smooth step from 0 V down to -kDischargeVoltage.
    const double s = 0.5 * (1.0 + std::tanh((z - zDrop) / width));
    return -kDischargeVoltage * s;
}

int main() {
    // --- mesh ---------------------------------------------------------------
    Mesh<kDim> mesh({kNodesRadial, kNodesAxial},
                    {0.0, 0.0},
                    {kChannelWidth / (kNodesRadial - 1),
                     kChannelLength / (kNodesAxial - 1)});

    std::cout << "mesh: " << mesh.count(0) << " x " << mesh.count(1)
              << " = " << mesh.numNodes() << " nodes\n";
    std::cout << "spacing: " << mesh.spacing(0) << " x " << mesh.spacing(1) << " m\n";

    // --- initial state ------------------------------------------------------
    PlasmaState<kDim> state(mesh);

    for (int i = 0; i < mesh.numNodes(); ++i) {
        const auto k = static_cast<std::size_t>(i);
        const Vec<kDim> pos = mesh.position(i);
        const double z = pos[1];

        // Neutrals injected at the anode, decaying downstream.
        state.n_n[k] = kAnodeNeutralDensity * std::exp(-3.0 * z / kChannelLength);
        state.v_n[k] = Vec<kDim>{0.0, kNeutralInjectSpeed};

        // Seed plasma, weighted toward the ionization region near the exit,
        // on top of a uniform background. The background matters: without it
        // the far upstream nodes sit at the numerical density floor, the
        // Boltzmann log goes to large negative values, and the potential
        // gradient across one cell becomes enormous.
        const double zPeak = 0.7 * kChannelLength;
        const double width = 0.25 * kChannelLength;
        constexpr double kBackgroundDensity = 0.05 * kSeedPlasmaDensity;
        state.n_i[k] = kBackgroundDensity
                     + kSeedPlasmaDensity * std::exp(-std::pow((z - zPeak) / width, 2.0));
        state.n_e[k] = state.n_i[k];

        // Reverse choking at the anode, per Eq. 13.67: ions start moving
        // backward at just under the ion sonic speed.
        const double ionSound = std::sqrt(5.0 / 3.0
                                          * constants::kTe_J(kSeedElectronTemp)
                                          / constants::kXenonMass);
        state.v_i[k] = Vec<kDim>{0.0, -0.999 * ionSound};

        state.T_e[k] = kSeedElectronTemp;
        state.B[k]   = placeholderMagneticField(pos);
        state.phi[k] = placeholderPotential(pos);
    }

    // --- solver -------------------------------------------------------------
    FluidSolver<kDim>::Parameters params;
    // Imported potential, matching the plan of reading field data from ANSYS.
    // See the PotentialClosure comment in FluidSolver.hpp for why Boltzmann is
    // not a drop-in alternative here.
    params.closure          = FluidSolver<kDim>::PotentialClosure::Imported;
    params.anomalousAlpha   = 1.0 / 16.0;
    params.referenceDensity = kSeedPlasmaDensity;
    params.referencePhi     = 0.0;   // anode at zero potential, Eq. 13.67

    FluidSolver<kDim> solver(mesh, params);

    // Prime the derived quantities so the first timestep estimate sees a real
    // electron velocity rather than the zeros the state was constructed with.
    solver.step(state, 0.0);

    double dt = solver.maxStableTimestep(state);
    std::cout << "initial stable dt: " << dt << " s\n";
    std::cout << "electron plasma period at seed density: "
              << FluidSolver<kDim>::plasmaPeriod(kSeedPlasmaDensity) << " s\n";
    std::cout << "  (only binding if you replace the Boltzmann closure with Poisson)\n\n";

    // --- main loop ----------------------------------------------------------
    int frame = 0;
    for (int step = 0; step < kNumSteps; ++step) {
        // Recompute the limit periodically, but only ever shrink it.
        //
        // Letting dt grow again when velocities dip is how this loop
        // originally destabilised: the CFL number bounds advection, but the
        // energy equation source terms have their own stiffer limit that this
        // estimate does not see. Monotone dt is the conservative choice.
        if (step % 100 == 0) dt = std::min(dt, solver.maxStableTimestep(state));

        solver.step(state, dt);

        if (!state.isPhysical()) {
            std::cerr << "\nstate went non-physical at step " << step
                      << " -- reduce dt or check boundary conditions\n";
            break;
        }

        if (step % kOutputEvery == 0) {
            VtiWriter<kDim> writer(mesh);
            writer.addScalarField("electron_density",     state.n_e);
            writer.addScalarField("ion_density",          state.n_i);
            writer.addScalarField("neutral_density",      state.n_n);
            writer.addScalarField("electron_temperature", state.T_e);
            writer.addScalarField("potential",            state.phi);
            writer.addVectorField("ion_velocity",         state.v_i);
            writer.addVectorField("electron_velocity",    state.v_e);
            writer.addVectorField("neutral_velocity",     state.v_n);
            writer.addVectorField("magnetic_field",       state.B);

            const std::string filename = makeVtiFilename("field", frame);
            writer.write(filename);
            std::cout << "step " << step << " -> " << filename
                      << "   max charge separation "
                      << state.maxChargeSeparation() << '\n';
            ++frame;
        }
    }

    std::cout << "\nwrote " << frame << " frames.\n"
              << "In ParaView: File > Open, select field_..vti, click Apply,\n"
              << "then colour by electron_density and press play.\n";
    return 0;
}
