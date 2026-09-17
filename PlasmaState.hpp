#pragma once
//
// PlasmaState.hpp
//
// The seven properties of interest named at the end of Section 13.4:
//
//     v_i, v_e, v_n, n_e, n_n, phi, T_e
//
// plus n_i, which the 1-D model folds into n_e by quasineutrality but which is
// worth carrying separately in 2-D so that quasineutrality becomes something
// you can check rather than something you assume silently.
//
// Storage is structure-of-arrays: one contiguous vector per quantity, indexed
// by flat node index. This keeps each field ready to hand straight to the
// .vti writer with no repacking, and it vectorises better than an array of
// per-node structs.
//
// DCFHT Simulation Team
//

#include "Mesh.hpp"
#include "Vec.hpp"

#include <cstddef>
#include <vector>

template <int D>
struct PlasmaState {
    // number densities, m^-3
    std::vector<double> n_e;
    std::vector<double> n_i;
    std::vector<double> n_n;

    // mean velocities, m/s
    std::vector<Vec<D>> v_e;
    std::vector<Vec<D>> v_i;
    std::vector<Vec<D>> v_n;

    // electron temperature, eV  (see Constants.hpp on the eV/J convention)
    std::vector<double> T_e;

    // electric potential, V
    std::vector<double> phi;

    // imported fields (read from file, not solved for)
    std::vector<Vec<D>> B;   // tesla
    std::vector<Vec<D>> E;   // V/m, if imported rather than derived from phi

    explicit PlasmaState(const Mesh<D>& mesh) { resize(mesh.numNodes()); }

    void resize(int numNodes) {
        const auto n = static_cast<std::size_t>(numNodes);
        n_e.assign(n, 0.0);
        n_i.assign(n, 0.0);
        n_n.assign(n, 0.0);
        v_e.assign(n, Vec<D>{});
        v_i.assign(n, Vec<D>{});
        v_n.assign(n, Vec<D>{});
        T_e.assign(n, 0.0);
        phi.assign(n, 0.0);
        B.assign(n, Vec<D>{});
        E.assign(n, Vec<D>{});
    }

    std::size_t size() const { return n_e.size(); }

    // Quasineutrality residual, max over the domain of |n_i - n_e| / n_e.
    // The fluid model assumes this stays small. If it grows, either the
    // model assumption is being violated or something upstream is wrong --
    // either way you want to know, so check it rather than assume it.
    double maxChargeSeparation() const {
        double worst = 0.0;
        for (std::size_t i = 0; i < n_e.size(); ++i) {
            if (n_e[i] <= 0.0) continue;
            const double d = std::abs(n_i[i] - n_e[i]) / n_e[i];
            if (d > worst) worst = d;
        }
        return worst;
    }

    // True if every density is finite and non-negative and T_e is positive.
    // Call this each step during development; a blown-up run shows up here
    // long before it shows up in ParaView.
    bool isPhysical() const {
        for (std::size_t i = 0; i < n_e.size(); ++i) {
            if (!std::isfinite(n_e[i]) || n_e[i] < 0.0) return false;
            if (!std::isfinite(n_i[i]) || n_i[i] < 0.0) return false;
            if (!std::isfinite(n_n[i]) || n_n[i] < 0.0) return false;
            if (!std::isfinite(T_e[i]) || T_e[i] <= 0.0) return false;
            if (!std::isfinite(phi[i])) return false;
        }
        return true;
    }
};
