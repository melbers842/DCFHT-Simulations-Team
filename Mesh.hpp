#pragma once
//
// Mesh.hpp
//
// Uniform structured Cartesian mesh, templated on dimension.
//
// For an axisymmetric DCFHT the natural 2-D interpretation of the axes is
// (r, z): axis 0 = radial, axis 1 = axial. The mesh itself is agnostic --
// it just stores node counts, an origin and a spacing per axis. Attaching
// physical meaning to the axes is the solver's job.
//
// Node ordering is x-fastest (axis 0 varies fastest), which matches the
// VTK convention, so mesh data can be written straight out to a .vti file
// without reordering.
//
// DCFHT Simulation Team
//

#include "Vec.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

template <int D>
class Mesh {

private:
    std::array<int, D>    nodeCount_;
    std::array<double, D> x0_; //X initial
    std::array<double, D> dx_; //spacing
    int                   total_ = 0;

public:
    Mesh(std::array<int, D> nodeCounts,
         std::array<double, D> origin,
         std::array<double, D> spacing)
        : nodeCount_(nodeCounts), x0_(origin), dx_(spacing) {
        total_ = 1;
        //
        for (int d = 0; d < D; ++d) {
            assert(n_[static_cast<std::size_t>(d)] >= 2 && "each axis needs at least 2 nodes");
            assert(dx_[static_cast<std::size_t>(d)] > 0.0 && "spacing must be positive");
            total_ *= nodeCount_[static_cast<std::size_t>(d)];
        }
    }

    int  numNodes()          const { return total_; }
    int  count(int d)        const { return nodeCount_[static_cast<std::size_t>(d)]; }
    double spacing(int d)    const { return dx_[static_cast<std::size_t>(d)]; }
    double origin(int d)     const { return x0_[static_cast<std::size_t>(d)]; }

    // --- index <-> multi-index ---------------------------------------------
    //Flatten from 3D to 2D
    int flatten(const std::array<int, D>& idx) const {
        int flat = 0;
        int stride = 1;
        for (int d = 0; d < D; ++d) {
            flat += idx[static_cast<std::size_t>(d)] * stride;
            stride *= nodeCount_[static_cast<std::size_t>(d)];
        }
        return flat;
    }

    std::array<int, D> unflatten(int flat) const {
        std::array<int, D> idx{};
        for (int d = 0; d < D; ++d) {
            idx[static_cast<std::size_t>(d)] = flat % nodeCount_[static_cast<std::size_t>(d)];
            flat /= nodeCount_[static_cast<std::size_t>(d)];
        }
        return idx;
    }

    // Physical position of a node.
    Vec<D> position(int flat) const {
        const auto idx = unflatten(flat);
        Vec<D> p;
        for (int d = 0; d < D; ++d) {
            p[d] = x0_[static_cast<std::size_t>(d)]
                 + idx[static_cast<std::size_t>(d)] * dx_[static_cast<std::size_t>(d)];
        }
        return p;
    }

    bool isBoundary(int flat) const {
        const auto idx = unflatten(flat);
        for (int d = 0; d < D; ++d) {
            if (idx[static_cast<std::size_t>(d)] == 0 ||
                idx[static_cast<std::size_t>(d)] == nodeCount_[static_cast<std::size_t>(d)] - 1) {
                return true;
            }
        }
        return false;
    }

    // Neighbour index along axis d, offset by +1 or -1. Returns -1 if outside.
    int neighbour(int flat, int d, int offset) const {
        auto idx = unflatten(flat);
        const int moved = idx[static_cast<std::size_t>(d)] + offset;
        if (moved < 0 || moved >= nodeCount_[static_cast<std::size_t>(d)]) return -1;
        idx[static_cast<std::size_t>(d)] = moved;
        return flatten(idx);
    }

    // --- finite-difference operators ---------------------------------------
    //
    // Central differences in the interior, one-sided at the boundaries.
    // Second-order accurate in the interior, first-order on the boundary --
    // adequate to get running, worth upgrading later if boundary accuracy
    // turns out to matter for the mirroring results.

    Vec<D> gradient(const std::vector<double>& field, int flat) const {
        Vec<D> g;
        for (int d = 0; d < D; ++d) {
            const int lo = neighbour(flat, d, -1);
            const int hi = neighbour(flat, d, +1);
            const double h = dx_[static_cast<std::size_t>(d)];
            if (lo >= 0 && hi >= 0) {
                g[d] = (field[static_cast<std::size_t>(hi)] -
                        field[static_cast<std::size_t>(lo)]) / (2.0 * h);
            } else if (hi >= 0) {
                g[d] = (field[static_cast<std::size_t>(hi)] -
                        field[static_cast<std::size_t>(flat)]) / h;
            } else if (lo >= 0) {
                g[d] = (field[static_cast<std::size_t>(flat)] -
                        field[static_cast<std::size_t>(lo)]) / h;
            } else {
                g[d] = 0.0;
            }
        }
        return g;
    }

    // Divergence of a vector field stored per node.
    double divergence(const std::vector<Vec<D>>& field, int flat) const {
        double div = 0.0;
        for (int d = 0; d < D; ++d) {
            const int lo = neighbour(flat, d, -1);
            const int hi = neighbour(flat, d, +1);
            const double h = dx_[static_cast<std::size_t>(d)];
            if (lo >= 0 && hi >= 0) {
                div += (field[static_cast<std::size_t>(hi)][d] -
                        field[static_cast<std::size_t>(lo)][d]) / (2.0 * h);
            } else if (hi >= 0) {
                div += (field[static_cast<std::size_t>(hi)][d] -
                        field[static_cast<std::size_t>(flat)][d]) / h;
            } else if (lo >= 0) {
                div += (field[static_cast<std::size_t>(flat)][d] -
                        field[static_cast<std::size_t>(lo)][d]) / h;
            }
        }
        return div;
    }

    // Divergence of the flux (scalar * vector), which is the form that appears
    // in every continuity equation. Computed on the product rather than as
    // scalar * div(vector), because the product form conserves mass far better.
    double divergenceOfFlux(const std::vector<double>& scalar,
                            const std::vector<Vec<D>>& velocity,
                            int flat) const {
        double div = 0.0;
        for (int d = 0; d < D; ++d) {
            const int lo = neighbour(flat, d, -1);
            const int hi = neighbour(flat, d, +1);
            const double h = dx_[static_cast<std::size_t>(d)];
            auto flux = [&](int i) {
                return scalar[static_cast<std::size_t>(i)] *
                       velocity[static_cast<std::size_t>(i)][d];
            };
            if (lo >= 0 && hi >= 0) {
                div += (flux(hi) - flux(lo)) / (2.0 * h);
            } else if (hi >= 0) {
                div += (flux(hi) - flux(flat)) / h;
            } else if (lo >= 0) {
                div += (flux(flat) - flux(lo)) / h;
            }
        }
        return div;
    }

    // Smallest cell dimension, used for CFL-style timestep limits.
    double minSpacing() const {
        double m = dx_[0];
        for (int d = 1; d < D; ++d) m = (dx_[static_cast<std::size_t>(d)] < m)
                                        ? dx_[static_cast<std::size_t>(d)] : m;
        return m;
    }
};
