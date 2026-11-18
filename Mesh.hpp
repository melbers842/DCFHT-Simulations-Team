#pragma once
//
// Mesh.hpp
//
// Uniform structured mesh covering the DCFHT channel, templated on dimension.
//
// AXES: axis 0 is r (radial, from the thruster axis), axis 1 is z (axial, from
// the anode face). In 3-D, axis 2 is azimuthal.
//
// The mesh is a rectangle in r-z that BOUNDS the cone. The cone itself is
// carved out of it by masking: isInsideChannel() tells you which nodes are
// plasma and which are outside the wall. Cells outside are still allocated --
// they just never get updated. That costs some memory and keeps every
// difference operator working on a regular grid, which is the whole point of
// the masked approach.
//
// COORDINATE SYSTEM: see the Geometry enum below. This is a real modelling
// decision, not a switch to flip casually.
//
// DCFHT Simulation Team
//

#include "ChannelGeometry.h"
#include "Vec.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

// How the divergence operator treats the radial direction.
//
//   CartesianSlab -- div(F) = dF_r/dr + dF_z/dz. Simple, and WRONG for a cone:
//       it ignores that a shell at larger r has more area, so radial fluxes do
//       not conserve mass. Kept so the two can be compared and the error
//       measured rather than assumed.
//
//   Axisymmetric  -- div(F) = (1/r) d(r F_r)/dr + dF_z/dz. Correct for a solid
//       of revolution. The 1/r is singular on the axis; since the cone
//       INCLUDES the axis, the r = 0 row needs the symmetry treatment in
//       divergence() below. An annular channel would not have this problem;
//       a cone does.
enum class Geometry { CartesianSlab, Axisymmetric };

template <int D>
class Mesh {
public:
    Mesh(std::array<int, D> nodeCounts,
         std::array<double, D> origin,
         std::array<double, D> spacing,
         Geometry coordinateSystem = Geometry::Axisymmetric,
         std::optional<ChannelGeometry> channel = std::nullopt)
        : n_(nodeCounts), x0_(origin), dx_(spacing),
          coords_(coordinateSystem), channel_(channel) {
        total_ = 1;
        for (int d = 0; d < D; ++d) {
            assert(n_[static_cast<std::size_t>(d)] >= 2 && "each axis needs at least 2 nodes");
            assert(dx_[static_cast<std::size_t>(d)] > 0.0 && "spacing must be positive");
            total_ *= n_[static_cast<std::size_t>(d)];
        }
        if (channel_) channel_->validateOrThrow();
    }

    // Build a mesh sized to a channel: radial extent covers the widest part of
    // the cone, axial extent covers its length. Takes the node counts and works
    // out origin and spacing, so callers never compute them by hand.
    static Mesh<D> forChannel(const ChannelGeometry& geom,
                              int radialNodes, int axialNodes,
                              Geometry coordinateSystem = Geometry::Axisymmetric) {
        std::array<int, D> counts{};
        std::array<double, D> origin{};
        std::array<double, D> spacing{};
        counts[0] = radialNodes;
        counts[1] = axialNodes;
        origin[0] = 0.0;              // the axis itself is a mesh node
        origin[1] = 0.0;              // anode face
        spacing[0] = geom.exitRadius() / (radialNodes - 1);
        spacing[1] = geom.length      / (axialNodes  - 1);
        for (int d = 2; d < D; ++d) { counts[d] = 2; origin[d] = 0.0; spacing[d] = 1.0; }
        return Mesh<D>(counts, origin, spacing, coordinateSystem, geom);
    }

    int    numNodes()     const { return total_; }
    int    count(int d)   const { return n_[static_cast<std::size_t>(d)]; }
    double spacing(int d) const { return dx_[static_cast<std::size_t>(d)]; }
    double origin(int d)  const { return x0_[static_cast<std::size_t>(d)]; }
    Geometry coordinateSystem() const { return coords_; }
    const std::optional<ChannelGeometry>& channel() const { return channel_; }

    // --- index <-> multi-index, x-fastest to match the VTK layout -----------
    int flatten(const std::array<int, D>& idx) const {
        int flat = 0, stride = 1;
        for (int d = 0; d < D; ++d) {
            flat += idx[static_cast<std::size_t>(d)] * stride;
            stride *= n_[static_cast<std::size_t>(d)];
        }
        return flat;
    }

    std::array<int, D> unflatten(int flat) const {
        std::array<int, D> idx{};
        for (int d = 0; d < D; ++d) {
            idx[static_cast<std::size_t>(d)] = flat % n_[static_cast<std::size_t>(d)];
            flat /= n_[static_cast<std::size_t>(d)];
        }
        return idx;
    }

    Vec<D> position(int flat) const {
        const auto idx = unflatten(flat);
        Vec<D> p;
        for (int d = 0; d < D; ++d) {
            p[d] = x0_[static_cast<std::size_t>(d)]
                 + idx[static_cast<std::size_t>(d)] * dx_[static_cast<std::size_t>(d)];
        }
        return p;
    }

    double radiusAt(int flat) const { return position(flat)[0]; }
    double axialAt(int flat)  const { return position(flat)[1]; }

    int neighbour(int flat, int d, int offset) const {
        auto idx = unflatten(flat);
        const int moved = idx[static_cast<std::size_t>(d)] + offset;
        if (moved < 0 || moved >= n_[static_cast<std::size_t>(d)]) return -1;
        idx[static_cast<std::size_t>(d)] = moved;
        return flatten(idx);
    }

    bool isOnAxis(int flat) const {
        return unflatten(flat)[0] == 0 && x0_[0] == 0.0;
    }

    // --- channel masking ----------------------------------------------------
    //
    // A node is plasma when it lies at or inside the cone wall. With no channel
    // attached every node counts as inside, so existing rectangular-domain code
    // keeps working unchanged.
    //
    // The comparison is <=, so a node sitting exactly on the wall is INSIDE.
    // Those are the nodes wall boundary conditions act on, and they are easier
    // to find if they are part of the domain rather than excluded from it.
    bool isInsideChannel(int flat) const {
        if (!channel_) return true;
        const Vec<D> p = position(flat);
        const double r = p[0];
        const double z = p[1];
        if (r < 0.0) return false;
        if (z < 0.0 || z > channel_->length) return false;
        return r <= channel_->wallRadiusAt(z);
    }

    // Inside the channel, but with at least one neighbour outside it. These are
    // the nodes the wall boundary conditions apply to.
    bool isWallAdjacent(int flat) const {
        if (!channel_ || !isInsideChannel(flat)) return false;
        for (int d = 0; d < D; ++d) {
            for (int off : {-1, +1}) {
                const int nb = neighbour(flat, d, off);
                if (nb < 0 || !isInsideChannel(nb)) return true;
            }
        }
        return false;
    }

    // Outward unit normal of the cone wall, pointing from plasma into the wall.
    // The wall is the line r = anodeRadius + z*tan(halfAngle), so a tangent
    // along it is (tan, 1)/|..| and the outward normal is (1, -tan)/|..|.
    // Independent of position because the wall is straight in r-z.
    Vec<D> wallNormal() const {
        Vec<D> nrm;
        if (!channel_) return nrm;
        const double t = std::tan(channel_->halfAngleRad());
        const double mag = std::sqrt(1.0 + t * t);
        nrm[0] = 1.0 / mag;
        nrm[1] = -t / mag;
        return nrm;
    }

    int countInside() const {
        int c = 0;
        for (int i = 0; i < total_; ++i) if (isInsideChannel(i)) ++c;
        return c;
    }

    // --- difference operators -----------------------------------------------
    // Central in the interior, one-sided at the edges. Second order inside,
    // first order on the boundary.

    Vec<D> gradient(const std::vector<double>& field, int flat) const {
        Vec<D> g;
        for (int d = 0; d < D; ++d) {
            const int lo = neighbour(flat, d, -1);
            const int hi = neighbour(flat, d, +1);
            const double h = dx_[static_cast<std::size_t>(d)];
            const auto f = [&](int i) { return field[static_cast<std::size_t>(i)]; };
            if (lo >= 0 && hi >= 0)      g[d] = (f(hi) - f(lo)) / (2.0 * h);
            else if (hi >= 0)            g[d] = (f(hi) - f(flat)) / h;
            else if (lo >= 0)            g[d] = (f(flat) - f(lo)) / h;
            else                         g[d] = 0.0;
        }
        return g;
    }

    // Divergence of a flux (scalar * velocity), the form every continuity
    // equation uses. Built on the product rather than scalar*div(v), because
    // the product form conserves mass much better.
    //
    // In Axisymmetric mode the radial term carries the r weighting:
    //     (1/r) d(r F_r)/dr
    // On the axis that is 0/0. By symmetry F_r = 0 there, and the limit of the
    // radial term is 2 dF_r/dr, which is what the axis branch below uses.
    double divergenceOfFlux(const std::vector<double>& scalar,
                            const std::vector<Vec<D>>& velocity,
                            int flat) const {
        double div = 0.0;
        for (int d = 0; d < D; ++d) {
            const int lo = neighbour(flat, d, -1);
            const int hi = neighbour(flat, d, +1);
            const double h = dx_[static_cast<std::size_t>(d)];
            const auto flux = [&](int i) {
                return scalar[static_cast<std::size_t>(i)]
                     * velocity[static_cast<std::size_t>(i)][d];
            };

            const bool radial = (d == 0 && coords_ == Geometry::Axisymmetric);
            if (!radial) {
                if (lo >= 0 && hi >= 0) div += (flux(hi) - flux(lo)) / (2.0 * h);
                else if (hi >= 0)       div += (flux(hi) - flux(flat)) / h;
                else if (lo >= 0)       div += (flux(flat) - flux(lo)) / h;
                continue;
            }

            const double r = radiusAt(flat);
            if (r <= 0.0) {
                // On the axis: lim (1/r) d(r F)/dr = 2 dF/dr, with F_r = 0 by
                // symmetry. Use the one-sided difference outward.
                if (hi >= 0) div += 2.0 * (flux(hi) - flux(flat)) / h;
                continue;
            }
            const auto rflux = [&](int i) { return radiusAt(i) * flux(i); };
            if (lo >= 0 && hi >= 0) div += (rflux(hi) - rflux(lo)) / (2.0 * h * r);
            else if (hi >= 0)       div += (rflux(hi) - rflux(flat)) / (h * r);
            else if (lo >= 0)       div += (rflux(flat) - rflux(lo)) / (h * r);
        }
        return div;
    }

    double minSpacing() const {
        double m = dx_[0];
        for (int d = 1; d < D; ++d)
            if (dx_[static_cast<std::size_t>(d)] < m) m = dx_[static_cast<std::size_t>(d)];
        return m;
    }

private:
    std::array<int, D>    n_;
    std::array<double, D> x0_;
    std::array<double, D> dx_;
    Geometry              coords_;
    std::optional<ChannelGeometry> channel_;
    int                   total_ = 0;
};
