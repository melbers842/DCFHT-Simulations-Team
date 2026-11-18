#pragma once
//
// ChannelGeometry.hpp
//
// Geometry of the DCFHT discharge channel: a cone of revolution about the
// thruster axis, opening downstream.
//
// GEOMETRY, stated once so nobody has to re-derive it from the code:
//
//   - The channel is a CONE, not an annulus. There is no centrebody. Plasma
//     fills the interior from the axis out to the wall.
//   - r is distance from the thruster axis; z is axial distance from the
//     anode face. The simulated domain is the r-z half-plane, and the wall is
//     a surface of revolution about r = 0.
//   - A point is inside the channel when 0 <= r <= wallRadiusAt(z).
//   - divergenceAngleDeg is the FULL opening angle of the cone. The wall
//     itself tilts by half that. Getting this wrong is a silent factor of two,
//     so halfAngleRad() is the only place the conversion happens.
//
// Dimensions come from the team's CATIA half-section of the thruster.
//
// DCFHT Simulation Team
//

#include "Constants.hpp"

#include <cmath>
#include <stdexcept>

struct ChannelGeometry {
    double anodeRadius;          // m, wall radius at z = 0
    double length;               // m, anode face to exit plane
    double divergenceAngleDeg;   // FULL opening angle of the cone

    // Largest radius the field data covers. The ANSYS export spans
    // x in [-60, 60] mm, so the wall must stay inside that or the import
    // has nothing to interpolate from at the exit plane.
    static constexpr double kMaxFieldRadius = 0.060;  // m

    double halfAngleRad() const {
        return 0.5 * divergenceAngleDeg * (constants::kPi / 180.0);
    }

    // The single wall. Radius grows with z because the cone opens downstream.
    double wallRadiusAt(double z) const {
        return anodeRadius + z * std::tan(halfAngleRad());
    }

    double exitRadius() const { return wallRadiusAt(length); }

    // Flow area at an axial station, for mass-flux bookkeeping later.
    double areaAt(double z) const {
        const double r = wallRadiusAt(z);
        return constants::kPi * r * r;
    }

    // Build a geometry from two measured radii instead of an angle, which is
    // how the CATIA drawing gives it: radius at the anode, radius at the exit.
    static ChannelGeometry fromRadii(double anodeR, double exitR, double len) {
        if (len <= 0.0) throw std::invalid_argument("ChannelGeometry: length must be positive");
        const double halfAngle = std::atan((exitR - anodeR) / len);
        return ChannelGeometry{anodeR, len,
                               2.0 * halfAngle * (180.0 / constants::kPi)};
    }

    // Largest FULL opening angle whose exit radius still fits the field data.
    double maxValidAngleDeg() const {
        const double usable = kMaxFieldRadius - anodeRadius;
        if (usable <= 0.0) return 0.0;
        return 2.0 * std::atan(usable / length) * (180.0 / constants::kPi);
    }

    bool isValid() const {
        return anodeRadius > 0.0
            && length > 0.0
            && divergenceAngleDeg > 0.0
            && divergenceAngleDeg < 180.0
            && exitRadius() <= kMaxFieldRadius;
    }

    void validateOrThrow() const {
        if (!isValid()) {
            throw std::invalid_argument(
                "ChannelGeometry: invalid channel -- check that the radius and "
                "length are positive and that the exit radius stays within the "
                "field data extent");
        }
    }
};

// Dimensions from the team's CATIA half-section: 10 mm wall radius at the
// anode, 38.579 mm at the exit, over a 49.5 mm channel. Those three imply a
// 30 degree half-angle, so 60 degrees full opening angle.
inline const ChannelGeometry kDcfhtGeometry{0.010, 0.0495, 60.0};
