//
// Created by Angel Madrigal on 10/21/26.
//

#ifndef URIMESHPRACTICE_CHANNELGEOMETRY_H
#define URIMESHPRACTICE_CHANNELGEOMETRY_H

#include <cmath>
#include <stdexcept>
#include "Constants.hpp"


struct ChannelGeometry {
    double innerRadius;             //m, at z = 0 
    double outerRadius;             //m, at z = 0
    double length;                  //m
    double totaldivergenceAngleDeg;      //half-angle, 30-75

    // Minimum clearance the inboard wall must keep from the thruster axis.
    // Not physics -- a guard so an impossible geometry fails loudly instead of
    // producing a mesh with negative radii.
    static constexpr double kMinAxisClearance = 0.002;  // m

    double halfangleRad() const {
        return 0.5 * totaldivergenceAngleDeg * ( constants::kPi/180.0);
    }

    // Inboard wall tilts toward the axis, outboard wall tilts away, so the
    // channel opens into a V with increasing z. Both tilt by the same angle.
    double innerWallAt(double z) const {
        return innerRadius - z*std::tan(halfangleRad());
    }     //innerRadius - z*tan
    double outerWallAt(double z) const {
        return outerRadius + z*std::tan(halfangleRad());
    }     //outerRadius + z*tan

    double channelWidthAt(double z) const {
        return outerWallAt(z) - innerWallAt(z);
    }

    // Largest divergence angle this channel can take before the inboard wall
    // reaches the axis clearance limit at the exit plane.
    //
    //   innerRadius - length*tan(theta) >= kMinAxisClearance
    //   tan(theta) <= (innerRadius - kMinAxisClearance) / length
    double maxValidAngleDeg() const {
        const double usable = innerRadius - kMinAxisClearance;
        if (usable <= 0.0) return 0.0;
        return std::atan(usable / length) * (180.0 / constants::kPi) * 2;
    }

    bool isValid() const {
        return innerRadius > 0.0
            && outerRadius > innerRadius
            && length > 0.0
            && totaldivergenceAngleDeg >= 0.0
            && innerWallAt(length) >= kMinAxisClearance;
    }

    void validateOrThrow() const {
        if (!isValid()) {
            throw std::invalid_argument(
                "ChannelGeometry: invalid channel -- inboard wall reaches the axis "
                "before the exit plane, or radii/length are nonpositive");
        }
    }

};

// PLACEHOLDER dimensions, roughly BHT-200 scale. NOT measured, NOT confirmed.
// Replace once Dr. White provides the real DCFHT geometry.
// inline constexpr ChannelGeometry kPlaceholderGeometry{
//     0.021,   // innerRadius, m
//     0.028,   // outerRadius, m
//     0.025,   // length, m
//     30.0     // divergenceAngleDeg
// };




#endif //URIMESHPRACTICE_CHANNELGEOMETRY_H
