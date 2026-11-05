//
// Created by Angel Madrigal on 11/4/26.
//

#include "ChannelGeometry.h"
#include <iostream>

int main() {

    ChannelGeometry placeholderGeometry{
        0.021,   // innerRadius, m
        0.028,   // outerRadius, m
        0.025,   // length, m
        30.0     // divergenceAngleDeg
    };

    std::cout << "max valid total angle: " << placeholderGeometry.maxValidAngleDeg() << " deg" << std::endl;

    std::cout << "inner at z = 0: " << placeholderGeometry.innerWallAt(0.0) * 1000<< "mm" << std::endl;
    std::cout << "outer at z = 0: " <<  placeholderGeometry.outerWallAt(0.0) * 1000 << "mm" << std::endl;
    std::cout << "channel width at z = 0: " <<  placeholderGeometry.channelWidthAt(0.0) * 1000 << "mm" << std::endl;

    std::cout << "inner at z = length/2: " << placeholderGeometry.innerWallAt(placeholderGeometry.length/2) * 1000 << "mm" << std::endl;
    std::cout << "outer at z = length/2: " << placeholderGeometry.outerWallAt(placeholderGeometry.length/2) * 1000 << "mm" << std::endl;
    std::cout << "channel width at z = length/2: " <<  placeholderGeometry.channelWidthAt(placeholderGeometry.length/2) * 1000 << "mm" << std::endl;

    std::cout << "inner at z = length: " << placeholderGeometry.innerWallAt(placeholderGeometry.length) * 1000 << "mm" << std::endl;
    std::cout << "outer at z = length: " << placeholderGeometry.outerWallAt(placeholderGeometry.length) * 1000 << "mm" << std::endl;
    std::cout << "channel width at z = length: " <<  placeholderGeometry.channelWidthAt(placeholderGeometry.length) * 1000 << "mm" << std::endl;

    placeholderGeometry.validateOrThrow();

    return 0;
}