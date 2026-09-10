//
// vtu_demo_main.cpp
//
// Proves the output pipeline works, with no physics involved.
// Writes 50 timesteps of a rotating ring of 100 particles, each carrying a
// scalar ("angle") and a vector ("velocity").
//
// Open all 50 files at once in ParaView and press play.
//
// DCFHT Simulation Team
//

#include "vtu_writer.hpp"

#include <cmath>
#include <iostream>

int main() {
    constexpr int    kNumParticles = 100;
    constexpr int    kNumSteps     = 50;
    constexpr double kRadius       = 0.05;   // metres
    constexpr double kPi           = 3.14159265358979323846;

    for (int step = 0; step < kNumSteps; ++step) {
        // Rotate the whole ring a little further each timestep.
        const double spin = 2.0 * kPi * step / kNumSteps;

        VtuWriter writer;
        std::vector<double> angles;
        std::vector<double> velocities;

        for (int i = 0; i < kNumParticles; ++i) {
            const double theta = 2.0 * kPi * i / kNumParticles + spin;

            const double x = kRadius * std::cos(theta);
            const double y = kRadius * std::sin(theta);
            writer.addPoint(x, y, 0.0);

            // One scalar per particle: its angle around the ring.
            angles.push_back(theta);

            // One vector per particle: tangential velocity.
            velocities.push_back(-std::sin(theta));
            velocities.push_back( std::cos(theta));
            velocities.push_back( 0.0);
        }

        writer.addScalarField("angle", angles);
        writer.addVectorField("velocity", velocities);

        const std::string filename = makeFilename("out", step);
        writer.write(filename);
        std::cout << "wrote " << filename << '\n';
    }

    std::cout << "\nDone. In ParaView: File > Open, select all out_*.vtu,\n"
              << "click Apply, then set the coloring dropdown to \"angle\"\n"
              << "and press the play button.\n";
    return 0;
}
