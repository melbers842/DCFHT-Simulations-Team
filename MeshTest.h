//
// Created by Angel Madrigal on 1/12/26.
//

#ifndef MESHTEST_H
#define MESHTEST_H

#include <Eigen/Dense>
#include <string>
//using Eigen::MatrixXd; or using namespace Eigen;

// This class represents a 2D (later 3D) mesh for the Hall thruster simulation
class Mesh {
private:
    // Grid dimensions
    int xN, yN;          // Number of nodes in X and Y
    double xD, yD;       // Grid spacing
    double minX, minY;   // Grid origin coordinates

    // Physical properties stored at each mesh node
    // Each Eigen::MatrixXd is Nx x Ny (2D)
    // Later can be extended to Nx x Ny x Nz for 3D

    Eigen::MatrixXd electronVelocityX;  // Electron velocity in X
    Eigen::MatrixXd electronVelocityY;  // Electron velocity in Y
    Eigen::MatrixXd ionVelocityX;       // Ion velocity in X
    Eigen::MatrixXd ionVelocityY;       // Ion velocity in Y
    Eigen::MatrixXd neutralVelocityX;   // Neutral velocity in X
    Eigen::MatrixXd neutralVelocityY;   // Neutral velocity in Y

    Eigen::MatrixXd electronDensity;    // Number density of electrons
    Eigen::MatrixXd ionDensity;         // Number density of ions
    Eigen::MatrixXd neutralDensity;     // Number density of neutrals

    Eigen::MatrixXd electronTemperature; // Electron temperature
    Eigen::MatrixXd voltagePotential;    // Electric potential at nodes

public:

    // Constructor: initializes the mesh and allocates storage
    Mesh(int nx, int ny, double dx, double dy, double xmin, double ymin);

    // Initialize all fields to default values
    // TODO: later replace with real initial conditions
    void static initializeFields();

    // Print mesh info for debugging
    void printInfo() const;

    // Placeholder: update velocities of particles
    // TODO: will be implemented in solver
    void static updateVelocities();

    // Placeholder: update densities
    // TODO: will be implemented in solver
    void static updateDensities();

    // Placeholder: update electric potential
    // TODO: solve Poisson's equation in solver
    void static updatePotential();

    // Placeholder: update electron temperature
    void static updateElectronTemperature();

    // Optional: save mesh state to file for later visualization
    // TODO: implement in io.hpp/io.cpp
    void saveMeshToFile(const std::string& filename) const;
};

#endif //MESHTEST_H