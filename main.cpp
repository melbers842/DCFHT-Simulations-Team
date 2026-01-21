// Gabriella Garner 01.20.26 with ChatGPT assistance

#include "mesh.hpp"
#include <iostream>

Mesh::Mesh(int nx, int ny, double dx_, double dy_, double xmin_, double ymin_)
    : Nx(nx), Ny(ny), dx(dx_), dy(dy_), xmin(xmin_), ymin(ymin_)
{
    // Allocate matrices for all fields
    electronVelocityX = Eigen::MatrixXd::Zero(Nx, Ny);
    electronVelocityY = Eigen::MatrixXd::Zero(Nx, Ny);
    ionVelocityX = Eigen::MatrixXd::Zero(Nx, Ny);
    ionVelocityY = Eigen::MatrixXd::Zero(Nx, Ny);
    neutralVelocityX = Eigen::MatrixXd::Zero(Nx, Ny);
    neutralVelocityY = Eigen::MatrixXd::Zero(Nx, Ny);

    electronDensity = Eigen::MatrixXd::Zero(Nx, Ny);
    ionDensity = Eigen::MatrixXd::Zero(Nx, Ny);
    neutralDensity = Eigen::MatrixXd::Zero(Nx, Ny);

    electronTemperature = Eigen::MatrixXd::Zero(Nx, Ny);
    voltagePotential = Eigen::MatrixXd::Zero(Nx, Ny);
}

void Mesh::initializeFields() {
    // TODO: set initial conditions for velocities, densities, potentials, temperatures
    // e.g., electrons start with small random velocity, uniform density, 0 potential
    std::cout << "Initializing mesh fields (placeholders)\n";
}

void Mesh::printInfo() const {
    std::cout << "Mesh info:\n";
    std::cout << "Nx = " << Nx << ", Ny = " << Ny << ", dx = " << dx << ", dy = " << dy << "\n";
    std::cout << "Sample values (first node):\n";
    std::cout << "Electron velocity X: " << electronVelocityX(0,0) << "\n";
    std::cout << "Electron density: " << electronDensity(0,0) << "\n";
    std::cout << "Voltage potential: " << voltagePotential(0,0) << "\n";
}

void Mesh::updateVelocities() {
    // TODO: implement physics here
    // e.g., particle motion based on E and B fields
}

void Mesh::updateDensities() {
    // TODO: implement particle source, sink, or collisions
}

void Mesh::updatePotential() {
    // TODO: solve Poisson equation for electric potential
}

void Mesh::updateElectronTemperature() {
    // TODO: compute electron energy balance
}

void Mesh::saveMeshToFile(const std::string& filename) const {
    // TODO: write mesh fields to file (CSV or VTK)
    std::cout << "Saving mesh to file (placeholder): " << filename << "\n";
}