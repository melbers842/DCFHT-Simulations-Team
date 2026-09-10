#pragma once
//
// vtu_writer.hpp
//
// Minimal, dependency-free writer for VTK XML unstructured grid files (.vtu)
// containing a point cloud (one VTK_VERTEX cell per point).
//
// Intended for writing particle data that ParaView can open directly.
// Nothing here knows anything about plasma physics -- it just takes numbers
// and arranges them in the layout ParaView expects.
//
// Usage:
//     VtuWriter w;
//     w.addPoint(x, y, z);                 // repeat for every particle
//     w.addScalarField("energy", energies);
//     w.addVectorField("velocity", vels);  // vels.size() == 3 * numPoints
//     w.write("out_0000.vtu");
//
// DCFHT Simulation Team
//

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

class VtuWriter {
public:
    // Add one particle position. Use z = 0.0 for a 2-D simulation;
    // ParaView always works in 3-D internally.
    void addPoint(double x, double y, double z = 0.0) {
        points_.push_back(x);
        points_.push_back(y);
        points_.push_back(z);
    }

    std::size_t numPoints() const { return points_.size() / 3; }

    // Attach one number per particle (temperature, energy, charge, ...).
    void addScalarField(const std::string& name, const std::vector<double>& values) {
        if (values.size() != numPoints()) {
            throw std::runtime_error(
                "VtuWriter: scalar field '" + name + "' has " +
                std::to_string(values.size()) + " values but there are " +
                std::to_string(numPoints()) + " points");
        }
        scalars_.push_back({name, values});
    }

    // Attach three numbers per particle (velocity, force, ...).
    // Layout is flat and interleaved: x0,y0,z0, x1,y1,z1, ...
    void addVectorField(const std::string& name, const std::vector<double>& values) {
        if (values.size() != 3 * numPoints()) {
            throw std::runtime_error(
                "VtuWriter: vector field '" + name + "' has " +
                std::to_string(values.size()) + " values but expected " +
                std::to_string(3 * numPoints()));
        }
        vectors_.push_back({name, values});
    }

    // Remove all points and fields so the writer can be reused next timestep.
    void clear() {
        points_.clear();
        scalars_.clear();
        vectors_.clear();
    }

    void write(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out) {
            throw std::runtime_error("VtuWriter: could not open '" + filename + "' for writing");
        }
        out << std::setprecision(9);

        const std::size_t n = numPoints();

        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "  <UnstructuredGrid>\n";
        out << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfCells=\"" << n << "\">\n";

        // --- particle positions ---------------------------------------------
        out << "      <Points>\n";
        out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (std::size_t i = 0; i < n; ++i) {
            out << "          " << points_[3 * i] << ' '
                << points_[3 * i + 1] << ' '
                << points_[3 * i + 2] << '\n';
        }
        out << "        </DataArray>\n";
        out << "      </Points>\n";

        // --- per-particle data ----------------------------------------------
        out << "      <PointData";
        if (!scalars_.empty()) out << " Scalars=\"" << scalars_.front().name << "\"";
        if (!vectors_.empty()) out << " Vectors=\"" << vectors_.front().name << "\"";
        out << ">\n";

        for (const auto& f : scalars_) {
            out << "        <DataArray type=\"Float64\" Name=\"" << f.name
                << "\" format=\"ascii\">\n          ";
            for (std::size_t i = 0; i < n; ++i) {
                out << f.values[i] << ' ';
            }
            out << "\n        </DataArray>\n";
        }

        for (const auto& f : vectors_) {
            out << "        <DataArray type=\"Float64\" Name=\"" << f.name
                << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (std::size_t i = 0; i < n; ++i) {
                out << "          " << f.values[3 * i] << ' '
                    << f.values[3 * i + 1] << ' '
                    << f.values[3 * i + 2] << '\n';
            }
            out << "        </DataArray>\n";
        }
        out << "      </PointData>\n";

        // --- cells ----------------------------------------------------------
        // Every point is its own VTK_VERTEX cell (type 1). ParaView needs each
        // point to belong to a cell before it will draw it.
        out << "      <Cells>\n";

        out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n          ";
        for (std::size_t i = 0; i < n; ++i) out << i << ' ';
        out << "\n        </DataArray>\n";

        out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
        for (std::size_t i = 0; i < n; ++i) out << (i + 1) << ' ';
        out << "\n        </DataArray>\n";

        out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
        for (std::size_t i = 0; i < n; ++i) out << "1 ";
        out << "\n        </DataArray>\n";

        out << "      </Cells>\n";
        out << "    </Piece>\n";
        out << "  </UnstructuredGrid>\n";
        out << "</VTKFile>\n";
    }

private:
    struct Field {
        std::string name;
        std::vector<double> values;
    };

    std::vector<double> points_;
    std::vector<Field> scalars_;
    std::vector<Field> vectors_;
};

// Build a zero-padded time-series filename, e.g. makeFilename("out", 42) -> "out_0042.vtu".
// ParaView groups files named this way into an animation automatically.
inline std::string makeFilename(const std::string& prefix, int step, int width = 4) {
    std::string s = std::to_string(step);
    while (static_cast<int>(s.size()) < width) s.insert(s.begin(), '0');
    return prefix + "_" + s + ".vtu";
}
