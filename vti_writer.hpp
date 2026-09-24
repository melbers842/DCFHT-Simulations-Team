#pragma once
//
// vti_writer.hpp
//
// Writes VTK XML image data (.vti) for a uniform structured mesh, the mesh
// counterpart to vtu_writer.hpp. ParaView opens a numbered series of these as
// an animation exactly the same way.
//
// Use .vti for node quantities (densities, temperature, potential, velocity
// fields). Use .vtu for scattered particles, if you later add a kinetic
// species.
//
// Node ordering must be x-fastest, which is what Mesh<D> already uses, so
// state vectors can be handed over with no repacking.
//
// DCFHT Simulation Team
//

#include "Mesh.hpp"
#include "Vec.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

template <int D>
class VtiWriter {
    static_assert(D == 2 || D == 3, "VtiWriter supports 2-D and 3-D meshes");

public:
    explicit VtiWriter(const Mesh<D>& mesh) : mesh_(mesh) {}

    void addScalarField(const std::string& name, const std::vector<double>& values) {
        if (static_cast<int>(values.size()) != mesh_.numNodes()) {
            throw std::runtime_error("VtiWriter: scalar field '" + name +
                                     "' size does not match node count");
        }
        scalars_.push_back({name, values});
    }

    // Vector fields are always written with three components, since that is
    // what ParaView expects. A 2-D field gets zero in the third slot.
    void addVectorField(const std::string& name, const std::vector<Vec<D>>& values) {
        if (static_cast<int>(values.size()) != mesh_.numNodes()) {
            throw std::runtime_error("VtiWriter: vector field '" + name +
                                     "' size does not match node count");
        }
        std::vector<double> flat;
        flat.reserve(values.size() * 3);
        for (const auto& v : values) {
            flat.push_back(v[0]);
            flat.push_back(v[1]);
            flat.push_back(D == 3 ? v[D - 1] : 0.0);
        }
        vectors_.push_back({name, flat});
    }

    void clear() {
        scalars_.clear();
        vectors_.clear();
    }

    void write(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("VtiWriter: cannot open " + filename);
        out << std::setprecision(9);

        const int nx = mesh_.count(0);
        const int ny = mesh_.count(1);
        const int nz = (D == 3) ? mesh_.count(D - 1) : 1;

        const double dz = (D == 3) ? mesh_.spacing(D - 1) : 1.0;
        const double z0 = (D == 3) ? mesh_.origin(D - 1) : 0.0;

        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "  <ImageData WholeExtent=\"0 " << nx - 1
            << " 0 " << ny - 1 << " 0 " << nz - 1 << "\"\n";
        out << "             Origin=\"" << mesh_.origin(0) << ' '
            << mesh_.origin(1) << ' ' << z0 << "\"\n";
        out << "             Spacing=\"" << mesh_.spacing(0) << ' '
            << mesh_.spacing(1) << ' ' << dz << "\">\n";
        out << "    <Piece Extent=\"0 " << nx - 1
            << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";

        out << "      <PointData";
        if (!scalars_.empty()) out << " Scalars=\"" << scalars_.front().name << "\"";
        if (!vectors_.empty()) out << " Vectors=\"" << vectors_.front().name << "\"";
        out << ">\n";

        for (const auto& f : scalars_) {
            out << "        <DataArray type=\"Float64\" Name=\"" << f.name
                << "\" format=\"ascii\">\n          ";
            for (double v : f.values) out << v << ' ';
            out << "\n        </DataArray>\n";
        }

        for (const auto& f : vectors_) {
            out << "        <DataArray type=\"Float64\" Name=\"" << f.name
                << "\" NumberOfComponents=\"3\" format=\"ascii\">\n          ";
            for (double v : f.values) out << v << ' ';
            out << "\n        </DataArray>\n";
        }

        out << "      </PointData>\n";
        out << "    </Piece>\n";
        out << "  </ImageData>\n";
        out << "</VTKFile>\n";
    }

private:
    struct Field {
        std::string name;
        std::vector<double> values;
    };

    const Mesh<D>& mesh_;
    std::vector<Field> scalars_;
    std::vector<Field> vectors_;
};

// Zero-padded time-series filename, e.g. "field_0042.vti".
inline std::string makeVtiFilename(const std::string& prefix, int step, int width = 4) {
    std::string s = std::to_string(step);
    while (static_cast<int>(s.size()) < width) s.insert(s.begin(), '0');
    return prefix + "_" + s + ".vti";
}
