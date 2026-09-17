#pragma once
//
// Vec.hpp
//
// A small fixed-size vector templated on dimension.
//
// This is the piece that makes "2-D now, 3-D later" a one-line change instead
// of a rewrite. Every physics routine is written against Vec<D>, so switching
// kDim from 2 to 3 recompiles the whole solver in three dimensions.
//
// Note on cross products: the magnetic force term v_e x B needs special care
// in 2-D. In an r-z simulation the azimuthal (theta) direction is not on the
// mesh, but B and the E x B drift both have components there. cross2d() below
// returns the out-of-plane scalar component; crossWithOutOfPlane() applies an
// out-of-plane field back onto in-plane velocity. Read the comments before
// using either -- this is where 2-D magnetised plasma codes usually go wrong.
//
// DCFHT Simulation Team
//

#include <array>
#include <cmath>
#include <cstddef>

template <int D>
struct Vec {
    std::array<double, D> c{};

    constexpr Vec() = default;

    // 2-D convenience constructor.
    constexpr Vec(double a, double b) requires(D == 2) : c{a, b} {}

    // 3-D convenience constructor.
    constexpr Vec(double a, double b, double d) requires(D == 3) : c{a, b, d} {}

    double&       operator[](int i)       { return c[static_cast<std::size_t>(i)]; }
    const double& operator[](int i) const { return c[static_cast<std::size_t>(i)]; }

    Vec& operator+=(const Vec& o) { for (int i = 0; i < D; ++i) c[i] += o.c[i]; return *this; }
    Vec& operator-=(const Vec& o) { for (int i = 0; i < D; ++i) c[i] -= o.c[i]; return *this; }
    Vec& operator*=(double s)     { for (int i = 0; i < D; ++i) c[i] *= s;      return *this; }

    double dot(const Vec& o) const {
        double s = 0.0;
        for (int i = 0; i < D; ++i) s += c[i] * o.c[i];
        return s;
    }

    double normSquared() const { return dot(*this); }
    double norm()        const { return std::sqrt(normSquared()); }

    static constexpr int dimension() { return D; }
};

template <int D> Vec<D> operator+(Vec<D> a, const Vec<D>& b) { a += b; return a; }
template <int D> Vec<D> operator-(Vec<D> a, const Vec<D>& b) { a -= b; return a; }
template <int D> Vec<D> operator*(Vec<D> a, double s)        { a *= s; return a; }
template <int D> Vec<D> operator*(double s, Vec<D> a)        { a *= s; return a; }

// --- magnetic-force helpers -----------------------------------------------
//
// In a 2-D (r,z) simulation of an axisymmetric thruster:
//   - B lies in the (r,z) plane   (this is the mirror field you care about)
//   - v x B for an in-plane velocity therefore points out of plane (azimuthal)
//   - the azimuthal Hall drift, crossed with in-plane B, pushes back in plane
//
// So you need both directions of the operation. Keeping them as two named
// functions makes the physics explicit at every call site.

// In-plane a crossed with in-plane b -> scalar out-of-plane (azimuthal) component.
inline double crossOutOfPlane(const Vec<2>& a, const Vec<2>& b) {
    return a[0] * b[1] - a[1] * b[0];
}

// Out-of-plane scalar w crossed with in-plane b -> in-plane result.
inline Vec<2> crossFromOutOfPlane(double w, const Vec<2>& b) {
    return Vec<2>{ w * b[1], -w * b[0] };
}

// Full 3-D cross product, used when kDim == 3.
inline Vec<3> cross(const Vec<3>& a, const Vec<3>& b) {
    return Vec<3>{ a[1] * b[2] - a[2] * b[1],
                   a[2] * b[0] - a[0] * b[2],
                   a[0] * b[1] - a[1] * b[0] };
}
