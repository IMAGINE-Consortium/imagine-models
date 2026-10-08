#include <cmath>

#include "ImagineModels/KST24.h"
#include "ImagineModels/units.h"

namespace imagine {

namespace {

template <typename T> T pow_abs(const T &u, const T &n) {
    if (value(u) == 0.)
        return T(0.);
    return pow(abs(u), n);
}

template <typename T> void add(Vec3<T> &a, const Vec3<T> &b) {
    for (int c = 0; c < 3; ++c)
        a[c] += b[c];
}

}

template <typename T>
Vec3<T> KST24MagneticField::arm(const double &x, const double &y, const double &z, const T &B, const T &pitch,
                                const T &phase, const T &x0, const T &rz, const T &rdisk, const T &n, int i) const {
    Vec3<T> b{{0., 0., 0.}};
    if (std::abs(z) > arm_zmax[i])
        return b;
    const T xs = x + x0;
    const T r = sqrt(xs * xs + y * y);
    if (r < arm_rmin[i] || r > arm_rmax[i])
        return b;

    const T k = tan(pitch * units::deg);
    const T cos_pitch = cos(pitch * units::deg);
    const T sin_pitch = sin(pitch * units::deg);
    const T phi0 = phase * units::deg;
    const T phi = atan2(T(y), xs);
    int nn = int(std::floor(value((log(r / spiral_a) - k * (phi + phi0)) / (units::twopi * k))));

    // closest arm crossing
    T r1 = spiral_a * exp(k * (phi + phi0 + units::twopi * nn));
    const T r2 = r1 * exp(k * units::twopi);
    if (abs(r - r1) > abs(r - r2)) {
        r1 = r2;
        nn++;
    }
    T xi = r1 * cos(phi);
    T yi = r1 * sin(phi);
    T d1 = k * xi - yi;
    T d2 = k * yi + xi;
    T norm = sqrt(d1 * d1 + d2 * d2);
    d1 /= norm;
    d2 /= norm;

    // first-order axis correction
    const T delta_phi = -((xi - xs) * d1 + (yi - y) * d2) / (r1 / cos_pitch);
    r1 = spiral_a * exp(k * (phi + phi0 + delta_phi + units::twopi * nn));
    xi = r1 * cos(phi + delta_phi);
    yi = r1 * sin(phi + delta_phi);
    d1 = k * xi - yi;
    d2 = k * yi + xi;
    norm = sqrt(d1 * d1 + d2 * d2);
    d1 /= norm;
    d2 /= norm;

    // squircle cross section
    const T distance = sqrt((xi - xs) * (xi - xs) + (yi - y) * (yi - y));
    const T widening = (r1 / sin_pitch - width_reference_r / sin_pitch) * arm_widening[i] * units::deg;
    T width_plane = rdisk * (1. + widening) > rdisk ? T(rdisk * (1. + widening)) : rdisk;
    T width_z = rz * (1. + widening) > rz ? T(rz * (1. + widening)) : rz;
    if (width_plane > width_max)
        width_plane = width_max;
    if (width_z > width_max)
        width_z = width_max;
    if (pow_abs(T(distance / width_plane), n) + pow_abs(T(z / width_z), n) > 1.)
        return b;

    const T scale = rz * rdisk / (width_plane * width_z);
    b[0] = B * d1 * scale;
    b[1] = B * d2 * scale;
    return b;
}

template <typename T>
Vec3<T> KST24MagneticField::torus(const double &x, const double &y, const double &z, const T &B, const T &zmin,
                                  const T &zmax, const T &rmax) const {
    Vec3<T> b{{0., 0., 0.}};
    const double r = std::sqrt(x * x + y * y);
    if (z <= zmin || zmax <= z || r < torus_rmin || r > rmax)
        return b;
    b[0] = -B * y / r;
    b[1] = B * x / r;
    return b;
}

template <typename T>
Vec3<T> KST24MagneticField::xfield(const double &x, const double &y, const double &z,
                                   const KST24Parameters<T> &p) const {
    Vec3<T> b{{0., 0., 0.}};
    if (std::abs(z) > X_zmax)
        return b;
    const double sign = z < 0. ? -1. : 1.;
    const double r = std::sqrt(x * x + y * y);
    const T theta = p.theta_X * units::deg;
    const T r0 = r - std::abs(z) * tan(theta);
    if (r0 < X_rmin || r0 > p.rmax_X)
        return b;
    const T strength = p.B_X * r0 / r;
    b[0] = strength * (x / r) * sign * sin(theta);
    b[1] = strength * (y / r) * sign * sin(theta);
    b[2] = strength * cos(theta);
    return b;
}

template <typename T>
Vec3<T> KST24MagneticField::bubble(const double &x, const double &y, const double &z,
                                   const KST24Parameters<T> &p) const {
    Vec3<T> b{{0., 0., 0.}};
    const Vec3<T> e{{x - p.x_LB, y - p.y_LB, z - p.z_LB}};
    const T distance = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
    if (distance < p.r_LB || distance > p.r_LB + p.dr_LB)
        return b;
    const T l = p.l_LB * units::deg;
    const T lat = p.b_LB * units::deg;
    const Vec3<T> direction{{cos(lat) * cos(l), cos(lat) * sin(l), sin(lat)}};
    T cos_theta = 0.;
    for (int c = 0; c < 3; ++c)
        cos_theta += e[c] / distance * direction[c];
    const T amplification = 1. + p.r_LB * p.r_LB / (2. * p.r_LB * p.dr_LB + p.dr_LB * p.dr_LB);
    for (int c = 0; c < 3; ++c)
        b[c] = abs(p.B_local) * amplification * (direction[c] - e[c] / distance * cos_theta);
    return b;
}

template <typename T>
Vec3<T> KST24MagneticField::field(const double &x, const double &y, const double &z,
                                  const KST24Parameters<T> &p) const {
    Vec3<T> b{{0., 0., 0.}};
    // northern and southern halo
    if (z > p.zmin_ntor || z < p.zmax_stor) {
        if (z > p.zmin_ntor)
            add(b, torus(x, y, z, p.B_ntor, p.zmin_ntor, p.zmax_ntor, p.rmax_ntor));
        else
            add(b, torus(x, y, z, p.B_stor, p.zmin_stor, p.zmax_stor, p.rmax_stor));
        add(b, xfield(x, y, z, p));
        return b;
    }
    // Local Bubble replaces the disk
    const T dx = x - p.x_LB, dy = y - p.y_LB, dz = z - p.z_LB;
    if (sqrt(dx * dx + dy * dy + dz * dz) < p.r_LB + p.dr_LB)
        return bubble(x, y, z, p);

    add(b, arm(x, y, z, p.B_local, p.pitch, p.b_local, p.x0_local, p.rz_local, p.rdisk_local, p.n_local, 0));
    add(b, arm(x, y, z, p.B_sagcar, p.pitch, p.b_sagcar, p.x0_sagcar, p.rz_sagcar, p.rdisk_sagcar, p.n_sagcar, 1));
    add(b, arm(x, y, z, p.B_scutum, p.pitch, p.b_scutum, p.x0_scutum, p.rz_scutum, p.rdisk_scutum, p.n_scutum, 2));
    add(b, arm(x, y, z, p.Bprime_perseus, p.pitch, p.b_perseus, p.x0_perseus, p.rzprime_perseus, p.rdisk_perseus,
               p.n_perseus, 3));
    add(b,
        arm(x, y, z, p.B_perseus, p.pitch, p.b_perseus, p.x0_perseus, p.rz_perseus, p.rdisk_perseus, p.n_perseus, 4));
    add(b, xfield(x, y, z, p));
    return b;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(KST24MagneticField)

}
