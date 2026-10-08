#include <array>
#include <cmath>

#include "ImagineModels/YT20.h"
#include "ImagineModels/units.h"

namespace imagine {

namespace {

constexpr double G = 6.67430e-8;
constexpr double M_sun = 1.98841e33;
constexpr double m_p = 1.67262192e-24;
constexpr double keV = 1.602176634e-9;
constexpr double kpc_cm = 3.0856775814913673e21;
constexpr int n_nodes = 64;

struct GaussLegendre {
    std::array<double, n_nodes> x, w;
    GaussLegendre() {
        for (int i = 0; i < n_nodes; ++i) {
            double t = std::cos(units::pi * (i + 0.75) / (n_nodes + 0.5)), dp = 0.;
            for (int iteration = 0; iteration < 100; ++iteration) {
                double p0 = 1., p1 = t;
                for (int k = 2; k <= n_nodes; ++k) {
                    const double p2 = ((2. * k - 1.) * t * p1 - (k - 1.) * p0) / k;
                    p0 = p1;
                    p1 = p2;
                }
                dp = n_nodes * (t * p1 - p0) / (t * t - 1.);
                const double dt = p1 / dp;
                t -= dt;
                if (std::abs(dt) < 1e-16)
                    break;
            }
            x[i] = t;
            w[i] = 2. / ((1. - t * t) * dp * dp);
        }
    }
};

const GaussLegendre &quadrature() {
    static const GaussLegendre q;
    return q;
}

// hydrostatic profile, eq. 3
template <typename T> T profile(const T &r, const T &r_s, const T &upsilon) {
    const T s = r / r_s;
    const T ratio = value(s) < 1e-8 ? T(1. - s / 2.) : T(log(1. + s) / s);
    return exp(-upsilon * (1. - ratio));
}

}

template <typename T>
T YT20::field(const double &x, const double &y, const double &z, const YT20Parameters<T> &p) const {
    const double r = std::sqrt(x * x + y * y + z * z);
    if (r > p.r_vir)
        return T(0.);

    // spherical component, eqs. 3-4
    const T r_s = p.r_vir / p.c_NFW;
    const T f_c = log(1. + p.c_NFW) - p.c_NFW / (1. + p.c_NFW);
    const T upsilon = G * p.M_vir * 1e12 * M_sun * mu * m_p / (r_s * kpc_cm * f_c * p.T_halo * keV);
    const auto &q = quadrature();
    T integral = 0.;
    for (int i = 0; i < n_nodes; ++i) {
        const T ri = p.r_vir * (q.x[i] + 1.) / 2.;
        integral += q.w[i] * 4. * units::pi * ri * ri * profile(ri, r_s, upsilon);
    }
    integral *= p.r_vir / 2.;
    const T n0_sphere = p.M_b * 1e11 * M_sun / (mu_e * m_p * integral * kpc_cm * kpc_cm * kpc_cm);
    const T sphere = n0_sphere * profile(T(r), r_s, upsilon);

    // disk-like component, eq. 2
    const double R = std::sqrt(x * x + y * y);
    const T disk = p.n0_disk / p.Z_halo * exp(-R / p.R0_disk - std::abs(z) / p.z0_disk);
    return sphere + disk;
}

IMAGINE_INSTANTIATE_SCALAR_MODEL(YT20)

}
