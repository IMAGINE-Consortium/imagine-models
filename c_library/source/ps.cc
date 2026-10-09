#include <algorithm>
#include <cmath>

#include "ImagineModels/PS.h"
#include "ImagineModels/helpers.h"
#include "ImagineModels/units.h"

namespace imagine {

template <typename T>
Vec3<T> PSMagneticField::field(const double &x, const double &y, const double &z, const PSParameters<T> &p) const {

    const double r = std::sqrt(x * x + y * y);
    const double R = std::sqrt(x * x + y * y + z * z);
    const double phi = std::atan2(y, x);

    Vec3<T> B{{0, 0, 0}};

    // BSS-S disk
    if (do_disk && r > 0. && r <= b_r_max) {
        auto pitch = p.b_p * units::deg;
        auto beta = 1. / tan(pitch);
        auto phase = beta * log(1. + p.b_d / p.b_Rsun) - units::pi / 2.;
        auto b_r = p.b_b0 * p.b_Rsun / (std::max(r, b_r_min) * cos(phase));
        const double theta = units::pi - phi;
        auto B_rt = b_r * cos(theta - beta * log(r / p.b_Rsun) + phase) * exp(-std::abs(z) / p.b_z0); // eqs. 3, 4
        Vec3<T> B_cyl{{B_rt * sin(pitch), -B_rt * cos(pitch), 0.}};
        B = Cyl2Cart<Vec3<T>>(phi, B_cyl);
    }

    // toroidal halo, eqs. 7-9
    if (do_halo && r > 0. && z != 0.) {
        T b_max = p.h_b0;
        if (r > h_R0) {
            b_max = p.h_b0 * std::exp((h_R0 - r) / h_R0);
        }
        auto u = (std::abs(z) - p.h_z0) / p.h_w;
        auto b_t = (z > 0. ? 1. : -1.) * b_max / (1. + u * u);
        B[0] += -b_t * std::sin(phi);
        B[1] += b_t * std::cos(phi);
    }

    // dipole, eq. 10
    if (do_dipole) {
        if (R < d_r_core) {
            B[2] += d_b_core;
        } else {
            const double R5 = R * R * R * R * R;
            B[0] += -3. * p.d_mu * z * x / R5;
            B[1] += -3. * p.d_mu * z * y / R5;
            B[2] += p.d_mu * (R * R - 3. * z * z) / R5;
        }
    }

    return B;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(PSMagneticField)

}
