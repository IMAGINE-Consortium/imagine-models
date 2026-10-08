#include "ImagineModels/units.h"
#include <cmath>

#include "ImagineModels/TT.h"
#include "ImagineModels/helpers.h"

namespace imagine {

template <typename T>
Vec3<T> TTMagneticField::field(const double &x, const double &y, const double &z, const TTParameters<T> &p) const {

    double r = sqrt(x * x + y * y);
    double phi = atan2(y, x);

    Vec3<T> B_vec3{{0, 0, 0}};
    if (r > b_r_max || r == 0.) {
        return B_vec3;
    }

    auto pitch = p.b_p * units::deg;

    auto beta = 1. / tan(pitch);

    auto phase = (beta * log(1. + p.b_d / p.b_Rsun)) - units::pi / 2.; // eq. 2

    double sign = 1.;
    if (z < 0) { // antisymmetric halo, eq. 5
        sign = -1.;
    }

    auto f_z = sign * exp(-(std::abs(z) / p.b_z0)); // eq. 5

    // 1/cos(phase) as in Kachelriess
    auto b_r = p.b_b0 * (p.b_Rsun / (b_r_min * cos(phase)));
    if (r > b_r_min) {
        b_r = p.b_b0 * (p.b_Rsun / (r * cos(phase)));
    }

    else {
        b_r = p.b_b0 * (p.b_Rsun / (b_r_min * cos(phase))); // eq. 3
    }

    const double theta = units::pi - phi;
    auto B_r_phi = b_r * cos(theta - beta * log(r / p.b_Rsun) + phase); // eq. 1

    // cylindrical components
    Vec3<T> B_cyl{{B_r_phi * sin(pitch) * f_z, -B_r_phi * cos(pitch) * f_z, 0.}};

    B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);

    return B_vec3;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(TTMagneticField)

}
