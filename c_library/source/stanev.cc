#include "ImagineModels/Stanev.h"
#include "ImagineModels/units.h"
#include <cmath>

#include "ImagineModels/helpers.h"

namespace imagine {

template <typename T>
Vec3<T> StanevBSSMagneticField::field(const double &x, const double &y, const double &z,
                                      const StanevBSSParameters<T> &p) const {

    Vec3<T> B_vec3{{0, 0, 0}};
    const double r = sqrt(x * x + y * y);
    const double phi = atan2(y, x);

    if (r > b_r_max || r == 0.) {
        return B_vec3;
    }

    auto phi_prime = p.b_phi0 - phi; // clockwise from negative x-axis
    auto beta = 1. / tan(p.b_p * (M_PI / 180.));

    auto B_0 = 3 * p.b_Rsun / b_r_min;
    if (r > b_r_min) {
        B_0 = 3 * p.b_Rsun / r;
    }

    auto z_0 = p.b_z01;
    if (std::abs(z) > p.b_z0_border) {
        z_0 = p.b_z02;
    }
    // eq. 1, 3, 4
    Vec3<T> B_cyl{
        {B_0 * cos(phi_prime - beta * log(r / p.b_r0)) * sin(p.b_p * (M_PI / 180.)) * exp(-std::abs(z) / z_0),
         -B_0 * cos(phi_prime - beta * log(r / p.b_r0)) * cos(p.b_p * (M_PI / 180.)) * exp(-std::abs(z) / z_0), 0.}};

    B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);
    return B_vec3;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(StanevBSSMagneticField)

}
