#include "ImagineModels/Fauvet.h"
#include "ImagineModels/helpers.h"
#include "ImagineModels/units.h"
#include <cmath>

namespace imagine {

template <typename T>
Vec3<T> FauvetMagneticField::field(const double &x, const double &y, const double &z,
                                   const FauvetParameters<T> &p) const {
    Vec3<T> B_vec3{{0, 0, 0}};
    const double r = sqrt(x * x + y * y);

    if (r > b_r_max || r < b_r_min) {
        return B_vec3;
    }

    double phi = atan2(y, x);
    auto chi_z = p.b_chi0 * units::deg * tanh(z / p.b_z0);
    auto beta = 1. / tan(p.b_p * units::deg);

    auto b_r = p.b_b0 * exp(-(r - p.b_Rsun) / p.b_RB);

    // cylindrical components
    Vec3<T> B_cyl{{b_r * cos(phi + beta * log(r / p.b_r0)) * sin(p.b_p * units::deg) * cos(chi_z),
                   -b_r * cos(phi + beta * log(r / p.b_r0)) * cos(p.b_p * units::deg) * cos(chi_z), b_r * sin(chi_z)}};

    B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);
    return B_vec3;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(FauvetMagneticField)

}
