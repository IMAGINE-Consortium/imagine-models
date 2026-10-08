#include "ImagineModels/WMAP.h"
#include "ImagineModels/units.h"
#include <cmath>

#include "ImagineModels/helpers.h"

namespace imagine {

template <typename T>
Vec3<T> WMAPMagneticField::field(const double &x, const double &y, const double &z, const WMAPParameters<T> &p) const {

    Vec3<T> B_vec3{{0, 0, 0}};
    double r = sqrt(x * x + y * y);

    if (r > b_r_max || r < b_r_min) {
        return B_vec3;
    }

    double phi = atan2(y, x);

    auto psi_r = p.b_psi0 * units::deg + p.b_psi1 * units::deg * log(r / p.b_r0);
    auto xsi_z = p.b_xsi0 * units::deg * tanh(z / p.b_z0);

    Vec3<T> B_cyl{{p.b_b0 * sin(psi_r) * cos(xsi_z), // eq. 9
                   p.b_b0 * cos(psi_r) * cos(xsi_z), p.b_b0 * sin(xsi_z)}};

    B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);

    // antisymmetric variant
    if (anti && z > 0) {
        B_vec3[0] *= (-1.);
        B_vec3[1] *= (-1.);
        B_vec3[2] *= (-1.);
    }
    return B_vec3;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(WMAPMagneticField)

}
