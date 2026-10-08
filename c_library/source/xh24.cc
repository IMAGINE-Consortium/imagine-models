#include "ImagineModels/XH24.h"

namespace imagine {

template <typename T>
Vec3<T> XH24MagneticField::field(const double &x, const double &y, const double &z, const XH24Parameters<T> &p) const {
    const double r = std::sqrt(x * x + y * y);
    if (r <= 0. || r > r_max || z == 0.)
        return Vec3<T>{{0., 0., 0.}};
    const double sign_z = z > 0. ? 1. : -1.;
    const double abs_z = std::abs(z);
    const T b_phi =
        sign_z * p.B0 * (abs_z / p.z0) * exp(-(abs_z - p.z0) / p.z0) * exp(-((r - p.R0) / p.RT) * ((r - p.R0) / p.RT));
    return Vec3<T>{{-b_phi * y / r, b_phi * x / r, 0.}};
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(XH24MagneticField)

}
