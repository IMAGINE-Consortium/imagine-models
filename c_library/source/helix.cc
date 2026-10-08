#include "ImagineModels/Helix.h"
#include "ImagineModels/units.h"
#include <cmath>

namespace imagine {

template <typename T>
Vec3<T> HelixMagneticField::field(const double &x, const double &y, const double &z,
                                  const HelixParameters<T> &p) const {

    const double phi = std::atan2(y, x);       // azimuthal angle in cylindrical coordinates
    const double r = std::sqrt(x * x + y * y); // radius in cylindrical coordinates
    Vec3<T> b{{0.0, 0.0, 0.0}};
    if ((r > rmin) && (r < rmax)) {
        b[0] = std::cos(phi) * p.ampx;
        b[1] = std::sin(phi) * p.ampy;
        b[2] = p.ampz;
    }
    return b;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(HelixMagneticField)

}
