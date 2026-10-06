#include <cmath>
#include "ImagineModels/units.h"
#include "ImagineModels/Fauvet.h"
#include "ImagineModels/helpers.h"

namespace imagine {

// ??????, implementation from Hammurabi (old)

template <typename T>
Vec3<T> FauvetMagneticField::field(const double &x, const double &y, const double &z,  const FauvetParameters<T> &p) const
{
    Vec3<T> B_vec3{{0, 0, 0}};
    const double r = sqrt(x * x + y * y);

    if (r > b_r_max || r < b_r_min)
    {
        return B_vec3;
    }

    double phi = atan2(y, x);
    auto chi_z = p.b_chi0 * (M_PI / 180.) * tanh(z / p.b_z0);
    auto beta = 1. / tan(p.b_p * (M_PI / 180.));

    // B-field in cylindrical coordinates:
    Vec3<T> B_cyl{{p.b_b0 * cos(phi + beta * log(r / p.b_r0)) * sin(p.b_p * (M_PI / 180.)) * cos(chi_z),
                  -p.b_b0 * cos(phi + beta * log(r / p.b_r0)) * cos(p.b_p * (M_PI / 180.)) * cos(chi_z),
                  p.b_b0 * sin(chi_z)}};

    // Taking into account the halo field
    T h_z1;
    if (std::abs(z) < p.h_z0)
    {
        h_z1 = p.h_z1a;
    }
    else
    {
        h_z1 = p.h_z1b;
    }

    auto hf_piece1 = (h_z1 * h_z1) / (h_z1 * h_z1 + (std::abs(z) - p.h_z0) * (std::abs(z) - p.h_z0));
    auto hf_piece2 = exp(-(r - p.h_r0) / (p.h_r0));

    auto halo_field = p.h_b0 * hf_piece1 * (r / p.b_r0) * hf_piece2;
    B_cyl[1] += halo_field;

    B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);
    return B_vec3;
}


IMAGINE_INSTANTIATE_VECTOR_MODEL(FauvetMagneticField)

}
