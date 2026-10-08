#include "ImagineModels/SVT22.h"
#include "ImagineModels/units.h"
#include <cassert>
#include <cmath>
#include <iostream>

namespace imagine {

template <typename T>
Vec3<T> SVT22MagneticField::field(const double &x, const double &y, const double &z,
                                  const SVT22Parameters<T> &p) const {
    const double r{sqrt(x * x + y * y)};
    const double rho{sqrt(x * x + y * y + z * z)};
    const double phi{atan2(y, x)};

    T B_cyl[3] = {0, 0, 0}; // cylindrical components

    // toroidal halo

    if (do_halo) {
        T b1, rh;
        T B_h = 0.;
        const double z_min = 0.1;
        if (z >= 0) { // North
            b1 = p.B_val;
        } else { // South
            b1 = -p.B_val;
        }

        B_h = b1 * (exp(-z_min / std::abs(z)) * exp(-std::abs(r) / p.r_cut) *
                    exp(-(std::abs(z)) / (p.z_cut))); // vertical exponential fall-off
        const T B_cyl_h[3] = {0., B_h * 1, 0.};
        // add fields together
        B_cyl[0] += B_cyl_h[0];
        B_cyl[1] += B_cyl_h[1];
        B_cyl[2] += B_cyl_h[2];
    }

    // convert field to cartesian coordinates
    Vec3<T> B_cart{{0.0, 0.0, 0.0}};
    B_cart[0] = B_cyl[0] * cos(phi) - B_cyl[1] * sin(phi);
    B_cart[1] = B_cyl[0] * sin(phi) + B_cyl[1] * cos(phi);
    B_cart[2] = B_cyl[2];
    return B_cart;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(SVT22MagneticField)

}
