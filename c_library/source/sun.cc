#include "ImagineModels/Sun.h"
#include "ImagineModels/units.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModels/helpers.h"

namespace imagine {

void SunMagneticField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown Sun model '" + model + "'.");
    active_model = model;
    parameters = SunParameters<double>{};
    if (model == "Sun10b")
        parameters.b_Bc = 0.5;
}

template <typename T>
Vec3<T> SunMagneticField::field(const double &x, const double &y, const double &z, const SunParameters<T> &p) const {

    double r = sqrt(x * x + y * y);
    double phi = atan2(y, x);

    // D2, ASS+RING, eq. 8
    double D2;
    if (r > 7.5) {
        D2 = 1.;
    } else if (r <= 7.5 && r > 6.) {
        D2 = -1.;
    } else if (r <= 6. && r > 5.) {
        D2 = 1.;
    } else {
        D2 = -1.;
    }

    // D1, eq. 7
    T D1;
    if (r > p.b_Rc) {
        D1 = p.b_B0 * exp(-((r - p.b_Rsun) / p.b_R0) - (std::abs(z) / p.b_z0));
    } else {
        D1 = p.b_Bc;
    }

    auto p_ang = p.b_p * units::deg;
    Vec3<T> B_cyl{{D1 * D2 * sin(p_ang), // eq. 6
                   -D1 * D2 * cos(p_ang), 0.}};

    // halo field
    T halo_field;

    T b3H_z1_actual;
    if (std::abs(z) < p.bH_z0) {
        b3H_z1_actual = p.bH_z1a;
    } else {
        b3H_z1_actual = p.bH_z1b;
    }
    auto hf_piece1 = (b3H_z1_actual * b3H_z1_actual) /
                     (b3H_z1_actual * b3H_z1_actual + (std::abs(z) - p.bH_z0) * (std::abs(z) - p.bH_z0));
    auto hf_piece2 = exp(-(r - p.bH_R0) / (p.bH_R0));

    halo_field = p.bH_B0 * hf_piece1 * (r / p.bH_R0) * hf_piece2; // eq. 10

    // sign(z), Sun & Reich eq. 1
    if (z < 0) {
        halo_field *= -1.;
    }

    B_cyl[1] += halo_field;

    Vec3<T> B_vec3;

    B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);

    return B_vec3;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(SunMagneticField)

}
