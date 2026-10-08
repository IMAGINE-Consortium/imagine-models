#include "ImagineModels/Han.h"

#include "ImagineModels/helpers.h"

#include <algorithm>
#include <stdexcept>

namespace imagine {

void HanMagneticField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown Han model '" + model + "'.");
    active_model = model;
    parameters = HanParameters<double>{};
    R_min = 3.;
    R_max = 15.;
    R_s = {3.0, 4.1, 4.9, 6.1, 7.5, 8.5, 10.5, 15.0};
    active_parameters = {"B_p", "A", "H", "B_s1", "B_s2", "B_s3", "B_s4", "B_s5", "B_s6"};
    if (model == "XH24") {
        R_s[5] = 8.16;
        R_max = 20.;
        parameters.B_s7 = 4.5;
        active_parameters.push_back("B_s7");
    }
}

// J. L. Han et al 2018 ApJS 234 11
template <typename T>
Vec3<T> HanMagneticField::field(const double &x, const double &y, const double &z, const HanParameters<T> &p) const {
    Vec3<T> B_cyl{{0., 0., 0.}};
    const double r = sqrt(x * x + y * y);
    const double phi = atan2(y, x);

    if (r < R_min || r > R_max)
        return B_cyl;

    T B_0 = 0.;

    auto p_ang = p.B_p * M_PI / 180.;
    const double phi_han = -(phi + M_PI); // nneeded to fix different coordinate system convention

    T R_0 = r * exp(phi_han * tan(p_ang)); // eq. 4 is wrong, need to change psi and phi!

    std::array<T, 7> B_s = {p.B_s1, p.B_s2, p.B_s3, p.B_s4, p.B_s5, p.B_s6, p.B_s7}; // table 5

    if (R_0 < R_s[0]) {
        R_0 = r * exp((phi_han + 2 * M_PI) * tan(p_ang)); // eq. 4
    }
    if (R_0 > R_s[6]) {
        R_0 = r * exp((phi_han - 2 * M_PI) * tan(p_ang)); // eq. 4
    }

    for (int i = 0; i < 7; i++) {
        if (R_s[i] < R_0) {
            if (R_0 < R_s[i + 1]) {
                B_0 = B_s[i];
                break;
            }
        }
    }

    T B_r = B_0 * exp(-r / p.A) * exp(-std::abs(z) / p.H); // eq. 3

    B_cyl[0] = B_r * sin(p_ang);
    B_cyl[1] = B_r * cos(p_ang);

    Vec3<T> B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);
    return B_vec3;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(HanMagneticField)

}
