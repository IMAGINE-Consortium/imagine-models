#include "ImagineModels/units.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModelsRandom/JF12Random.h"

namespace imagine {

void JF12RandomField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown JF12 model '" + model + "'.");
    active_model = model;
    regular_base.set_model(model);
    b0_1 = 10.81;
    b0_2 = 6.96;
    b0_3 = 9.59;
    b0_4 = 6.96;
    b0_5 = 1.96;
    b0_6 = 16.34;
    b0_7 = 37.29;
    b0_8 = 10.35;
    b0_int = 7.63;
    b0_halo = 4.68;
    arm_shift = 1.;
    if (model == "JF12")
        return;
    const double b_iso = 7.8;
    b0_1 = b0_3 = b0_5 = b0_7 = 0.4 * b_iso;
    b0_2 = b0_4 = b0_6 = b0_8 = 0.8 * b_iso;
    b0_int = 0.5 * b_iso;
    b0_halo = 0.94 * b_iso;
    if (model == "Planck12c") {
        b0_6 = 1.6 * b_iso;
        b0_1 = b0_3 = b0_7 = 0.1 * b_iso;
        b0_5 = 0.;
        arm_shift = 0.97;
    }
}

double JF12RandomField::spectrum(const double &abs_k) const {
    return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

Vec3<double> JF12RandomField::anisotropy_direction(const double &x, const double &y, const double &z) const {
    return regular_base.at_position(x, y, z);
}

double JF12RandomField::rms(const double &x, const double &y, const double &z) const {

    const double r{sqrt(x * x + y * y)};
    const double rho{sqrt(x * x + y * y + z * z)};
    const double phi{atan2(y, x)};

    const double b_arms[8] = {b0_1, b0_2, b0_3, b0_4, b0_5, b0_6, b0_7, b0_8};

    double scaling_disk = 0.0;
    double scaling_halo = 0.0;

    if (r > Rmax || rho < rho_GC) {
        return 0.0;
    }
    if (r < 5.) {
        scaling_disk = b0_int;
    } else {
        double r_negx = r * exp(-1 / tan(units::deg * (90 - inc)) * (phi - units::pi));
        if (r_negx > rc_B[7] * arm_shift) {
            r_negx = r * exp(-1 / tan(units::deg * (90 - inc)) * (phi + units::pi));
        }
        if (r_negx > rc_B[7] * arm_shift) {
            r_negx = r * exp(-1 / tan(units::deg * (90 - inc)) * (phi + 3 * units::pi));
        }
        for (int i = 7; i >= 0; i--) {
            if (r_negx < rc_B[i] * arm_shift) {
                scaling_disk = b_arms[i] * (5.) / r;
            }
        } // "region 8,7,6,..,2"
    }

    scaling_disk = scaling_disk * exp(-0.5 * z * z / (z0_spiral * z0_spiral));
    scaling_halo = b0_halo * exp(-r / r0_halo) * exp(-0.5 * z * z / (z0_halo * z0_halo));

    return std::sqrt(scaling_disk * scaling_disk + scaling_halo * scaling_halo);
}

}
