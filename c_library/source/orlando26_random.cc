#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModelsRandom/Orlando26Random.h"

namespace imagine {

void Orlando26RandomField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown Orlando26 model '" + model + "'.");
    active_model = model;
    regular_base = XH24MagneticField();
    independent_ordered = true;
    b_ran = 4.9;
    r0_ran = 30.;
    z0_ran = 4.;
    r_sun = 8.5;
    b_ordered = model == "halo4kpc" ? 4.3 : 3.2; // Table 1
}

double Orlando26RandomField::spectrum(const double &abs_k) const {
    return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

double Orlando26RandomField::isotropic_rms(const double &x, const double &y, const double &z) const {
    const double r = std::sqrt(x * x + y * y);
    return b_ran * std::exp(-(r - r_sun) / r0_ran) * std::exp(-std::abs(z) / z0_ran);
}

double Orlando26RandomField::ordered_amplitude(const double &x, const double &y, const double &z) const {
    const double b_h = regular_base.parameters.B0;
    const double b_or = std::sqrt(std::max(b_ordered * b_ordered - b_h * b_h, 0.));
    const Vec3<double> b = regular_base.at_position(x, y, z);
    const double shape = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]) / b_h; // eq. 2
    return std::sqrt(3.) * b_or * shape;
}

Vec3<double> Orlando26RandomField::anisotropy_direction(const double &x, const double &y, const double &z) const {
    return regular_base.at_position(x, y, z);
}

double Orlando26RandomField::rms(const double &x, const double &y, const double &z) const {
    return combined_rms(isotropic_rms(x, y, z), apply_anisotropy ? ordered_amplitude(x, y, z) : 0.);
}

}
