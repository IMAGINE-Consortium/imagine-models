#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModelsRandom/JaffeRandom.h"

namespace imagine {

void JaffeRandomField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown Jaffe random model '" + model + "'.");
    active_model = model;
    regular_base.set_model(model);
    b_rms = 3.5;
    h_rms = 2.;
    r_grf = 20.;
    f_ord = 0.15;
    spectral_offset = 0.;
    spectral_slope = 0.37;
    k_min = 10.; // 1/D_co
}

double JaffeRandomField::spectrum(const double &abs_k) const {
    return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

bool JaffeRandomField::outside(const double &x, const double &y, const double &z) const {
    return regular_base.r_max > 0. && std::sqrt(x * x + y * y + z * z) > regular_base.r_max;
}

double JaffeRandomField::compression(const double &x, const double &y, const double &z) const {
    const auto &p = regular_base.parameters;
    const std::vector<double> arm = regular_base.arm_compress<double>(x, y, z, p);
    const double b_r = regular_base.radial_scaling<double>(x, y, p);
    double sum = 0.;
    for (std::size_t i = 0; i < arm.size(); ++i) {
        const double weight = i + 1 < arm.size() ? 1. : std::abs(p.ring_amp);
        sum += weight * arm[i] / b_r; // rho_c without B(r)
    }
    return sum;
}

double JaffeRandomField::isotropic_rms(const double &x, const double &y, const double &z) const {
    if (outside(x, y, z))
        return 0.;
    const double r2 = x * x + y * y;
    const double sech = 1. / std::cosh(z / h_rms);
    return b_rms * (sech * sech * std::exp(-r2 / (r_grf * r_grf)) + compression(x, y, z));
}

double JaffeRandomField::ordered_amplitude(const double &x, const double &y, const double &z) const {
    if (outside(x, y, z))
        return 0.;
    return b_rms * f_ord * compression(x, y, z);
}

Vec3<double> JaffeRandomField::anisotropy_direction(const double &x, const double &y, const double &z) const {
    return regular_base.orientation<double>(x, y, z, regular_base.parameters);
}

double JaffeRandomField::rms(const double &x, const double &y, const double &z) const {
    const Vec3<double> e = anisotropy_direction(x, y, z);
    const bool ordered = apply_anisotropy && (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]) > 1e-20;
    return combined_rms(isotropic_rms(x, y, z), ordered ? ordered_amplitude(x, y, z) : 0.);
}

}
