#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModelsRandom/SunRandom.h"

namespace imagine {

void SunRandomField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown Sun random model '" + model + "'.");
    active_model = model;
    r_sun = 8.5;
    r0 = 30.;
    h_disk = 1.;
    h_halo = 3.;
    f_disk = 0.5;
    b_iso = (model == "Sun10b") ? 6.4 : 3.;
}

double SunRandomField::spectrum(const double &abs_k) const {
    return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

double SunRandomField::rms(const double &x, const double &y, const double &z) const {
    if (active_model == "Sun10")
        return b_iso;
    const double r = std::sqrt(x * x + y * y);
    const auto sech2 = [](double u) { return 1. / (std::cosh(u) * std::cosh(u)); };
    const double vertical = (1. - f_disk) * sech2(z / h_halo) + f_disk * sech2(z / h_disk);
    return b_iso * std::exp(-(r - r_sun) / r0) * vertical; // Table C.1
}

}
