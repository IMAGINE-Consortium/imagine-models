#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModelsRandom/UF26Random.h"

namespace imagine {

namespace {

double gaussian(const double x, const double sigma) {
    return std::exp(-0.5 * x * x / (sigma * sigma));
}

}

void UF26RandomField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown UF26 model '" + model + "'.");
    active_model = model;
    r_max = 18.;
    w_max = 2.;
    r_sun = 8.178;
    r_ring = 4.5;
    w_ring = 2.0;
    z_ring = 1.9;
    if (model == "expDisk") {
        b_disk = 4.4;
        z_disk = 1.0;
        l_r = 15.;
        r_c = 1.;
        b_ring = 0.;
    } else {
        b_disk = 4.1;
        z_disk = 0.9;
        l_r = 100.;
        r_c = 0.;
        b_ring = 4.1;
    }
}

double UF26RandomField::spectrum(const double &abs_k) const {
    return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

double UF26RandomField::rms(const double &x, const double &y, const double &z) const {
    const double r = std::sqrt(x * x + y * y);
    const double vertical = 1. / std::cosh(std::acosh(std::exp(1.)) * z / z_disk);
    const double sigmoid = 1. / (1. + std::exp((r - r_max) / w_max));
    const auto radial = [&](double rr) { return std::exp(-std::sqrt(rr * rr + r_c * r_c) / l_r); };
    double b_d = b_disk * vertical * sigmoid * radial(r) / radial(r_sun);
    if (active_model == "expDisk")
        return b_d;
    b_d *= gaussian(std::min(r - r_ring, 0.), w_ring);
    const double b_a = b_ring * gaussian(z, z_ring) * gaussian(r - r_ring, w_ring);
    return std::sqrt(b_d * b_d + b_a * b_a);
}

}
