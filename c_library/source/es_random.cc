#include <cmath>

#include "ImagineModelsRandom/ESRandom.h"

namespace imagine {

double ESRandomField::spectrum(const double &abs_k) const {
    return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

double ESRandomField::rms(const double &x, const double &y, const double &z) const {
    const double r_cyl{std::sqrt(x * x + y * y) - std::sqrt(observer[0] * observer[0] + observer[1] * observer[1])};
    const double zz{std::fabs(z) - std::fabs(observer[2])};
    return b0 * std::sqrt(std::exp(-r_cyl / r0) * std::exp(-zz / z0));
}

}
