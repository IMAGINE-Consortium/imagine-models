#include <cmath>

#include "ImagineModelsRandom/LogNormal.h"

namespace imagine {

double LogNormalScalarField::transform(const double &g, const double &x, const double &y, const double &z) const {
    return std::exp(log_mu + log_sigma * g);
}

double LogNormalScalarField::mean(const double &x, const double &y, const double &z) const {
    return std::exp(log_mu + 0.5 * log_sigma * log_sigma);
}

double LogNormalScalarField::rms(const double &x, const double &y, const double &z) const {
    return std::sqrt(std::expm1(log_sigma * log_sigma)) * mean(x, y, z);
}

double LogNormalScalarField::spectrum(const double &abs_k) const {
    return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

}
