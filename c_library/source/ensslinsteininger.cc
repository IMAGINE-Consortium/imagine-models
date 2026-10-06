#include <cmath>

#include "ImagineModelsRandom/EnsslinSteininger.h"

namespace imagine {

double ESRandomField::spectrum(const double &abs_k) const {
  return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

double ESRandomField::rms(const double &x, const double &y, const double &z) const {
  const double r_cyl{std::sqrt(x * x + y * y)};
  const double zz{std::fabs(z)};
  return std::exp(-r_cyl / r0) * std::exp(-zz / z0);
}

}
