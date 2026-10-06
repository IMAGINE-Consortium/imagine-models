#include "ImagineModelsRandom/GaussianScalar.h"

namespace imagine {

double GaussianScalarField::spectrum(const double &abs_k) const {
  return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

}
