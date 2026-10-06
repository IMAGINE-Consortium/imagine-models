#include <cmath>
#include <cassert>
#include <iostream>
#include "ImagineModels/units.h"
#include "ImagineModelsRandom/LogNormal.h"

namespace imagine {

void LogNormalScalarField::_sample(FFTWWorkspace &ws, const RegularGrid &grid, const int seed, ScalarGridData &out) const {

      seed_complex_random_numbers(ws.complex(), grid.shape, grid.increment, seed);
      
      ws.backward();
      double* val = out.component(0);
      ws.copy_unpadded(val);
      // normalize, add mean and exponentiate
      int gs = ws.size();
      for (int s = 0; s < gs; ++s)
        val[s] = std::exp(val[s]/std::sqrt(gs) + log_mean);  
}


double LogNormalScalarField::calculate_fourier_sigma(const double &abs_k, const double &dk) const {
  double sigma = simple_spectrum(abs_k, dk, spectral_offset, spectral_slope);
  return sigma;
}

}
