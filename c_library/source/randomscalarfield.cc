#include <cmath>
#include <iostream>

#include "ImagineModelsRandom/RandomScalarField.h"

namespace imagine {

ScalarGridData RandomScalarField::sample(const RegularGrid &grid, const int seed) const {
  FFTWWorkspace ws(grid.shape);
  ScalarGridData out(grid.shape);
  _sample(ws, grid, seed, out);
  return out;
}

ScalarGridData RandomScalarField::random_numbers(const RegularGrid &grid, const int seed) const {
  FFTWWorkspace ws(grid.shape);
  ScalarGridData out(grid.shape);
  double* val = ws.real();
  int gs = ws.size();
  double sqrt_gs = std::sqrt(gs);

  seed_complex_random_numbers(ws.complex(), grid.shape, grid.increment, seed);
  ws.backward();
  for (std::size_t s = 0; s < ws.padded_size(); ++s)  {
    val[s] /= sqrt_gs;
  }
  ws.copy_unpadded(out.component(0));
  return out;
}

void RandomScalarField::_sample(FFTWWorkspace &ws, const RegularGrid &grid, const int seed, ScalarGridData &out) const {
  double* val = ws.real();
  int gs = ws.size();
  double sqrt_gs = std::sqrt(gs);
  // Step 1: draw random numbers with variance 1, possibly correlated

  seed_complex_random_numbers(ws.complex(), grid.shape, grid.increment, seed);
  ws.backward();

  // Step 2: apply spatial amplitude, possibly introduce anisotropy depending on regular field.
  if (!no_profile) {
    for_each_point(RegularGrid(ws.padded_shape(), grid.reference_point, grid.increment), [&](std::size_t idx, double xx, double yy, double zz) {
      double sp = spatial_profile(xx, yy, zz);
      // apply profile
      val[idx] *= sp;
    });
  }

  for (std::size_t s = 0; s < ws.padded_size(); ++s)  {
    val[s] /= sqrt_gs;
  }
  ws.copy_unpadded(out.component(0));
}

}
