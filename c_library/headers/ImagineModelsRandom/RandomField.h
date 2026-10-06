#ifndef RANDOMFIELD_H
#define RANDOMFIELD_H

#include <array>

#include <fftw3.h>

#include "ImagineModels/types.h"
#include "ImagineModels/Grid.h"
#include "ImagineModelsRandom/fftw.h"

namespace imagine {

class RandomField {
protected:
  bool no_profile = false;

  void seed_complex_random_numbers(fftw_complex* vec,  const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed) const;

  double simple_spectrum(const double &abs_k, const double &dk, const double &k0, const double &s) const;

  double hammurabi_spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const;

public:
  virtual ~RandomField() = default;

  bool apply_spectrum = true;

  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  virtual double calculate_fourier_sigma(const double &abs_k, const double &dk) const = 0;

  ScalarGridData profile(const Grid &grid) const;
};

}

#endif /* RANDOMFIELD_H */
