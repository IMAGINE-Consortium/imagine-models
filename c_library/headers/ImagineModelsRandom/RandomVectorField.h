#ifndef RANDOMVECTORFIELD_H
#define RANDOMVECTORFIELD_H

#include "ImagineModelsRandom/RandomField.h"

namespace imagine {

class RandomVectorField : public RandomField {
protected:
  void _sample(std::array<FFTWWorkspace*, 3> ws, const RegularGrid &grid, const int seed) const;

  void unit_random_numbers(std::array<FFTWWorkspace*, 3> ws, const RegularGrid &grid, const int seed) const;

public:
  bool clean_divergence = true;
  bool apply_anisotropy = true;

  double anisotropy_rho = 1.;

  virtual Vec3<double> anisotropy_direction(const double &x, const double &y, const double &z) const {
    return {0., 0., 0.};
  }

  VectorGridData sample(const RegularGrid &grid, const int seed) const;

  VectorGridData random_numbers(const RegularGrid &grid, const int seed) const;

  void divergence_cleaner(fftw_complex* bx, fftw_complex* by, fftw_complex* bz,  const std::array<int, 3> &shp, const std::array<double, 3> &inc) const;
};

}

#endif /* RANDOMVECTORFIELD_H */
