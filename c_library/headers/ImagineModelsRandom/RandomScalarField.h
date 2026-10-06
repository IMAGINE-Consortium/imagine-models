#ifndef RANDOMSCALARFIELD_H
#define RANDOMSCALARFIELD_H

#include "ImagineModelsRandom/RandomField.h"

namespace imagine {

class RandomScalarField : public RandomField {
protected:
  virtual double transform(const double &g, const double &x, const double &y, const double &z) const {
    return mean(x, y, z) + rms(x, y, z) * g;
  }

public:
  virtual double mean(const double &x, const double &y, const double &z) const { return 0.; }

  ScalarGridData sample(const RegularGrid &grid, const int seed) const;

  ScalarGridData random_numbers(const RegularGrid &grid, const int seed) const;
};

}

#endif /* RANDOMSCALARFIELD_H */
