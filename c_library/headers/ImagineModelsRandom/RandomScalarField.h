#ifndef RANDOMSCALARFIELD_H
#define RANDOMSCALARFIELD_H

#include "ImagineModelsRandom/RandomField.h"

namespace imagine {

class RandomScalarField : public RandomField {
protected:
  virtual void _sample(FFTWWorkspace &ws, const RegularGrid &grid, const int seed, ScalarGridData &out) const;

public:
  ScalarGridData sample(const RegularGrid &grid, const int seed) const;

  ScalarGridData random_numbers(const RegularGrid &grid, const int seed) const;
};

}

#endif /* RANDOMSCALARFIELD_H */
