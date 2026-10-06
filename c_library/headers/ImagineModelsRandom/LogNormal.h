#ifndef LOGNORMAL_H
#define LOGNORMAL_H

#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "ImagineModelsRandom/RandomScalarField.h"

namespace imagine {

class LogNormalScalarField : public RandomScalarField {
  protected:
    bool DEBUG = false;
  public:
    using RandomScalarField :: RandomScalarField;

    double log_mean = 0;
    double spectral_amplitude = 1.;
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    void _sample(FFTWWorkspace &ws, const RegularGrid &grid, const int seed, ScalarGridData &out) const override;

    double calculate_fourier_sigma(const double &abs_k, const double &dk) const override;

    double spatial_profile(const double &x, const double &y, const double &z) const override {
        return 1.;
    }; 

};

}

#endif
