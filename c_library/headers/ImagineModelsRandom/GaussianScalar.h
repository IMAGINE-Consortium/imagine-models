#ifndef GAUSSIANSCALAR_H
#define GAUSSIANSCALAR_H

#include "ImagineModelsRandom/RandomScalarField.h"

namespace imagine {

class GaussianScalarField : public RandomScalarField {
  public:
    double mu = 0.;
    double sigma = 1.;
    double spectral_offset = .001;
    double spectral_slope = 2.;

    double spectrum(const double &abs_k) const override;
    double mean(const double &x, const double &y, const double &z) const override { return mu; }
    double rms(const double &x, const double &y, const double &z) const override { return sigma; }
};

}

#endif
