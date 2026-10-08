#pragma once

#include "ImagineModelsRandom/RandomScalarField.h"

namespace imagine {

class LogNormalScalarField : public RandomScalarField {
protected:
    double transform(const double &g, const double &x, const double &y, const double &z) const override;

public:
    double log_mu = 0.;
    double log_sigma = 1.;
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    double spectrum(const double &abs_k) const override;
    double mean(const double &x, const double &y, const double &z) const override;
    double rms(const double &x, const double &y, const double &z) const override;
};

}
