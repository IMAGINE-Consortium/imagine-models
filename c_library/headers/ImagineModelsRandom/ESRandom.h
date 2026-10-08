// Reference: none
// Based on: hammurabiX (brnd_es)
// Deviations:
// - no publication; rms profile as in hammurabiX, b0 * sqrt(exp(-(r - r_obs)/r0) exp(-(|z| - |z_obs|)/z0))

#pragma once

#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

class ESRandomField : public RandomVectorField {
public:
    double b0 = 0.8; // muG
    double r0 = 8.;
    double z0 = 1.;
    std::array<double, 3> observer{-8.3, 0., 0.006};
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;
};

}
