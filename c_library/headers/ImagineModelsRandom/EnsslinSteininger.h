#ifndef ENSSLINSTEININGER_H
#define ENSSLINSTEININGER_H

#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

class ESRandomField : public RandomVectorField {
  public:
    double r0 = 8.5;
    double z0 = 1.5;
    std::array<double, 3> observer{8.5, 0, 0};
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;
};

}

#endif
