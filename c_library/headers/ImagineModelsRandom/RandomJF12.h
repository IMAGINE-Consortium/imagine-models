#ifndef RANDOMJF12_H
#define RANDOMJF12_H

#include <cmath>

#include "ImagineModelsRandom/RandomVectorField.h"
#include "ImagineModels/RegularJF12.h"

namespace imagine {

class JF12RandomField : public RandomVectorField {
  public:
    double b0_1 = 10.81; // uG
    double b0_2 = 6.96; // uG
    double b0_3 = 9.59; // uG
    double b0_4 = 6.96; // uG
    double b0_5 = 1.96; // uG
    double b0_6 = 16.34; // uG
    double b0_7 = 37.29; // uG
    double b0_8 = 10.35; // uG
    double b0_int = 7.63; // uG
    double z0_spiral = 0.61; // kpc
    double b0_halo = 4.68; // uG
    double r0_halo = 10.97; // kpc
    double z0_halo = 2.84; // kpc
    double Rmax = 20.;
    double rho_GC = 1.;
    const double rc_B[8] = {
          5.1, 6.3,  7.1,  8.3,
          9.8, 11.4, 12.7, 15.5}; // neg x crossings of spiral arms
    const double inc = 11.5; // inclination, in degrees
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    JF12MagneticField regular_base = JF12MagneticField();

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;
    Vec3<double> anisotropy_direction(const double &x, const double &y, const double &z) const override;
};

}

#endif
