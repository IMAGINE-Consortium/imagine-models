#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "RandomVectorField.h"


class ESRandomField : public RandomVectorField {
  protected:
    bool DEBUG = false;
  public:
    using RandomVectorField :: RandomVectorField;

    double r0 = 8.5;
    double z0 = 1.5;
    std::array<double, 3> observer{8.5, 0, 0};

    double spectral_amplitude = 1.; 
    double spectral_offset = 1.; 
    double spectral_slope = 2.;


    //void _on_grid(std::array<double*, 3> val, const std::array<int, 3> &shp, const std::array<double, 3> &zpt, const std::array<double, 3> &inc, const int seed) override;

    double calculate_fourier_sigma(const double &abs_k, const double &dk) const override;

    double spatial_profile(const double &x, const double &y, const double &z) const override; 

};
