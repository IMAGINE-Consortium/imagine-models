#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "Field.h"
#include "RandomScalarField.h"

class LogNormalScalarField : public RandomScalarField {
  protected:
    bool DEBUG = false;
  public:
    using RandomScalarField :: RandomScalarField;

    double log_mean = 0;
    double spectral_amplitude = 1.;
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    void _on_grid(double* val, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) override;

    double calculate_fourier_sigma(const double &abs_k, const double &dk) const override;

    double spatial_profile(const double &x, const double &y, const double &z) const override {
        return 1.;
    }; 

};
