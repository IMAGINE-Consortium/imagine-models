#include <cmath>
#include <cassert>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/EnsslinSteininger.h"


double ESRandomField::calculate_fourier_sigma(const double &abs_k) const {
  double sigma = simple_spectrum(abs_k, spectral_amplitude, spectral_offset, spectral_slope);
  return sigma;
};

double ESRandomField::spatial_profile(const double &x, const double &y, const double &z) const {
      const double r_cyl{std::sqrt(x * x + y * y) -
                        std::fabs(observer[0])};
      const double zz{std::fabs(z) - std::fabs(observer[2])};
  return std::exp(-r_cyl / r0) * std::exp(-zz / z0);
};


void ESRandomField::_on_grid(std::array<double*, 3> val, const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      int gs = grid_size(grid_shape);
      
      std::array<fftw_complex*, 3> val_comp = construct_plans(val, shape); 
        
      for (int i =0; i<3; ++i) {
        draw_random_numbers_complex(val_comp[i], grid_shape, grid_increment, seed);
        fftw_execute(c2r[i]);
      }

      auto multiply_profile = [&](double xx, double yy, double zz) {
          int _nx = (int)((xx - grid_zeropoint[0])/grid_increment[0]); 
          int _ny = (int)((yy - grid_zeropoint[1])/grid_increment[1]);
          int _nz = (int)((zz - grid_zeropoint[2])/grid_increment[2]);
          double sp = spatial_profile(xx, yy, zz);
          int idx = _nz + grid_shape[2]*(_ny + grid_shape[1]*_nx);
          std::array<double, 3> v = {
            (val[0])[idx]*sp, 
            (val[1])[idx]*sp, 
            (val[2])[idx]*sp
            };
          return v;
        };

      evaluate_function_on_grid(val, grid_shape, grid_zeropoint, grid_increment, multiply_profile);

      for (int i =0; i<3; ++i) {
        fftw_execute(r2c[i]);
      }
      
      divergence_cleaner(val_comp[0], val_comp[1], val_comp[2], grid_shape, grid_increment);

      for (int i =0; i<3; ++i) {
        fftw_execute(c2r[i]);
        for (int s = 0; s < gs; ++s)
          (val[i])[s] /= gs;  
      }

};