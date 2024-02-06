#include <cmath>
#include <cassert>
#include <iostream>
#include "hamunits.h"
#include "LogNormal.h"

void LogNormalScalarField::_on_grid(double* val, const std::array<int, 3> &shp, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed) {

      fftw_complex* val_comp = construct_plans(val, shp);;
        
      seed_complex_random_numbers(val_comp, shp, grid_increment, seed);
      
      fftw_execute(c2r);
      std::array<int, 3> padded_shp = {shp[0],  shp[1],  2*(shp[2]/2 + 1)}; 
      int pad =  padded_shp[2] - shp[2];
      remove_padding(val, shp, pad);
      // normalize, add mean and exponentiate
      int gs = grid_size(shp);
      for (int s = 0; s < gs; ++s)
        val[s] = std::exp(val[s]/std::sqrt(gs) + log_mean);  
}


double LogNormalScalarField::calculate_fourier_sigma(const double &abs_k, const double &dk) const {
  double sigma = simple_spectrum(abs_k, dk, spectral_offset, spectral_slope);
  return sigma;
}