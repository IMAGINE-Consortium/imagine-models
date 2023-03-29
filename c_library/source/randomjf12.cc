#include <cmath>
#include <cassert>
#include <iostream>
#include "../headers/hamunits.h"
#include "../headers/RandomJF12.h"


double JF12RandomField::calculate_fourier_sigma(const double &abs_k) const {
  double sigma = simple_spectrum(abs_k, spectral_amplitude, spectral_offset, spectral_slope);
  return sigma;
}

double JF12RandomField::spatial_profile(const double &x, const double &y, const double &z) const {

      const double r{sqrt(x * x + y * y)};
      const double rho{
          sqrt(x*x + y*y + z*z)};
      const double phi{atan2(y, z)};

      const double rc_B[8] = {
          5.1, 6.3,  7.1,  8.3,
          9.8, 11.4, 12.7, 15.5}; // neg x crossings of spiral arms
      const double inc = 11.5; // inclination, in degrees

      const double b_arms[8] = {b0_1, b0_2, b0_3, b0_4, b0_5, b0_6, b0_7, b0_8};

      double scaling_disk = 0.0;
      double scaling_halo = 0.0;

      // boundaries outside which B is zero, not sure if this works?
      if (r > Rmax || rho < rho_GC) {
        return 0.0;
      }
      if (r < 5.) {
        scaling_disk = b0_int;
      } else {
        double r_negx =
            r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi - M_PI));
        if (r_negx > rc_B[7]) {
          r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + M_PI));
        }
        if (r_negx > rc_B[7]) {
          r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + 3 * M_PI));
        }
        for (int i = 7; i >= 0; i--) {
          if (r_negx < rc_B[i] ) {
            scaling_disk = b_arms[i] * (5.) / r;
          }
        } // "region 8,7,6,..,2"
      }

      scaling_disk = scaling_disk * exp(-0.5 * z * z / (z0_spiral * z0_spiral));
      scaling_halo =
          b0_halo * exp(-std::fabs(r / r0_halo)) * exp(-std::fabs(z / z0_halo));

      return (scaling_disk * scaling_disk + scaling_halo * scaling_halo);

};

void JF12RandomField::_on_grid(std::array<double*, 3> val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {
  std::array<fftw_complex*, 3> val_comp = construct_plans(val, shp); 
  auto gen_int = std::mt19937(seed);
  std::uniform_int_distribution<int> uni(0, 100000000000);
  for (int i =0; i<3; ++i) {
    int sub_seed = uni(gen_int); 
    std::cout<< "sub seed " << sub_seed << std::endl;
    draw_random_numbers_complex(val_comp[i], shp, inc, sub_seed);
    fftw_execute(c2r[i]);
  }
  auto apply_profile = [&](double xx, double yy, double zz) {
      int _nx = (int)((xx - rpt[0])/inc[0]); 
      int _ny = (int)((yy - rpt[1])/inc[1]);
      int _nz = (int)((zz - rpt[2])/inc[2]);
      int idx = _nz + shp[2]*(_ny + shp[1]*_nx);
      std::array<double, 3> b_reg_val = regular_base.at_position(xx, yy, zz);
      std::cout<<"b_reg_val: " << b_reg_val[0] << ", " << b_reg_val[1] << ", " << b_reg_val[2] << std::endl;
      double b_reg_length = std::sqrt(std::pow(b_reg_val[0], 2) + std::pow(b_reg_val[1], 2) + std::pow(b_reg_val[2], 2));
      double sp = std::sqrt(spatial_profile(xx, yy, zz));
        // assemble b_Re
      std::cout<<"sp: " << sp << std::endl;
      std::cout<<"val: " << *val[0] << ", " << *val[1] << ", " << *val[2] << std::endl;
      std::array<double, 3> b_rand_val{(val[0])[idx] * sp,
                                       (val[1])[idx] * sp,
                                       (val[2])[idx] * sp};
      std::cout<<"b_rand 1: " << b_rand_val[0] << ", " << b_rand_val[1] << ", " << b_rand_val[2] << std::endl;
      if (b_reg_length < 1e-10) // zero regular field, no prefered anisotropy
        return b_rand_val;
      // impose anisotropy
      b_reg_val[0] /= b_reg_length;
      const double rho2 = anisotropy_rho * anisotropy_rho;
      const double rhonorm = 1. / std::sqrt(0.33333333 * rho2 + 0.66666667 / rho2);
      double reg_dot_rand  = b_reg_val[0]*b_rand_val[0] + b_reg_val[1]*b_rand_val[1] + b_reg_val[2]*b_rand_val[2];

      for (int ii=0; ii==3; ++ii) {
        double b_rand_par = b_reg_val[ii] * reg_dot_rand;
        double b_rand_perp = b_rand_val[ii]  - b_rand_par;
        b_rand_val[ii] = (b_rand_par * anisotropy_rho + b_rand_perp / anisotropy_rho) * rhonorm;
      }
      std::cout<<"b_rand 2: " << b_rand_val[0] << ", " << b_rand_val[1] << ", " << b_rand_val[2] << std::endl;
      return b_rand_val;
    };

  evaluate_function_on_grid(val, shp, rpt, inc, apply_profile);
  
  for (int i =0; i<3; ++i) {
    fftw_execute(r2c[i]);
  }
  
  divergence_cleaner(val_comp[0], val_comp[1], val_comp[2], shp, inc);

  int grid_size = shp[0]*shp[1]*shp[2];
  for (int i =0; i<3; ++i) {
    fftw_execute(c2r[i]);
    for (int s = 0; s < grid_size; ++s)
      (val[i])[s] /= grid_size;  
  }
};