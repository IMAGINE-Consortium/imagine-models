#include <cmath>
#include <iostream>
#include <random>

#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

void RandomVectorField::unit_random_numbers(std::array<FFTWWorkspace*, 3> ws, const RegularGrid &grid, const int seed) const {
  auto gen_int = std::mt19937(seed);
  std::uniform_int_distribution<int> uni(0, 1215752192);
  const double norm = 1. / std::sqrt(3. * ws[0]->size());
  for (int i = 0; i < 3; ++i) {
    seed_complex_random_numbers(ws[i]->complex(), grid.shape, grid.increment, uni(gen_int));
    ws[i]->backward();
    double* val = ws[i]->real();
    for (std::size_t s = 0; s < ws[i]->padded_size(); ++s)
      val[s] *= norm;
  }
}

VectorGridData RandomVectorField::random_numbers(const RegularGrid &grid, const int seed) const {
  FFTWWorkspace w0(grid.shape), w1(grid.shape), w2(grid.shape);
  std::array<FFTWWorkspace*, 3> ws{&w0, &w1, &w2};
  unit_random_numbers(ws, grid, seed);
  VectorGridData out(grid.shape);
  for (int i = 0; i < 3; ++i)
    ws[i]->copy_unpadded(out.component(i));
  return out;
}

VectorGridData RandomVectorField::sample(const RegularGrid &grid, const int seed) const {
  FFTWWorkspace w0(grid.shape), w1(grid.shape), w2(grid.shape);
  std::array<FFTWWorkspace*, 3> ws{&w0, &w1, &w2};
  _sample(ws, grid, seed);
  VectorGridData out(grid.shape);
  for (int i = 0; i < 3; ++i)
    ws[i]->copy_unpadded(out.component(i));
  return out;
}

void RandomVectorField::_sample(std::array<FFTWWorkspace*, 3> ws, const RegularGrid &grid, const int seed) const {

  const std::array<int, 3> &shp = grid.shape;
  const std::array<double, 3> &inc = grid.increment;
  std::array<double*, 3> val{ws[0]->real(), ws[1]->real(), ws[2]->real()};

  unit_random_numbers(ws, grid, seed);

  auto apply_profile = [&](std::array<double, 3> &b_rand_val, const double xx, const double yy, const double zz) {

    double sp = rms(xx, yy, zz);
    b_rand_val[0] *= sp;
    b_rand_val[1] *= sp;
    b_rand_val[2] *= sp;

    if (apply_anisotropy) {
      Vec3<double> b_reg_val = anisotropy_direction(xx, yy, zz);
      double b_reg_x = b_reg_val[0];
      double b_reg_y = b_reg_val[1];
      double b_reg_z = b_reg_val[2];

      double b_reg_length = std::sqrt(std::pow(b_reg_x, 2) + std::pow(b_reg_y, 2) + std::pow(b_reg_z, 2));

      if (b_reg_length > 1e-10) { // non zero regular field, -> prefered anisotropy

        b_reg_x /= b_reg_length;
        b_reg_y /= b_reg_length;
        b_reg_z /= b_reg_length;
        const double rho2 = anisotropy_rho * anisotropy_rho;
        const double rhonorm = 1. / std::sqrt(0.33333333 * rho2 + 0.66666667 / rho2);
        double reg_dot_rand  = b_reg_x*b_rand_val[0] + b_reg_y*b_rand_val[1] + b_reg_z*b_rand_val[2];

        for (int ii=0; ii==3; ++ii) {
          double b_rand_par = b_rand_val[ii] / reg_dot_rand;
          double b_rand_perp = b_rand_val[ii]  - b_rand_par;
          b_rand_val[ii] = (b_rand_par * anisotropy_rho + b_rand_perp / anisotropy_rho) * rhonorm;
        }
      }
    }
    return b_rand_val;
  };
  for_each_point(RegularGrid(ws[0]->padded_shape(), grid.reference_point, inc), [&](std::size_t idx, double xx, double yy, double zz) {
    std::array<double, 3> b{val[0][idx], val[1][idx], val[2][idx]};
    std::array<double, 3> eval = apply_profile(b, xx, yy, zz);
    val[0][idx] = eval[0];
    val[1][idx] = eval[1];
    val[2][idx] = eval[2];
  });

  if (clean_divergence) {
    for (int i = 0; i < 3; ++i)
      ws[i]->forward();
    divergence_cleaner(ws[0]->complex(), ws[1]->complex(), ws[2]->complex(), shp, inc);
    const double norm = 1. / double(ws[0]->size());
    for (int i = 0; i < 3; ++i) {
      ws[i]->backward();
      for (std::size_t s = 0; s < ws[i]->padded_size(); ++s)
        (val[i])[s] *= norm;
    }
  }
}

// this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
// original author: https://github.com/gioacchinowang
void RandomVectorField::divergence_cleaner(fftw_complex* bx, fftw_complex* by, fftw_complex* bz,  const std::array<int, 3> &shp, const std::array<double, 3> &inc) const {
    double lx = shp[0]*inc[0];
    double ly = shp[1]*inc[1];
    double lz = shp[2]*inc[2];
  
    #ifdef _OPENMP
      #pragma omp parallel for schedule(static)
    #endif
      for (int i = 0; i < shp[0]; ++i) {
        double kx = i / lx;
        if (i >= (shp[0] + 1) / 2)
          kx -= 1. /  inc[0];
          // it's faster to calculate indices manually
        const int idx_lv1 = i * shp[1] * shp[2];
        for (int j = 0; j < shp[1]; ++j) {
          double ky = j /  ly;
          if (j >= (shp[1] + 1) / 2)
            ky -= 1. /  inc[1];
          const int idx_lv2 = idx_lv1 + j * shp[2];
          for (int l = 0; l < (int)shp[2]/2 + 1; ++l) {
            // 0th term is fixed to zero in allocation
            if (i == 0 and j == 0 and l == 0)
              continue;
            double kz = l /  lz;
            const int idx = idx_lv2 + l;
            double k_length = 0;
            double b_length = 0;
            double b_dot_k = 0;
            std::array<double, 3> k{kx, ky, kz};
            std::array<double, 3> b{(*bx)[idx], (*by)[idx], (*bz)[idx]};
            b_length = static_cast<double>(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
            k_length = static_cast<double>(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);

            if (k_length == 0 or b_length == 0) {
              continue;
              }
            k_length = std::sqrt(k_length);
            b_dot_k = (b[0]*k[0] + b[1]*k[1] + b[2]*k[2]);

            const double bk_over_k = b_dot_k / k_length;
            // multiply \sqrt(3) for preserving spectral power statistically
            (*bx)[idx] = 1.73205081 * ((*bx)[idx] - k[0] * bk_over_k);
            (*by)[idx] = 1.73205081 * ((*by)[idx] - k[1] * bk_over_k);
            (*bz)[idx] = 1.73205081 * ((*bz)[idx] - k[2] * bk_over_k);
          } // l
        } // j
      } // i
    }

}
