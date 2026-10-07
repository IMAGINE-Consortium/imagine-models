#include <cmath>
#include <iostream>
#include <random>

#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

void RandomVectorField::unit_random_numbers(std::array<FFTWWorkspace*, 3> ws, const RegularGrid &grid, const int seed) const {
  auto gen_int = std::mt19937(seed);
  const double norm = 1. / std::sqrt(3. * ws[0]->size());
  for (int i = 0; i < 3; ++i) {
    seed_complex_random_numbers(ws[i]->complex(), grid.shape, grid.increment, static_cast<int>(gen_int() >> 1));
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
      Vec3<double> e = anisotropy_direction(xx, yy, zz);
      const double e_length = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
      if (e_length > 1e-10) {
        for (double &c : e)
          c /= e_length;
        const double rho2 = anisotropy_rho * anisotropy_rho;
        const double rhonorm = 1. / std::sqrt(rho2 / 3. + 2. / (3. * rho2));
        const double b_dot_e = b_rand_val[0] * e[0] + b_rand_val[1] * e[1] + b_rand_val[2] * e[2];
        for (int ii = 0; ii < 3; ++ii) {
          const double b_par = b_dot_e * e[ii];
          const double b_perp = b_rand_val[ii] - b_par;
          b_rand_val[ii] = (b_par * anisotropy_rho + b_perp / anisotropy_rho) * rhonorm;
        }
      }
    }
    return b_rand_val;
  };
  const std::size_t nz = shp[2];
  const std::size_t padded_nz = ws[0]->padded_shape()[2];
  for_each_point(grid, [&](std::size_t point, double xx, double yy, double zz) {
    const std::size_t idx = point / nz * padded_nz + point % nz;
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

void RandomVectorField::divergence_cleaner(fftw_complex* bx, fftw_complex* by, fftw_complex* bz,  const std::array<int, 3> &shp, const std::array<double, 3> &inc) const {
  const double lx = shp[0] * inc[0];
  const double ly = shp[1] * inc[1];
  const double lz = shp[2] * inc[2];
  const int size_z = shp[2] / 2 + 1;
  const double power_correction = std::sqrt(1.5);

  for (int i = 0; i < shp[0]; ++i) {
    const double kx = (i > shp[0] / 2. ? i - shp[0] : i) / lx;
    for (int j = 0; j < shp[1]; ++j) {
      const double ky = (j > shp[1] / 2. ? j - shp[1] : j) / ly;
      for (int l = 0; l < size_z; ++l) {
        const double kz = l / lz;
        const int idx = (i * shp[1] + j) * size_z + l;
        const double k2 = kx * kx + ky * ky + kz * kz;
        const bool nyquist = (i == shp[0] / 2. or j == shp[1] / 2. or l == shp[2] / 2.);
        for (int part = 0; part < 2; ++part) {
          if (k2 == 0. or nyquist) {
            bx[idx][part] = by[idx][part] = bz[idx][part] = 0.;
            continue;
          }
          const double k_dot_b = (kx * bx[idx][part] + ky * by[idx][part] + kz * bz[idx][part]) / k2;
          bx[idx][part] = power_correction * (bx[idx][part] - kx * k_dot_b);
          by[idx][part] = power_correction * (by[idx][part] - ky * k_dot_b);
          bz[idx][part] = power_correction * (bz[idx][part] - kz * k_dot_b);
        }
      }
    }
  }
}

}
