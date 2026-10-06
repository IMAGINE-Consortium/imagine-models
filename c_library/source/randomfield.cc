#include <cmath>
#include <random>

#include "ImagineModelsRandom/RandomField.h"

namespace imagine {

double RandomField::hammurabi_spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const {
  // this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
  // original author: https://github.com/gioacchinowang
  const double p0 = rms*rms;
  double pi = 3.141592653589793;
  const double unit = 1. / (4 * pi * abs_k * abs_k);   // units fixing, wave vector in 1/kpc units
  // power laws
  const double band1 = double(abs_k < k1);
  const double band2 = double(abs_k > k1) * double(abs_k < k0);
  const double band3 = double(abs_k > k0);
  const double P = band1 * std::pow(k0 / k1, a1) * std::pow(abs_k / k1, 6.0) +
                  band2 / std::pow(abs_k / k0, a1) +
                  band3 / std::pow(abs_k / k0, a0);
  return P * p0 * unit;
  }

double RandomField::simple_spectrum(const double &abs_k, const double &k0, const double &s) const {
  double pi = 3.141592653589793;
  const double unit = 1. / (4 * pi * abs_k * abs_k);
  return unit / std::pow(abs_k + k0, s);
}

double RandomField::mode_power(const double &abs_k) const {
  return apply_spectrum ? spectrum(abs_k) : 1.;
}

void RandomField::seed_complex_random_numbers(fftw_complex* vec,  const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed) const {

  const double lx = shp[0]*inc[0];
  const double ly = shp[1]*inc[1];
  const double lz = shp[2]*inc[2];

  const double nyquist_x = shp[0]/2.;
  const double nyquist_y = shp[1]/2.;
  const double nyquist_z = shp[2]/2.;

  const int size_z = shp[2]/2 + 1;
  const double n = double(shp[0]) * shp[1] * shp[2];

  auto wave_vector_length = [&](int i, int j, int l) {
    const double kx = (i > nyquist_x ? i - shp[0] : i) / lx;
    const double ky = (j > nyquist_y ? j - shp[1] : j) / ly;
    const double kz = l / lz;
    return std::sqrt(kx * kx + ky * ky + kz * kz);
  };

  double total_power = 0.;
  for (int i = 0; i < shp[0]; ++i)
    for (int j = 0; j < shp[1]; ++j)
      for (int l = 0; l < size_z; ++l) {
        if (i == 0 and j == 0 and l == 0)
          continue;
        const double multiplicity = (l == 0 or l == nyquist_z) ? 1. : 2.;
        total_power += multiplicity * mode_power(wave_vector_length(i, j, l));
      }

  auto gen = std::mt19937(seed);
  std::normal_distribution<double> nd{0., 1.};
  const double half = std::sqrt(0.5);

  for (int i = 0; i < shp[0]; ++i) {
    const int idx_lv1 = i * shp[1] * size_z;
    for (int j = 0; j < shp[1]; ++j) {
      const int idx_lv2 = idx_lv1 + j * size_z;
      for (int l = 0; l < size_z; ++l) {
        const int idx = idx_lv2 + l;
        if (l == 0 and j == 0 and i == 0) {
          vec[0][0] = 0.;
          vec[0][1] = 0.;
          continue;
        }
        const double amplitude = std::sqrt(n * mode_power(wave_vector_length(i, j, l)) / total_power);

        bool l_is_zero_or_nyquist = (l == 0 or l == nyquist_z);
        bool j_is_zero_or_nyquist = (j == 0 or j == nyquist_y);
        bool i_is_zero_or_nyquist = (i == 0 or i == nyquist_x);
        int cg_idx = -1;

        if (l_is_zero_or_nyquist) {
          if (j_is_zero_or_nyquist) {
            if (i_is_zero_or_nyquist) {
              vec[idx][0] = amplitude * nd(gen);
              vec[idx][1] = 0.;
              continue;
            }
            else if (i > nyquist_x)
              cg_idx = (shp[0] - i) * shp[1] * size_z + j * size_z + l;
          }
          else if (i_is_zero_or_nyquist) {
            if (j > nyquist_y)
              cg_idx = i * shp[1] * size_z + (shp[1] - j) * size_z + l;
          }
          else if (i > nyquist_x)
            cg_idx = (shp[0] - i) * shp[1] * size_z + ((shp[1] - j) % shp[1]) * size_z + l;
        }
        if (cg_idx >= 0) {
          vec[idx][0] = vec[cg_idx][0];
          vec[idx][1] = - vec[cg_idx][1];
        }
        else {
          vec[idx][0] = amplitude * half * nd(gen);
          vec[idx][1] = amplitude * half * nd(gen);
        }
      }
    }
  }
}

ScalarGridData RandomField::evaluate_rms(const Grid &grid) const {
  ScalarGridData out(grid_shape(grid));
  for_each_point(grid, [&](std::size_t idx, double x, double y, double z)
                 { out(0, idx) = rms(x, y, z); });
  return out;
}

}
