#include <cmath>
#include <random>
#include <stdexcept>

#include "ImagineModelsRandom/RandomField.h"

namespace imagine {

namespace {

double uniform_53(std::mt19937 &gen) {
    const double a = gen() >> 5;
    const double b = gen() >> 6;
    return (a * 67108864. + b) / 9007199254740992.;
}

class StandardNormal {
    bool has_spare = false;
    double spare = 0.;

public:
    double operator()(std::mt19937 &gen) {
        if (has_spare) {
            has_spare = false;
            return spare;
        }
        const double r = std::sqrt(-2. * std::log(1. - uniform_53(gen)));
        const double phi = 6.283185307179586 * uniform_53(gen);
        spare = r * std::sin(phi);
        has_spare = true;
        return r * std::cos(phi);
    }
};

}

double RandomField::simple_spectrum(const double &abs_k, const double &k0, const double &s) const {
    double pi = 3.141592653589793;
    const double unit = 1. / (4 * pi * abs_k * abs_k);
    return unit / std::pow(abs_k + k0, s);
}

double RandomField::mode_power(const double &abs_k) const {
    if (!apply_spectrum)
        return 1.;
    return abs_k < k_min ? 0. : spectrum(abs_k);
}

void RandomField::seed_complex_random_numbers(fftw_complex *vec, const std::array<int, 3> &shp,
                                              const std::array<double, 3> &inc, const int seed) const {

    const double lx = shp[0] * inc[0];
    const double ly = shp[1] * inc[1];
    const double lz = shp[2] * inc[2];

    const double nyquist_x = shp[0] / 2.;
    const double nyquist_y = shp[1] / 2.;
    const double nyquist_z = shp[2] / 2.;

    const int size_z = shp[2] / 2 + 1;
    const double n = double(shp[0]) * shp[1] * shp[2];

    auto wave_vector_length = [&](int i, int j, int l) {
        const double kx = (i > nyquist_x ? i - shp[0] : i) / lx;
        const double ky = (j > nyquist_y ? j - shp[1] : j) / ly;
        const double kz = l / lz;
        return std::sqrt(kx * kx + ky * ky + kz * kz);
    };

    auto is_nyquist = [&](int i, int j, int l) { return i == nyquist_x or j == nyquist_y or l == nyquist_z; };

    double total_power = 0.;
    for (int i = 0; i < shp[0]; ++i)
        for (int j = 0; j < shp[1]; ++j)
            for (int l = 0; l < size_z; ++l) {
                if ((i == 0 and j == 0 and l == 0) or is_nyquist(i, j, l))
                    continue;
                const double multiplicity = (l == 0 or l == nyquist_z) ? 1. : 2.;
                total_power += multiplicity * mode_power(wave_vector_length(i, j, l));
            }
    if (!(total_power > 0.))
        throw std::invalid_argument("RandomField: no Fourier modes with power on this grid (k_min too large?).");

    auto gen = std::mt19937(seed);
    StandardNormal nd;
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
                const double amplitude =
                    is_nyquist(i, j, l) ? 0. : std::sqrt(n * mode_power(wave_vector_length(i, j, l)) / total_power);

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
                        } else if (i > nyquist_x)
                            cg_idx = (shp[0] - i) * shp[1] * size_z + j * size_z + l;
                    } else if (i_is_zero_or_nyquist) {
                        if (j > nyquist_y)
                            cg_idx = i * shp[1] * size_z + (shp[1] - j) * size_z + l;
                    } else if (i > nyquist_x)
                        cg_idx = (shp[0] - i) * shp[1] * size_z + ((shp[1] - j) % shp[1]) * size_z + l;
                }
                if (cg_idx >= 0) {
                    vec[idx][0] = vec[cg_idx][0];
                    vec[idx][1] = -vec[cg_idx][1];
                } else {
                    vec[idx][0] = amplitude * half * nd(gen);
                    vec[idx][1] = amplitude * half * nd(gen);
                }
            }
        }
    }
}

ScalarGridData RandomField::evaluate_rms(const Grid &grid) const {
    ScalarGridData out(grid_shape(grid));
    for_each_point(grid, [&](std::size_t idx, double x, double y, double z) { out(0, idx) = rms(x, y, z); });
    return out;
}

}
