#pragma once

#include <array>

#include <fftw3.h>

#include "ImagineModels/Grid.h"
#include "ImagineModels/types.h"
#include "ImagineModelsRandom/fftw.h"

namespace imagine {

class RandomField {
protected:
    void seed_complex_random_numbers(fftw_complex *vec, const std::array<int, 3> &shp, const std::array<double, 3> &inc,
                                     const int seed) const;

    double mode_power(const double &abs_k) const;

    double simple_spectrum(const double &abs_k, const double &k0, const double &s) const;

public:
    virtual ~RandomField() = default;

    bool apply_spectrum = true;

    virtual double spectrum(const double &abs_k) const = 0;

    virtual double rms(const double &x, const double &y, const double &z) const = 0;

    double variance(const double &x, const double &y, const double &z) const {
        const double r = rms(x, y, z);
        return r * r;
    }

    ScalarGridData evaluate_rms(const Grid &grid) const;
};

}
