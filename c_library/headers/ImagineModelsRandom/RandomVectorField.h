#pragma once

#include <cmath>

#include "ImagineModelsRandom/RandomField.h"

namespace imagine {

class RandomVectorField : public RandomField {
protected:
    double combined_rms(const double &a_iso, const double &a_ord) const {
        if (a_ord == 0.)
            return a_iso;
        const double cross = independent_ordered ? 0. : 2. * a_iso * a_ord;
        return std::sqrt(a_iso * a_iso + (cross + a_ord * a_ord) / 3.);
    }

    void _sample(std::array<FFTWWorkspace *, 3> ws, const RegularGrid &grid, const int seed) const;

    void unit_random_numbers(std::array<FFTWWorkspace *, 3> ws, const RegularGrid &grid, const int seed) const;

public:
    bool clean_divergence = true;
    bool apply_anisotropy = true;

    double anisotropy_rho = 1.;
    bool independent_ordered = false; // second field for ordered part

    virtual Vec3<double> anisotropy_direction(const double &x, const double &y, const double &z) const {
        return {0., 0., 0.};
    }

    // isotropic part
    virtual double isotropic_rms(const double &x, const double &y, const double &z) const { return rms(x, y, z); }

    // ordered part along anisotropy_direction
    virtual double ordered_amplitude(const double &x, const double &y, const double &z) const { return 0.; }

    VectorGridData sample(const RegularGrid &grid, const int seed) const;

    VectorGridData random_numbers(const RegularGrid &grid, const int seed) const;

    void divergence_cleaner(fftw_complex *bx, fftw_complex *by, fftw_complex *bz, const std::array<int, 3> &shp,
                            const std::array<double, 3> &inc) const;
};

}
