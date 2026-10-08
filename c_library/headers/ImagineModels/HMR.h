#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

// Harari, Mollerach, Roulet (HMR) see https://arxiv.org/abs/astro-ph/9906309, implementation of
// https://arxiv.org/pdf/astro-ph/0510444.pdf

#define HMR_PARAMETERS(X)             \
    X(b_Rsun, 8.5)       /* kpc */    \
    X(b_z1, 0.3)         /* kpc */    \
    X(b_z2, 4.)          /* kpc */    \
    X(b_r1, 2.)          /* kpc */    \
    X(b_p, -10)          /* degree */ \
    X(b_epsilon0, 10.55) /* kpc */

IMAGINE_PARAMETERS(HMRParameters, HMR_PARAMETERS)

class HMRMagneticField : public RegularVectorModel<HMRMagneticField, HMRParameters> {
public:
    double b_r_max = 20.; // kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const HMRParameters<T> &p) const;
};

}
