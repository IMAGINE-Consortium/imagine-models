// Reference: Stanev 1997, arXiv:astro-ph/9607086 (bisymmetric model)
// Based on: hammurabi v3.01 (old hammurabi)
// Deviations:
// - eq. 4 used with exp(-|z|/z0) (sign missing in the paper)
// - field cut at cylindrical r = 20 kpc (the paper: 20 kpc in all directions)

#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define STANEV_PARAMETERS(X)         \
    X(b_z01, 1.)        /* kpc */    \
    X(b_z02, 4.)        /* kpc */    \
    X(b_z0_border, 0.5) /* kpc */    \
    X(b_r0, 10.55)      /* kpc */    \
    X(b_p, -10)         /* degree */ \
    X(b_Rsun, 8.5)      /* kpc */    \
    X(b_phi0, M_PI)     /* radians */

IMAGINE_PARAMETERS(StanevBSSParameters, STANEV_PARAMETERS)

class StanevBSSMagneticField : public RegularVectorModel<StanevBSSMagneticField, StanevBSSParameters> {
public:
    double b_r_max = 20.; // kpc
    double b_r_min = 4.;  // kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const StanevBSSParameters<T> &p) const;
};

}
