// Reference: Tinyakov & Tkachev 2002, arXiv:astro-ph/0111305; Kachelriess et al. 2007, arXiv:astro-ph/0510444
// Based on: hammurabi v3.01 (old hammurabi)

#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define TT_PARAMETERS(X)                     \
    X(b_Rsun, 8.5) /* kpc */                 \
    X(b_b0, 1.4)   /* muG */                 \
    X(b_d, -0.5)   /* kpc */                 \
    X(b_z0, 1.5)   /* kpc, h in the paper */ \
    X(b_p, -8)     /* degree */

IMAGINE_PARAMETERS(TTParameters, TT_PARAMETERS)

class TTMagneticField : public RegularVectorModel<TTMagneticField, TTParameters> {
public:
    double b_r_max = 20.; // kpc
    double b_r_min = 4.;  // kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const TTParameters<T> &p) const;
};

}
