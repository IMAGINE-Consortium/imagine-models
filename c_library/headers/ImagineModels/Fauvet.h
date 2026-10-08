// Reference: Fauvet et al. 2012, arXiv:1201.5742 (Sec. 2.1)
// Based on: hammurabi v3.01 (old hammurabi), GPL-3.0

#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define FAUVET_PARAMETERS(X)     \
    X(b_b0, 2.1)    /* muG */    \
    X(b_RB, 8.5)    /* kpc */    \
    X(b_Rsun, 8.)   /* kpc */    \
    X(b_z0, 1.)     /* kpc */    \
    X(b_r0, 7.1)    /* kpc */    \
    X(b_p, -30.)    /* degree */ \
    X(b_chi0, 22.4) /* degree */

IMAGINE_PARAMETERS(FauvetParameters, FAUVET_PARAMETERS)

class FauvetMagneticField : public RegularVectorModel<FauvetMagneticField, FauvetParameters> {
public:
    double b_r_max = 20.; // kpc
    double b_r_min = 3.;  // kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const FauvetParameters<T> &p) const;
};

}
