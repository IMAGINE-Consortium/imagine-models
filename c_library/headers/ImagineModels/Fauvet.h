#ifndef FAUVET_H
#define FAUVET_H

#include <functional>
#include <cmath>

#include "ImagineModels/RegularModel.h"

namespace imagine {

// Fauvet magnetic field

#define FAUVET_PARAMETERS(X)     \
    X(b_b0, 2.1) /* muG */       \
    X(b_RB, 8.5) /* kpc */       \
    X(b_Rsun, 8.) /* kpc */      \
    X(b_z0, 1.) /* kpc */        \
    X(b_r0, 7.1) /* kpc */       \
    X(b_p, -30.) /* degree */    \
    X(b_chi0, 22.4) /* degree */

IMAGINE_PARAMETERS(FauvetParameters, FAUVET_PARAMETERS)

class FauvetMagneticField : public RegularVectorModel<FauvetMagneticField, FauvetParameters>
{
public:
    double b_r_max = 20.; // kpc
    double b_r_min = 3.;  // kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const FauvetParameters<T> &p) const;
};

}

#endif