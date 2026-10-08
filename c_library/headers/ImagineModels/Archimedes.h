// Reference: Jokipii, Levy & Hubbard 1977, ApJ 213, 861
// Based on: CRPropa (ArchimedeanSpiralField), GPL-3.0
// Deviations:
// - no fitted model; dimensionless parameters as in CRPropa (R_0 in kpc, Omega / v_w in 1/kpc, B_0 the radial field at R_0)

#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define ARCHIMEDES_PARAMETERS(X) \
    X(R_0, 3)                    \
    X(Omega, 1.)                 \
    X(v_w, 0.4)                  \
    X(B_0, 1.)

IMAGINE_PARAMETERS(ArchimedeanParameters, ARCHIMEDES_PARAMETERS)

class ArchimedeanMagneticField : public RegularVectorModel<ArchimedeanMagneticField, ArchimedeanParameters> {
public:
    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const ArchimedeanParameters<T> &p) const;
};

}
