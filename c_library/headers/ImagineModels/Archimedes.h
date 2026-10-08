#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

// simple archimdeean sprial, implementation based on CRPropa

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
