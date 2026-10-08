#ifndef HELIX_H
#define HELIX_H

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define HELIX_PARAMETERS(X) \
    X(ampx, 1.)             \
    X(ampy, 1.)             \
    X(ampz, 1.)

IMAGINE_PARAMETERS(HelixParameters, HELIX_PARAMETERS)

class HelixMagneticField : public RegularVectorModel<HelixMagneticField, HelixParameters> {
public:
    // non_differentiable parameters
    double rmax = 20.;
    double rmin = 1.;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const HelixParameters<T> &p) const;
};

}

#endif
