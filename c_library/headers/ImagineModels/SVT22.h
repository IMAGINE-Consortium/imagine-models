#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define SVT22_PARAMETERS(X) \
    X(B_val, 3.72)          \
    X(r_cut, 5)             \
    X(z_cut, 6)

IMAGINE_PARAMETERS(SVT22Parameters, SVT22_PARAMETERS)

class SVT22MagneticField : public RegularVectorModel<SVT22MagneticField, SVT22Parameters> {
public:
    bool do_halo = true;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const SVT22Parameters<T> &p) const;
};

}
