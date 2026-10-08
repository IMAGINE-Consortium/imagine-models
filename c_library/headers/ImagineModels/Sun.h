#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

// Sun et al. A&A V.477 2008 ASS+RING model magnetic field

#define SUN_PARAMETERS(X)                                                                      \
    X(b_Rsun, 8.5)                                                                             \
    X(b_R0, 10.)                                                                               \
    X(b_B0, 2.)                                                                                \
    X(b_z0, 1.)                                                                                \
    X(b_Rc, 5.)                                                                                \
    X(b_Bc, 2.)                                                                                \
    X(b_p, -12.)                                                                               \
    X(bH_B0, 2.) /* 10 in original publication, 2 in update https://arxiv.org/abs/1010.4394 */ \
    X(bH_R0, 4.)                                                                               \
    X(bH_z0, 1.5)                                                                              \
    X(bH_z1a, 0.2)                                                                             \
    X(bH_z1b, 4.)

IMAGINE_PARAMETERS(SunParameters, SUN_PARAMETERS)

class SunMagneticField : public RegularVectorModel<SunMagneticField, SunParameters> {
public:
    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const SunParameters<T> &p) const;
};

}
