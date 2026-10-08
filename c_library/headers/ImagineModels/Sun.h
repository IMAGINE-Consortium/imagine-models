// Reference: Sun et al. 2008, arXiv:0711.1572 (ASS+RING); halo: Sun & Reich 2010, arXiv:1010.4394
// Based on: hammurabi v3.01 (old hammurabi), GPL-3.0
// Deviations:
// - halo with the parameters of Sun & Reich 2010: bH_B0 = 2 muG, bH_z1a/bH_z1b = 0.2/4 kpc

#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define SUN_PARAMETERS(X)                    \
    X(b_Rsun, 8.5)                           \
    X(b_R0, 10.)                             \
    X(b_B0, 2.)                              \
    X(b_z0, 1.)                              \
    X(b_Rc, 5.)                              \
    X(b_Bc, 2.)                              \
    X(b_p, -12.)                             \
    X(bH_B0, 2.) /* muG, Sun & Reich 2010 */ \
    X(bH_R0, 4.)                             \
    X(bH_z0, 1.5)                            \
    X(bH_z1a, 0.2)                           \
    X(bH_z1b, 4.)

IMAGINE_PARAMETERS(SunParameters, SUN_PARAMETERS)

class SunMagneticField : public RegularVectorModel<SunMagneticField, SunParameters> {
public:
    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const SunParameters<T> &p) const;
};

}
