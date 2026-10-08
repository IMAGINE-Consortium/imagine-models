// Reference: Page et al. 2007, arXiv:astro-ph/0603450; corrected form: Jansson et al. 2009, arXiv:0905.2228 (Sec. 5.3.6)
// Based on: hammurabi v3.01 (old hammurabi), GPL-3.0
// Deviations:
// - psi0 = 27 deg and sin(psi) on the radial component, as corrected in Jansson et al. 2009
// - the paper gives no amplitude, b_b0 = 6 muG is a default
// - anti (field reversed for z > 0) from hammurabi, not in the paper

#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define WMAP_PARAMETERS(X)      \
    X(b_Rsun, 8.)  /* kpc */    \
    X(b_b0, 6.)    /* muG */    \
    X(b_z0, 1.)    /* kpc */    \
    X(b_r0, 8.)    /* kpc */    \
    X(b_psi0, 27)  /* degree */ \
    X(b_psi1, 0.9) /* degree */ \
    X(b_xsi0, 25)  /* degree */

IMAGINE_PARAMETERS(WMAPParameters, WMAP_PARAMETERS)

class WMAPMagneticField : public RegularVectorModel<WMAPMagneticField, WMAPParameters> {
public:
    double b_r_max = 20.; // kpc
    double b_r_min = 3.;  // kpc

    bool anti = false;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const WMAPParameters<T> &p) const;
};

}
