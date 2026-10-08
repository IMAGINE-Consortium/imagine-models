#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

// WMAP magnetic field

#define WMAP_PARAMETERS(X)                                                                              \
    X(b_Rsun, 8.)  /* kpc */                                                                            \
    X(b_b0, 6.)    /* muG  -> not given in original paper? Could also be 3 according to                 \
                      https://www.aanda.org/articles/aa/full_html/2010/14/aa12733-09/aa12733-09.html */ \
    X(b_z0, 1.)    /* kpc */                                                                            \
    X(b_r0, 8.)    /* kpc */                                                                            \
    X(b_psi0, 27)  /* degree */                                                                         \
    X(b_psi1, 0.9) /* degree */                                                                         \
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
