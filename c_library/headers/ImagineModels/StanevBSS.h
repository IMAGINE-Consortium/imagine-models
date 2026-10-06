#ifndef STANEVBSS_H
#define STANEVBSS_H


#include <functional>
#include <cmath>

#include "ImagineModels/RegularModel.h"

namespace imagine {


// StanevBSS (HMR) see https://arxiv.org/abs/astro-ph/9607086

#define STANEV_PARAMETERS(X)      \
    X(b_z01, 1.) /* kpc */        \
    X(b_z02, 4.) /* kpc */        \
    X(b_z0_border, 0.5) /* kpc */ \
    X(b_r0, 10.55) /* kpc */      \
    X(b_p, -10) /* degree */      \
    X(b_Rsun, 8.5) /* kpc */      \
    X(b_phi0, M_PI) /* radians */

IMAGINE_PARAMETERS(StanevBSSParameters, STANEV_PARAMETERS)

class StanevBSSMagneticField : public RegularVectorModel<StanevBSSMagneticField, StanevBSSParameters>
{
public:
    double b_r_max = 20.; // kpc
    double b_r_min = 4.; // kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const StanevBSSParameters<T> &p) const;
};

}

#endif
