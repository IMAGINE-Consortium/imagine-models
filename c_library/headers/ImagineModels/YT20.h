// Reference: Yamasaki & Totani 2020, ApJ 888, 105, arXiv:1909.00849 (eqs. 2-4)
// Deviations:
// - physical constants at full precision (the authors' script and pygedm use three digits), Upsilon = 2.61
// - density set to zero beyond r_vir, where the paper ends the DM integration

#pragma once

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define YT20_PARAMETERS(X)         \
    X(n0_disk, 7.4e-3)             \
    X(R0_disk, 4.9)                \
    X(z0_disk, 2.4)                \
    X(Z_halo, 0.3) /* Z_sun */     \
    X(T_halo, 0.3) /* keV */       \
    X(M_b, 1.2)    /* 1e11 Msun */ \
    X(M_vir, 1.)   /* 1e12 Msun */ \
    X(r_vir, 260.)                 \
    X(c_NFW, 12.)

IMAGINE_PARAMETERS(YT20Parameters, YT20_PARAMETERS)

class YT20 : public RegularScalarModel<YT20, YT20Parameters> {
public:
    double mu = 0.62;
    double mu_e = 1.18;

    template <typename T> T field(const double &x, const double &y, const double &z, const YT20Parameters<T> &p) const;
};

}
