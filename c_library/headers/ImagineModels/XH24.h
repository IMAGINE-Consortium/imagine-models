#ifndef XH24_H
#define XH24_H

#include <cmath>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define XH24_PARAMETERS(X) \
    X(B0, 0.73) /* muG */  \
    X(z0, 3.0)  /* kpc */  \
    X(R0, 7.97) /* kpc */  \
    X(RT, 5.31) /* kpc */

IMAGINE_PARAMETERS(XH24Parameters, XH24_PARAMETERS)

class XH24MagneticField : public RegularVectorModel<XH24MagneticField, XH24Parameters> {
public:
    double r_max = 20.; // kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const XH24Parameters<T> &p) const;
};

}

#endif
