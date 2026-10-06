#ifndef HAN_H
#define HAN_H


#include <cmath>

#include "ImagineModels/RegularModel.h"

namespace imagine {

//J. L. Han et al 2018 ApJS 234 11

#define HAN_PARAMETERS(X) \
    X(B_p, 11) /* pitch angle */ \
    X(A, 5.)              \
    X(H, 0.4)             \
    X(B_s1, 4.5)          \
    X(B_s2, -3.0)         \
    X(B_s3, 6.3)          \
    X(B_s4, -4.7)         \
    X(B_s5, 3.3)          \
    X(B_s6, -8.7)

IMAGINE_PARAMETERS(HanParameters, HAN_PARAMETERS)

class HanMagneticField : public RegularVectorModel<HanMagneticField, HanParameters> {
    public:
        double R_min = 3.;
        double R_max = 15.;
        std::array<double, 7> R_s{3.0, 4.1, 4.9, 6.1, 7.5, 8.5, 10.5};

        template <typename T>
        Vec3<T> evaluate(const double &x, const double &y, const double &z, const HanParameters<T> &p) const;
 };

}

 #endif
