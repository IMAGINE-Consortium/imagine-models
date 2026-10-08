// Reference: Han et al. 2018, arXiv:1712.01997; XH24 variant: Xu & Han 2024, arXiv:2404.02038
// Based on: XH24 variant: authors' code (GMFcal.c, diskB)
// Deviations:
// - XH24 variant: disk as in the authors' code, not in the papers: R_s(6) = 8.16 kpc, an extra zone 10.5-15 kpc with
// B_s7 = 4.5 muG, disk to 20 kpc

#pragma once

#include <array>
#include <cmath>
#include <string>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define HAN_PARAMETERS(X)             \
    X(B_p, 11) /* pitch angle, deg */ \
    X(A, 5.)                          \
    X(H, 0.4)                         \
    X(B_s1, 4.5)                      \
    X(B_s2, -3.0)                     \
    X(B_s3, 6.3)                      \
    X(B_s4, -4.7)                     \
    X(B_s5, 3.3)                      \
    X(B_s6, -8.7)                     \
    X(B_s7, 0.)

IMAGINE_PARAMETERS(HanParameters, HAN_PARAMETERS)

class HanMagneticField : public RegularVectorModel<HanMagneticField, HanParameters> {
public:
    const std::array<std::string, 2> available_models{"Han2018", "XH24"};
    explicit HanMagneticField(const std::string &model = "Han2018") { set_model(model); }
    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    double R_min = 3.;
    double R_max = 15.;
    std::array<double, 8> R_s{3.0, 4.1, 4.9, 6.1, 7.5, 8.5, 10.5, 15.0};

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const HanParameters<T> &p) const;

private:
    std::string active_model = "Han2018";
};

}
