#include "ImagineModels/Pshirkov.h"
#include "ImagineModels/units.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace imagine {

void PshirkovMagneticField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown Pshirkov model '" + model + "'.");
    active_model = model;
    parameters = PshirkovParameters<double>{};
    if (model == "ASS") {
        parameters.pitch = -5.;
        parameters.B0_Hs = 2.;
    }
}

template <typename T>
Vec3<T> PshirkovMagneticField::field(const double &x, const double &y, const double &z,
                                     const PshirkovParameters<T> &p) const {
    const double r = std::sqrt(x * x + y * y);
    const double phi = atan2(y, x);
    Vec3<T> b{{0.0, 0.0, 0.0}};

    if ((x == 0.) && (y == 0.)) {
        return b;
    }

    auto pitch = p.pitch * M_PI / 180;

    auto cos_pitch = cos(pitch);
    auto sin_pitch = sin(pitch);
    auto PHI = cos_pitch / sin_pitch * log(1. + p.d / p.R_sun) - M_PI / 2;
    auto cos_PHI = cos(PHI);

    // disk field
    if (useDisk) {
        auto theta = M_PI - phi; // PT11 azimuth convention
        double cos_theta = -x / r;
        double sin_theta = y / r;

        b[0] = -sin_pitch * cos_theta + cos_pitch * sin_theta;
        b[1] = sin_pitch * sin_theta + cos_pitch * cos_theta;
        auto bMag = cos(theta - cos_pitch / sin_pitch * log(r / p.R_sun) + PHI); // eq. 3 / 4
        if ((active_model == "ASS") and (bMag < 0))
            bMag *= -1.;
        bMag *= p.B0_D * p.R_sun / std::max(r, R_c) / cos_PHI * exp(-fabs(z) / p.z0_D); // eq. 5, eq. 4
        b[0] *= bMag;
        b[1] *= bMag;
        b[2] *= bMag;
    }

    // halo field
    if (useHalo) {
        auto bMag = (z > 0 ? p.B0_Hn : -p.B0_Hs);
        auto z1 = (fabs(z) < p.z0_H ? p.z11_H : p.z12_H);
        bMag *= r / p.R0_H * exp(1 - r / p.R0_H) / (1 + pow((fabs(z) - p.z0_H) / z1, 2.));
        // eq. 8
        b[0] += -y / r * bMag;
        b[1] += x / r * bMag;
    }

    return b;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(PshirkovMagneticField)

}
