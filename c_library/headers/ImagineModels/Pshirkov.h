#pragma once

#include <array>
#include <cmath>
#include <string>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define PSHIRKOV_PARAMETERS(X)                                                                        \
    X(pitch, -6)   /* pitch angle parameters, deg (paper usues -5 for ASS and -6 for BSS) */          \
    X(d, -0.6)     /* distance to first field reversal, kpc */                                        \
    X(R_sun, 8.5)  /* distance between sun and galactic center, kpc */                                \
    X(z0_D, 1.0)   /* vertical thickness in the galactic disk, kpc */                                 \
    X(B0_D, 2.0)   /* magnetic field scale, muG */                                                    \
    X(z0_H, 1.3)   /* halo vertical position, kpc */                                                  \
    X(R0_H, 8.0)   /* halo radial position, kpc */                                                    \
    X(B0_Hn, 4.0)  /* halo magnetic field scale (north), muG */                                       \
    X(B0_Hs, 4.0)  /* halo magnetic field scale (south), muG (paper usues 2 for ASS and 4 for BSS) */ \
    X(z11_H, 0.25) /* halo vertical thickness towards disc, kpc */                                    \
    X(z12_H, 0.4)  /* halo vertical thickness off the disk, kpc */

IMAGINE_PARAMETERS(PshirkovParameters, PSHIRKOV_PARAMETERS)

class PshirkovMagneticField : public RegularVectorModel<PshirkovMagneticField, PshirkovParameters> {
public:
    const std::array<std::string, 2> available_models{"ASS", "BSS"};

    explicit PshirkovMagneticField(const std::string &model = "BSS") { set_model(model); }

    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    bool useDisk = true; // switch for disk field
    bool useHalo = true; // switch for halo field

    // disk parameters
    double R_c = 5.0; // radius of central region, kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const PshirkovParameters<T> &p) const;

private:
    std::string active_model = "BSS";
};

}
