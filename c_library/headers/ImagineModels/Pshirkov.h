// Reference: Pshirkov et al. 2011, arXiv:1103.0814 (ASS and BSS of Table 3)
// Based on: CRPropa (PT11Field), GPL-3.0

#pragma once

#include <array>
#include <cmath>
#include <string>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define PSHIRKOV_PARAMETERS(X)                        \
    X(pitch, -6)   /* pitch angle, deg */             \
    X(d, -0.6)     /* first reversal distance, kpc */ \
    X(R_sun, 8.5)  /* Sun distance, kpc */            \
    X(z0_D, 1.0)   /* disk scale height, kpc */       \
    X(B0_D, 2.0)   /* disk field, muG */              \
    X(z0_H, 1.3)   /* halo height, kpc */             \
    X(R0_H, 8.0)   /* halo radius, kpc */             \
    X(B0_Hn, 4.0)  /* northern halo field, muG */     \
    X(B0_Hs, 4.0)  /* southern halo field, muG */     \
    X(z11_H, 0.25) /* inner halo thickness, kpc */    \
    X(z12_H, 0.4)  /* outer halo thickness, kpc */

IMAGINE_PARAMETERS(PshirkovParameters, PSHIRKOV_PARAMETERS)

class PshirkovMagneticField : public RegularVectorModel<PshirkovMagneticField, PshirkovParameters> {
public:
    const std::array<std::string, 2> available_models{"ASS", "BSS"};

    explicit PshirkovMagneticField(const std::string &model = "BSS") { set_model(model); }

    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    bool useDisk = true; // disk switch
    bool useHalo = true; // halo switch

    double R_c = 5.0; // central region radius, kpc

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const PshirkovParameters<T> &p) const;

private:
    std::string active_model = "BSS";
};

}
