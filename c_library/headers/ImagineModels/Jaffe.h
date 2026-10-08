// Reference: Jaffe et al. 2010, arXiv:0907.3994; Jaffe13 variant: Jaffe et al. 2013, arXiv:1302.0143 (Table A1)
// Based on: hammurabiX (breg_jaffe), GPL-3.0; Jaffe13 variant: hammurabi v3.01 (field_bb), GPL-3.0
// Deviations:
// - 3D form and default parameters from the hammurabiX template, not from a publication (the 2010 model is 2D, with R1 = 3 kpc and an arm cutoff at 15 kpc)
// - Jaffe13: Table A1 azimuths read as clockwise, so the arms are at 350, 260, 170, 80 deg (amplitudes 3, 0.5, -4, 1.2 muG); this matches the NE2001 arms and puts the -4 muG arm at the Sagittarius-Carina arm
// - Jaffe13: arm height h_c = 2 kpc from Table A1 (Planck XLII quotes 0.5 kpc as the original value)
// - Jaffe13: field directions, arm distances and the reversal inside a negative ring as in hammurabi v3.01; no reversal above the disk (introduced only for Jaffe13b)

#pragma once

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define JAFFE_PARAMETERS(X)                            \
    X(disk_amp, 0.167) /* disk amplitude, muG */       \
    X(disk_z0, 0.1)    /* disk height scale, kpc */    \
    X(halo_amp, 1.38)  /* halo amplitude, muG */       \
    X(halo_z0, 3.0)    /* halo height scale, kpc */    \
    X(r_inner, 0.5)    /* inner R scale, kpc */        \
    X(r_scale, 20.)    /* R scale, kpc */              \
    X(r_peak, 0.)      /* R peak, kpc */               \
    X(ring_amp, 0.023) /* ring field amplitude, muG */ \
    X(ring_r, 5.0)     /* ring radius, kpc */          \
    X(bar_amp, 0.023)  /* bar field amplitude, muG */  \
    X(bar_a, 5.0)      /* major scale, kpc */          \
    X(bar_b, 3.0)      /* minor scale, kpc */          \
    X(bar_phi0, 45.0)  /* bar major direction */       \
    X(arm_r0, 7.1)     /* arm ref radius, kpc */       \
    X(arm_z0, 0.1)     /* arm height scale, kpc */     \
    X(arm_phi1, 70)    /* arm ref angles, deg */       \
    X(arm_phi2, 160)                                   \
    X(arm_phi3, 250)                                   \
    X(arm_phi4, 340)                                   \
    X(arm_amp1, 2) /* arm field amplitudes, muG */     \
    X(arm_amp2, 0.133)                                 \
    X(arm_amp3, -3.78)                                 \
    X(arm_amp4, 0.32)                                  \
    X(arm_pitch, 11.5) /* pitch angle, deg */          \
    X(comp_c, 0.5)     /* compress factor */           \
    X(comp_d, 0.3)     /* arm cross-sec scale, kpc */  \
    X(comp_r, 12)      /* radial cutoff scale, kpc */  \
    X(comp_p, 3)       /* cutoff power */

IMAGINE_PARAMETERS(JaffeParameters, JAFFE_PARAMETERS)

class JaffeMagneticField : public RegularVectorModel<JaffeMagneticField, JaffeParameters> {
public:
    const std::array<std::string, 2> available_models{"hammurabiX", "Jaffe13"};
    explicit JaffeMagneticField(const std::string &model = "hammurabiX") { set_model(model); }
    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    bool quadruple = false; // quadruple pattern in halo
    bool bss = false;       // bi-symmetric

    bool ring = false; // molecular ring
    bool bar = true;   // elliptical bar, replaces ring

    int arm_num = 4; // # of spiral arms

    bool hammurabi_v3 = false; // hammurabi v3.01 conventions
    double r_max = 0.;         // kpc, 0 means none

    template <typename T>
    Vec3<T> orientation(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const;

    template <typename T> T radial_scaling(const double &x, const double &y, const JaffeParameters<T> &p) const;

    template <typename T>
    std::vector<T> arm_compress(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    std::vector<T> arm_compress_dust(const double &x, const double &y, const double &z,
                                     const JaffeParameters<T> &p) const;

    template <typename T> std::vector<T> dist2arm(const double &x, const double &y, const JaffeParameters<T> &p) const;

    template <typename T> T arm_scaling(const double &z, const JaffeParameters<T> &p) const;

    template <typename T> T disk_scaling(const double &z, const JaffeParameters<T> &p) const;

    template <typename T> T halo_scaling(const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const;

private:
    std::string active_model = "hammurabiX";
};

}
