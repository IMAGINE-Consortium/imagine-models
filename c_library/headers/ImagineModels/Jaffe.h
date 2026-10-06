#ifndef JAFFE_H
#define JAFFE_H

#include <cmath>
#include <vector>
#include <iostream>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define JAFFE_PARAMETERS(X)                               \
    X(disk_amp, 0.167) /* disk amplitude, microG */       \
    X(disk_z0, 0.1) /* disk height scale, kpc */          \
    X(halo_amp, 1.38) /* halo amplitude, microG */        \
    X(halo_z0, 3.0) /* halo height scale, kpc */          \
    X(r_inner, 0.5) /* inner R scale, kpc */              \
    X(r_scale, 20.) /* R scale, kpc */                    \
    X(r_peak, 0.) /* R peak, kpc */                       \
    X(ring_amp, 0.023) /* ring field amplitude, microG */ \
    X(ring_r, 5.0) /* ring radius, kpc */                 \
    X(bar_amp, 0.023) /* bar field amplitude, microG */   \
    X(bar_a, 5.0) /* major scale, kpc */                  \
    X(bar_b, 3.0) /* minor scale, kpc */                  \
    X(bar_phi0, 45.0) /* bar major direction */           \
    X(arm_r0, 7.1) /* arm ref radius, kpc */              \
    X(arm_z0, 0.1) /* arm heigth scale, kpc */            \
    X(arm_phi1, 70) /* arm ref angles, deg */             \
    X(arm_phi2, 160)                                      \
    X(arm_phi3, 250)                                      \
    X(arm_phi4, 340)                                      \
    X(arm_amp1, 2) /* arm field amplitudes, microG */     \
    X(arm_amp2, 0.133)                                    \
    X(arm_amp3, -3.78)                                    \
    X(arm_amp4, 0.32)                                     \
    X(arm_pitch, 11.5) /* pitch angle, deg */             \
    X(comp_c, 0.5) /* compress factor */                  \
    X(comp_d, 0.3) /* arm cross-sec scale, kpc */         \
    X(comp_r, 12) /* radial cutoff scale, kpc */          \
    X(comp_p, 3) /* cutoff power */

IMAGINE_PARAMETERS(JaffeParameters, JAFFE_PARAMETERS)

class JaffeMagneticField : public RegularVectorModel<JaffeMagneticField, JaffeParameters>
{
public:
    bool quadruple = false; // quadruple pattern in halo
    bool bss = false;       // bi-symmetric

    bool ring = false;       // molecular ring
    bool bar = true;         // elliptical bar
                             // either ring or bar!

    int arm_num = 4;      // # of spiral arms

    template <typename T>
    Vec3<T> orientation(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    T radial_scaling(const double &x, const double &y, const JaffeParameters<T> &p) const;

    template <typename T>
    std::vector<T> arm_compress(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    std::vector<T> arm_compress_dust(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    std::vector<T> dist2arm(const double &x, const double &y, const JaffeParameters<T> &p) const;

    template <typename T>
    T arm_scaling(const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    T disk_scaling(const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    T halo_scaling(const double &z, const JaffeParameters<T> &p) const;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const;
};

}

#endif