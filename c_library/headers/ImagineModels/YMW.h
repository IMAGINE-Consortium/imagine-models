#ifndef YMW16_H
#define YMW16_H

#include <cassert>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define YMW16_PARAMETERS(X)                                                             \
    X(r0, 8.3)             /* kpc, Galactic earth position */                           \
    X(t1_ad, 2.5)          /* kpc scale length of cutoff */                             \
    X(t1_bd, 15.)          /* kpc radius of begin of cutoff */                          \
    X(t1_n1, 0.01132)      /* cm^{-3}, normalization fitted by YMW */                   \
    X(t1_h1, 1.673)        /* kpc, scale-height, fitted by YMW */                       \
    X(t2_a2, 1.2)          /* kpc, scale length */                                      \
    X(t2_b2, 4.)           /* kpc, molecular ring central radius */                     \
    X(t2_n2, 0.404)        /* cm^{-3}, fitted by YMW */                                 \
    X(t2_k2, 1.54)         /* unitless, rescaling of scale height, fitted by YMW */     \
    X(t3_b2s, 4.)          /* kpc, is actually t2_b2 in paper */                        \
    X(t3_ka, 5.015)        /* spiral arm scale factor */                                \
    X(t3_aa, 11.680)       /* kpc, scale_length */                                      \
    X(t3_ncn, 2.4)         /* unitless, Carina relative over density, fitted by YMW */  \
    X(t3_wcn, 8.2)         /* // degree, theta scaling, fitted by YMW */                \
    X(t3_thetacn, 109.)    /* degree, theta correction, fitted by YMW */                \
    X(t3_nsg, 0.626)       /* unitless, Carina relative under density, fitted by YMW */ \
    X(t3_wsg, 20)          /* degree, theta scaling, fitted by YMW */                   \
    X(t3_thetasg, 75.8)    /* degree, theta correction, fitted by YMW */                \
    X(t4_ngc, 6.2)         /* cm^{-3}, normalization, fitted by YMW */                  \
    X(t4_agc, 0.160)       /* kpc, scale length, fixed by YMW based on CO */            \
    X(t4_hgc, 0.035)       /* kpc, scale height, fixed by YMW based on CO */            \
    X(t5_kgn, 1.4)         /* unitless, ellipsolloidal correction */                    \
    X(t5_ngn, 1.84)        /* cm^{-3}, normalization, fitted by YMW */                  \
    X(t5_wgn, 0.0151)      /* kpc, width of shell , fitted by YMW */                    \
    X(t5_agn, 0.1258)      /* kpc, mid line radius of shell, fitted by YMW */           \
    X(t6_offset, 0.040)    /* kpc, zylinder offset */                                   \
    X(t6_j_lb, 0.480)      /* unitless, scale factor, fitted by YMW */                  \
    X(t6_nlb1, 1.094)      /* cm^{-3}, normalization, fitted by YMW */                  \
    X(t6_detlb1, 28.4)     /* degree, longitude scaling, fitted by YMW */               \
    X(t6_wlb1, 0.0142)     /* kpc, scale length, fitted by YMW */                       \
    X(t6_hlb1, 0.1129)     /* kpc, scale height, fitted by YMW */                       \
    X(t6_thetalb1, 195.4)  /* degree, longitude position, fitted by YMW */              \
    X(t6_nlb2, 2.33)       /* cm^{-3}, normalization, fitted by YMW */                  \
    X(t6_detlb2, 14.7)     /* degree, longitude scaling, fitted by YMW */               \
    X(t6_wlb2, 0.0156)     /* kpc, scale length, fitted by YMW */                       \
    X(t6_hlb2, 0.0436)     /* kpc, scale height, fitted by YMW */                       \
    X(t6_thetalb2, 278.2)  /* degree, longitude position, fitted by YMW */              \
    X(t7_nli, 1.907)       /* cm{-3}, normalization, fitted by YMW */                   \
    X(t7_rli, 0.080)       /* kpc, loop distance */                                     \
    X(t7_wli, 0.015)       /* kpc, loop width */                                        \
    X(t7_detthetali, 30.0) /* degree, extent of cap */                                  \
    X(t7_thetali,                                                                       \
      40.0) /* degree, angle between the direction of the center of the spherical cap and the +x direction */

IMAGINE_PARAMETERS(YMW16Parameters, YMW16_PARAMETERS)

class YMW16 : public RegularScalarModel<YMW16, YMW16Parameters> {
public:
    double max_radius = 30.; // kpc

    // warp
    double t0_r_warp = 8.4;   // kpc
    double t0_theta0 = 0.;    // deg
    double t0_gamma_w = 0.14; // unitless

    // z scaling fit
    double h0 = 32.;
    double h1 = 1.6e-3;
    double h2 = 4.e-7;

    double localbubble_boundary = 0.110; // kpc

    // Thick disc
    bool do_thick_disc = true;

    // Thin disc
    bool do_thin_disc = true;

    // spiralarms
    bool do_spiral_arms = true;

    // arms are Norma-Outer, Perseus, Carina - Sagittarius, Crux-Scutum, Local
    std::array<double, 5> t3_rmin{3.35, 3.707, 3.56, 3.67, 8.21};  // initial radius, kpc
    std::array<double, 5> t3_thmin{0.77, 2.093, 3.81, 5.76, 0.96}; // initial azimuth angle, rad
    std::array<double, 5> t3_tan_pitch{0.202, 0.173, 0.183, 0.186, 0.0483};
    std::array<double, 5> t3_cos_pitch{0.98, 0.985, 0.9836, 0.983, 0.9988};
    std::array<double, 5> t3_narm{0.135, 0.129, 0.103, 0.116, 0.0057};
    // cm^{-3}, density where arm joins thin disc, fitted by YMW
    std::array<double, 5> t3_warm{.3, .5, .3, .5, .3}; // kpc, arm widths, fitted by YMW via "preliminary global fits"

    // Galactic Center
    bool do_galactic_center = true;
    double Xgc = 0.050;  // kpc, X-position
    double Ygc = 0.;     // kpc, Y-position
    double Zgc = -0.007; // kpc, Z-position

    // gum
    bool do_gum = true;
    double t5_lc = 264.;        // degree, longitude of gum center
    const double t5_bc = -4.;   // degree, latitude of gum center
    const double t5_dc = 0.450; // kpc, distance of gum center

    // local bubble
    bool do_local_bubble = true;
    double t6_zyl1 = 0.94; // unitless, zylinder scaling
    double t6_zyl2 = 0.34; // unitless, zylinder scaling
    // first region of enhanced density

    // loop
    bool do_loop = true;
    const double x_c = -0.010156;
    const double y_c = 8.106206;
    const double z_c = 0.010467;

    template <typename T>
    auto _z_scaling(const double &rr, const T &k, const double &h0, const double &h1, const double &h2) const;
    template <typename T> auto _cosh_scaling(const double &s, const T &a, const T &b = 0.) const;

    template <typename T> T thick(const double &zz, const double &rr, T &gd, const YMW16Parameters<T> &p) const;
    template <typename T> T thin(const double &zz, const double &rr, const T &gd, const YMW16Parameters<T> &p) const;
    template <typename T>
    T spiral(const double &xx, const double &yy, const double &zz, const double &rr, const T &gd,
             const YMW16Parameters<T> &p) const;
    template <typename T>
    T galcen(const double &xx, const double &yy, const double &zz, const YMW16Parameters<T> &p) const;
    template <typename T>
    T gum(const double &xx, const double &yy, const double &zz, const YMW16Parameters<T> &p) const;
    template <typename T>
    T localbubble(const double &xx, const double &yy, const double &zz, const double &ll, const double &Rlb,
                  const YMW16Parameters<T> &p) const;
    template <typename T>
    T nps(const double &xx, const double &yy, const double &zz, const YMW16Parameters<T> &p) const;

    template <typename T> T field(const double &x, const double &y, const double &z, const YMW16Parameters<T> &p) const;
};

}

#endif
