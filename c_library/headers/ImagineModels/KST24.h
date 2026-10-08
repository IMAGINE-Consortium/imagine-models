// Reference: Korochkin, Semikoz & Tinyakov 2025, A&A 693, A284, arXiv:2407.02148 (Table 2)
// Based on: CRPropa (KST24Field), GPL-3.0; authors' code (Zenodo 14743599), CC-BY-4.0
// Deviations:
// - from the authors' code, not in the paper: Sagittarius-Carina arm widening by 3 deg along the arm, radial (3-17 kpc) and vertical arm cut-offs, arm widths capped at 1.2 kpc, spiral scale a = 3 kpc
// - Sagittarius-Carina rdisk = 0.79 kpc as in the code (Table 2: 0.8)
// - outer Perseus field -3.5 muG as in Table 2 and CRPropa; the authors' Zenodo code uses -2.5 muG

#pragma once

#include <array>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define KST24_PARAMETERS(X)      \
    X(pitch, 20.) /* deg */      \
    X(B_local, -3.5)             \
    X(b_local, -2.2) /* deg */   \
    X(n_local, 1.45)             \
    X(rz_local, 0.73)            \
    X(rdisk_local, 1.)           \
    X(x0_local, -0.15)           \
    X(B_sagcar, 1.3)             \
    X(b_sagcar, -80.) /* deg */  \
    X(n_sagcar, 2.3)             \
    X(rz_sagcar, 1.)             \
    X(rdisk_sagcar, 0.79)        \
    X(x0_sagcar, 1.37)           \
    X(B_scutum, 4.9)             \
    X(b_scutum, -134.) /* deg */ \
    X(n_scutum, 2.)              \
    X(rz_scutum, 0.8)            \
    X(rdisk_scutum, 1.)          \
    X(x0_scutum, 1.)             \
    X(B_perseus, -3.5)           \
    X(Bprime_perseus, -2.)       \
    X(b_perseus, 46.) /* deg */  \
    X(n_perseus, 2.)             \
    X(rz_perseus, 1.2)           \
    X(rzprime_perseus, 0.4)      \
    X(rdisk_perseus, 1.1)        \
    X(x0_perseus, -1.)           \
    X(B_ntor, 3.2)               \
    X(zmin_ntor, 1.185)          \
    X(zmax_ntor, 2.1)            \
    X(rmax_ntor, 9.1)            \
    X(B_stor, -3.2)              \
    X(zmin_stor, -2.5)           \
    X(zmax_stor, -1.22)          \
    X(rmax_stor, 14.)            \
    X(B_X, 1.8)                  \
    X(rmax_X, 6.2)               \
    X(theta_X, 28.) /* deg */    \
    X(x_LB, -8.2)                \
    X(y_LB, 0.095)               \
    X(z_LB, -0.05)               \
    X(r_LB, 0.2)                 \
    X(dr_LB, 0.03)               \
    X(l_LB, 230.) /* deg */      \
    X(b_LB, -2.)  /* deg */

IMAGINE_PARAMETERS(KST24Parameters, KST24_PARAMETERS)

class KST24MagneticField : public RegularVectorModel<KST24MagneticField, KST24Parameters> {
public:
    // local, sagcar, scutum, perseus x2
    std::array<double, 5> arm_rmin{3.3, 3., 3., 3., 3.};
    std::array<double, 5> arm_rmax{14., 15., 16., 17., 17.};
    std::array<double, 5> arm_zmax{1., 1.2, 1., 1., 2.};
    std::array<double, 5> arm_widening{0., 3., 0., 0., 0.}; // deg
    double spiral_a = 3.;
    double width_reference_r = 5.;
    double width_max = 1.2;
    double torus_rmin = 1.;
    double X_rmin = 1.;
    double X_zmax = 10.;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const KST24Parameters<T> &p) const;

private:
    template <typename T>
    Vec3<T> arm(const double &x, const double &y, const double &z, const T &B, const T &pitch, const T &phase,
                const T &x0, const T &rz, const T &rdisk, const T &n, int i) const;
    template <typename T>
    Vec3<T> torus(const double &x, const double &y, const double &z, const T &B, const T &zmin, const T &zmax,
                  const T &rmax) const;
    template <typename T>
    Vec3<T> xfield(const double &x, const double &y, const double &z, const KST24Parameters<T> &p) const;
    template <typename T>
    Vec3<T> bubble(const double &x, const double &y, const double &z, const KST24Parameters<T> &p) const;
};

}
