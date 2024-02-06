#ifndef YMW16_H
#define YMW16_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>
#include <cassert>

#include "Field.h"
#include "RegularField.h"

class YMW16 : public RegularScalarField
{
protected:
  number _at_position(const double &x, const double &y, const double &z, const YMW16 &p) const;

#if autodiff_FOUND
  Eigen::VectorXd _jac(const double &x, const double &y, const double &z, YMW16 &p) const;
#endif
public:
  using RegularScalarField ::RegularScalarField;

  number r0 = 8.3;     // kpc, Galactic earth position

  // warp
  double t0_r_warp = 8.4; // kpc
  double t0_theta0 = 0.; // deg
  double t0_gamma_w = 0.14; // unitless

  // z scaling fit
  double h0 = 32.;
  double h1 = 1.6e-3;
  double h2 = 4.e-7;

  double localbubble_boundary = 0.110; // kpc 

  // Thick disc
  bool do_thick_disc = true;
  number t1_ad = 2.5; // kpc scale length of cutoff
  number t1_bd = 15.; // kpc radius of begin of cutoff
  number t1_n1 = 0.01132;  // cm^{-3}, normalization fitted by YMW
  number t1_h1 = 1.673; // kpc, scale-height, fitted by YMW

  // Thin disc
  bool do_thin_disc = true;
  number t2_a2 = 1.2; // kpc, scale length
  number t2_b2 = 4.; // kpc, molecular ring central radius
  number t2_n2 = 0.404; // cm^{-3}, fitted by YMW
  number t2_k2 = 1.54;  // unitless, rescaling of scale height, fitted by YMW

  // spiralarms
  bool do_spiral_arms = true;
  number t3_b2s = 4.; // kpc, is actually t2_b2 in paper
  number t3_ka = 5.015; // spiral arm scale factor
  number t3_aa = 11.680; // kpc, scale_length

  // Carina corrections
  number t3_ncn = 2.4; // unitless, Carina relative over density, fitted by YMW
  number t3_wcn = 8.2; // // degree, theta scaling, fitted by YMW
  number t3_thetacn = 109.; // degree, theta correction, fitted by YMW

  // Sagittarius correction
  number t3_nsg = 0.626;// unitless, Carina relative under density, fitted by YMW
  number t3_wsg = 20;// degree, theta scaling, fitted by YMW
  number t3_thetasg = 78.8; // degree, theta correction, fitted by YMW
  // arms are Norma-Outer, Perseus, Carina - Sagittarius, Crux-Scutum, Local 
  std::array<double, 5> t3_rmin{3.35, 3.707, 3.56, 3.670, 8.21};  // initial radius, kpc, Hou & Han fit
  std::array<double, 5> t3_phimin{44.4, 120.0, 218.6, 330.3, 55.1};
  // intial azimuth angle, degree, Hou & Han fit
  std::array<double, 5> t3_tpitch{11.43, 9.84, 10.38, 10.54, 2.77}; // pitch angle, degree, Hou & Han fit
  std::array<double, 5> t3_narm{0.135, 0.129, 0.103, 0.116, 0.0057};
  // cm^{-3}, density where arm joins thin disc, fitted by YMW
  std::array<double, 5> t3_warm{.3, .5, .3, .5, .3}; // kpc, arm widths, fitted by YMW via "preliminary global fits"


  // Galactic Center
  bool do_galactic_center = true;
  double Xgc = 0.050; // kpc, X-position  
  double Ygc = 0.; // kpc, Y-position  
  double Zgc = -0.007; // kpc, Z-position  
  number t4_ngc = 6.2;  // cm^{-3}, normalization, fitted by YMW
  number t4_agc = 0.160; // kpc, scale length, fixed by YMW based on CO
  number t4_hgc = 0.035; // kpc, scale height, fixed by YMW based on CO

  // gum
  bool do_gum = true;
  double t5_lc = 264.; // degree, longitude of gum center
  const double t5_bc = -4.; // degree, latitude of gum center
  const double t5_dc = 0.450; // kpc, distance of gum center
  number t5_kgn = 1.4; // unitless, ellipsolloidal correction
  number t5_ngn = 1.84; // cm^{-3}, normalization, fitted by YMW
  number t5_wgn = 0.0151; // kpc, width of shell , fitted by YMW
  number t5_agn = 0.1258; // kpc, mid line radius of shell, fitted by YMW 

  // local bubble
  bool do_local_bubble = true;
  double t6_zyl1 = 0.94; // unitless, zylinder scaling
  double t6_zyl2 = 0.34; // unitless, zylinder scaling
  number t6_offset = 0.004; // kpc, zylinder offset
  number t6_j_lb = 0.480; // unitless, scale factor, fitted by YMW
  // first region of enhanced density
  number t6_nlb1 = 1.094; // cm^{-3}, normalization, fitted by YMW
  number t6_detlb1 = 28.4; // degree, longitude scaling, fitted by YMW
  number t6_wlb1 = 0.0142; // kpc, scale length, fitted by YMW
  number t6_hlb1 = 0.1129; // kpc, scale height, fitted by YMW
  number t6_thetalb1 = 195.4; // degree, longitude position, fitted by YMW
  // second region of enhanced density
  number t6_nlb2 = 2.33; // cm^{-3}, normalization, fitted by YMW
  number t6_detlb2 = 14.7; // degree, longitude scaling, fitted by YMW
  number t6_wlb2 = 0.0156; // kpc, scale length, fitted by YMW
  number t6_hlb2 = 0.0436; // kpc, scale height, fitted by YMW
  number t6_thetalb2 = 278.2; // degree, longitude position, fitted by YMW

  // loop
  bool do_loop = true;
  const double x_c = -0.010156 ;
  const double y_c = 8.106207;
  const double z_c = 0.010467;
  number t7_nli = 1.907;  // cm{-3}, normalization, fitted by YMW 
  number t7_rli = 0.080; // kpc, loop distance
  number t7_wli = 0.015; // kpc, loop width
  number t7_detthetali = 30.0; //degree, extent of cap
  number t7_thetali = 40.0; // degree, angle between the direction of the center of the spherical cap and the +x direction
#if autodiff_FOUND
  const std::set<std::string> all_diff{
      "r0", "t1_ad", "t1_bd", "t1_n1", "t1_h1", "t2_a2", "t2_b2", "t2_n2", "t2_k2",
      "t3_b2s", "t3_ka", "t3_aa", "t3_ncn", "t3_wcn", "t3_thetacn", "t3_nsg", "t3_wsg", "t3_thetasg",
      "t4_ngc", "t4_agc", "t4_hgc", "t5_kgn", "t5_ngn", "t5_wgn", "t5_agn", "t6_j_lb", "t6_nlb1",
      "t6_detlb1", "t6_wlb1", "t6_hlb1", "t6_thetalb1", "t6_nlb2", "t6_detlb2", "t6_wlb2", "t6_hlb2", "t6_thetalb2",
      "t7_nli", "t7_rli", "t7_wli", "t7_detthetali", "t7_thetali"};
  std::set<std::string> active_diff{
      "r0", "t1_ad", "t1_bd", "t1_n1", "t1_h1", "t2_a2", "t2_b2", "t2_n2", "t2_k2",
      "t3_b2s", "t3_ka", "t3_aa", "t3_ncn", "t3_wcn", "t3_thetacn", "t3_nsg", "t3_wsg", "t3_thetasg",
      "t4_ngc", "t4_agc", "t4_hgc", "t5_kgn", "t5_ngn", "t5_wgn", "t5_agn", "t6_j_lb", "t6_nlb1",
      "t6_detlb1", "t6_wlb1", "t6_hlb1", "t6_thetalb1", "t6_nlb2", "t6_detlb2", "t6_wlb2", "t6_hlb2", "t6_thetalb2",
      "t7_nli", "t7_rli", "t7_wli", "t7_detthetali", "t7_thetali"};

  Eigen::VectorXd derivative(const double &x, const double &y, const double &z)
  {
    return _jac(x, y, z, *this);
  }
#endif

  auto _z_scaling(const double &rr, const number &k, const double &h0, const double &h1, const double &h2) const;
  auto _cosh_scaling(const double &s, const number &a, const number &b = 0. ) const;

  number thick(const double &zz, const double &rr, const YMW16 &p) const;
  number thin(const double &zz, const double &rr, const YMW16 &p) const;
  number spiral(const double &xx, const double &yy, const double &zz,
                const double &rr, const YMW16 &p) const;
  number galcen(const double &xx, const double &yy, const double &zz, const YMW16 &p) const;
  number gum(const double &xx, const double &yy, const double &zz, const YMW16 &p) const;
  number localbubble(const double &xx, const double &yy,
                     const double &zz, const double &ll,
                     const double &Rlb, const YMW16 &p) const;
  number nps(const double &xx, const double &yy, const double &zz, const YMW16 &p) const;

  number at_position(const double &x, const double &y, const double &z) const
  {
    return _at_position(x, y, z, *this);
  }
};

#endif
