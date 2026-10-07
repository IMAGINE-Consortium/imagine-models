#include <algorithm>
#include <cmath>

#include "ImagineModels/units.h"
#include "ImagineModels/YMW.h"

namespace imagine {

namespace {
const double cutoff = 6.;
}

template <typename T>
T YMW16::field(const double &x, const double &y, const double &z, const YMW16Parameters<T> &p) const
{
  // YMW16 using a different Cartesian frame from our default one
  const double xx = y;
  const double yy = -x;
  const double zz = z;
  if (sqrt(xx * xx + yy * yy + zz * zz) > max_radius)
    return 0.;
  // cylindrical r
  const double rr = sqrt(xx * xx + yy * yy);
  // warp, applied to the disk components only
  double zz_w = zz;
  if (rr >= t0_r_warp)
    zz_w -= t0_gamma_w * (rr - t0_r_warp) * cos(atan2(yy, xx) - t0_theta0 / 180 * M_PI);

  T ne_comp[8]{0.};
  T gd = 0.;
  const T ne_thick = thick(zz_w, rr, gd, p);
  // longitude, in deg
  const double ec_l{atan2(xx, p.r0 - yy) * 180 / M_PI};
  // since in YMW16, Fermi Bubble is not actually contributing, we ignore FB
  if (do_thick_disc)
    ne_comp[1] = ne_thick;
  if (do_thin_disc)
    ne_comp[2] = thin(zz_w, rr, gd, p);
  if (do_spiral_arms)
    ne_comp[3] = spiral(xx, yy, zz_w, rr, gd, p);
  if (do_galactic_center)
    ne_comp[4] = galcen(xx, yy, zz, p);
  if (do_gum)
    ne_comp[5] = gum(xx, yy, zz, p);
  if (do_local_bubble)
    ne_comp[6] = localbubble(xx, yy, zz, ec_l, localbubble_boundary, p);
  if (do_loop)
    ne_comp[7] = nps(xx, yy, zz, p);

  // adding up rules
  double weight_localbubble = 0.;
  double weight_gum = 0.;
  double weight_loop = 0.;
  ne_comp[0] = ne_comp[1] + std::max(ne_comp[2], ne_comp[3]);
  // distance to local bubble
  const T rlb = sqrt(pow((yy - p.r0 - p.t6_offset) * t6_zyl1 - t6_zyl2 * zz, 2) + xx * xx);
  if (rlb > localbubble_boundary)
  { // outside local bubble
    if (ne_comp[6] > ne_comp[0] and ne_comp[6] > ne_comp[5])
      weight_localbubble = 1;
  }
  else
  { // inside local bubble
    if (ne_comp[6] > ne_comp[0])
      weight_localbubble = 1;
    else
      ne_comp[0] = p.t6_j_lb * ne_comp[1] + std::max(ne_comp[2], ne_comp[3]);
  }
  if (ne_comp[7] > ne_comp[0])
    weight_loop = 1;
  if (ne_comp[5] > ne_comp[0])
    weight_gum = 1;
  // final density
  return (1 - weight_localbubble) *
             ((1 - weight_gum) * ((1 - weight_loop) * (ne_comp[0] + ne_comp[4]) +
                                  weight_loop * ne_comp[7]) +
              weight_gum * ne_comp[5]) +
         weight_localbubble * ne_comp[6];
}

// convenience function

template <typename T>
auto YMW16::_z_scaling(const double &rr, const T &k, const double &h0, const double &h1, const double &h2) const {
double rr_pc = rr * 1000;  // temporarily converting to pc, then back
return k * (h0  + h1 * rr_pc + h2 * rr_pc * rr_pc) * 0.001;
}

template <typename T>
auto YMW16::_cosh_scaling(const double &s, const T &a, const T &b) const {
  return pow(1. / cosh((s - b) / a), 2);
}

// thick disk, also sets the radial cutoff gd of the disk components
template <typename T>
T YMW16::thick(const double &zz, const double &rr, T &gd, const YMW16Parameters<T> &p) const {
  if (abs(zz) > cutoff * p.t1_h1 or (rr - p.t1_bd) > cutoff * p.t1_ad) {
    gd = 0.;
    return 0.;
  }
  if (rr < p.t1_bd)
    gd = 1.;
  else
    gd = _cosh_scaling(rr, p.t1_ad, p.t1_bd);
  return p.t1_n1 * gd * _cosh_scaling(zz, p.t1_h1);
}

// thin disk
template <typename T>
T YMW16::thin(const double &zz, const double &rr, const T &gd, const YMW16Parameters<T> &p) const
{
  // z scaling, K_2*h0 in ref
  auto k2h = _z_scaling(rr, p.t2_k2, h0, h1, h2);
  if ((rr - p.t2_b2) > cutoff * p.t2_a2 or abs(zz) > cutoff * k2h)
    return 0.;
  return p.t2_n2 * gd * _cosh_scaling(rr, p.t2_a2, p.t2_b2) * _cosh_scaling(zz, k2h);
}

// spiral arms
template <typename T>
T YMW16::spiral(const double &xx, const double &yy,
                     const double &zz, const double &rr, const T &gd, const YMW16Parameters<T> &p) const
{
  // z scaling, K_a*h0 in ref
  auto k3h = _z_scaling(rr, p.t3_ka, h0, h1, h2);
  if (abs(zz) >= 3. or gd == 0. or abs(zz) > cutoff * k3h)
    return 0.;
  const T scaling = gd * _cosh_scaling(zz, k3h) * _cosh_scaling(rr, p.t3_aa, p.t3_b2s);
  double theta = atan2(yy, xx);
  if (theta < 0)
    theta += 2 * M_PI;
  auto arm_distance = [&](int i, double d_phi) { return std::abs(rr - t3_rmin[i] * exp(d_phi * t3_tan_pitch[i])); };
  T ne3s = 0.;
  // looping through arms
  for (int i = 0; i < 5; ++i) {
    double detrr;
    if (i != 4) {
      // Norma-Outer and Perseus have one winding before theta_min, the other arms two
      const double d_phi = theta - t3_thmin[i];
      detrr = arm_distance(i, d_phi + 2 * M_PI);
      if (d_phi >= 0)
        detrr = std::min(detrr, arm_distance(i, d_phi));
      else if (i >= 2)
        detrr = std::min(detrr, arm_distance(i, d_phi + 4 * M_PI));
    }
    else if (theta >= t3_thmin[i] and theta < 2.) { // Local arm
      detrr = arm_distance(i, theta - t3_thmin[i]);
    }
    else {
      continue;
    }
    if (detrr > cutoff * t3_warm[i])
      continue;
    const double sech2 = pow(1. / cosh(detrr * t3_cos_pitch[i] / t3_warm[i]), 2);
    const double theta_deg = theta * 180 / M_PI;
    if (i != 2) {
      ne3s += t3_narm[i] * scaling * sech2;
    }
    else if (rr > 6 and theta_deg > p.t3_thetacn)
    { // correction for Carina-Sagittarius
      const T ga = (1. - p.t3_nsg * exp(-pow((theta_deg - p.t3_thetasg) / p.t3_wsg, 2))) * (1. + p.t3_ncn) * sech2;
      ne3s += t3_narm[i] * scaling * ga;
    }
    else
    {
      const T ga = (1. - p.t3_nsg * exp(-pow((theta_deg - p.t3_thetasg) / p.t3_wsg, 2))) *
                   (1. + p.t3_ncn * exp(-pow((theta_deg - p.t3_thetacn) / p.t3_wcn, 2))) * sech2;
      ne3s += t3_narm[i] * scaling * ga;
    }
  } // end of looping through arms
  return ne3s;
}

// galactic center
template <typename T>
T YMW16::galcen(const double &xx, const double &yy, const double &zz, const YMW16Parameters<T> &p) const
{
  const double Rgc = sqrt((xx - Xgc) * (xx - Xgc) + (yy - Ygc) * (yy - Ygc));
  if (Rgc > cutoff * p.t4_agc or abs(zz) > cutoff * p.t4_hgc)
    return 0.;
  return p.t4_ngc * exp(-Rgc * Rgc / (p.t4_agc * p.t4_agc)) * pow(1. / cosh((zz - Zgc) / p.t4_hgc), 2);
}

// gum nebula
template <typename T>
T YMW16::gum(const double &xx, const double &yy, const double &zz, const YMW16Parameters<T> &p) const
{
  // center of Gum Nebula
  const double rgalc = t5_dc * cos(t5_bc * M_PI / 180);
  const double xc = rgalc * sin(t5_lc * M_PI / 180);
  const T yc = p.r0 - rgalc * cos(t5_lc * M_PI / 180);
  const double zc = t5_dc * sin(t5_bc * M_PI / 180);
  const T theta = abs(atan((zz - zc) / sqrt((xx - xc) * (xx - xc) + (yy - yc) * (yy - yc))));
  const T RR = sqrt((xx - xc) * (xx - xc) + (yy - yc) * (yy - yc) + (zz - zc) * (zz - zc));
  T Dmin;
  if (theta == 0.) {
    // limit theta -> 0 (the original divides 0 by 0 here)
    Dmin = abs(RR - p.t5_agn);
  }
  else {
    const T tantheta = tan(theta);
    const T zp = p.t5_agn * p.t5_kgn / sqrt(1. + p.t5_kgn * p.t5_kgn / (tantheta * tantheta));
    const T xyp = zp / tantheta;
    T alpha;
    if (p.t5_agn - abs(xyp) < 1e-15)
      alpha = M_PI / 2;
    else
      alpha = atan(p.t5_kgn * xyp / sqrt(p.t5_agn * p.t5_agn - xyp * xyp));
    Dmin = abs((RR - sqrt(zp * zp + xyp * xyp)) * sin(theta + alpha));
  }
  if (Dmin > cutoff * p.t5_wgn)
    return 0.;
  return p.t5_ngn * exp(-Dmin * Dmin / (p.t5_wgn * p.t5_wgn));
}

// local bubble
template <typename T>
T YMW16::localbubble(const double &xx, const double &yy, const double &zz, const double &ll,
                          const double &Rlb, const YMW16Parameters<T> &p) const
{
  const T y_lb = p.r0 + p.t6_offset;
  // r_LB in ref
  const T rLB = sqrt(pow((yy - y_lb) * t6_zyl1 - t6_zyl2 * zz, 2) + xx * xx);
  // the first region is limited to a cone of 45 deg half opening angle
  const T y_axis = y_lb + t6_zyl2 / t6_zyl1 * zz;
  const T cos_a = (xx * xx + (p.r0 - yy) * (y_axis - yy)) /
                  (sqrt(xx * xx + (p.r0 - yy) * (p.r0 - yy)) * sqrt(xx * xx + (y_axis - yy) * (y_axis - yy)));
  // l-l_LB in ref
  const T dl1 = std::min(abs(ll + 360. - p.t6_thetalb1), abs(p.t6_thetalb1 - ll));
  const T dl2 = std::min(abs(ll + 360. - p.t6_thetalb2), abs(p.t6_thetalb2 - ll));
  T nel1 = 0.;
  if ((rLB - Rlb) <= cutoff * p.t6_wlb1 and abs(zz) <= cutoff * p.t6_hlb1 and dl1 <= cutoff * p.t6_detlb1)
    nel1 = p.t6_nlb1 * pow(1. / cosh(dl1 / p.t6_detlb1), 2) * pow(1. / cosh((rLB - Rlb) / p.t6_wlb1), 2) *
           pow(1. / cosh(zz / p.t6_hlb1), 2);
  // as in the original, the cuts of the second region only apply if the first one vanishes
  if (nel1 == 0. and ((rLB - Rlb) > cutoff * p.t6_wlb2 or abs(zz) > cutoff * p.t6_hlb2))
    return 0.;
  T nel2 = 0.;
  if (dl2 <= cutoff * p.t6_detlb2)
    nel2 = p.t6_nlb2 * pow(1. / cosh(dl2 / p.t6_detlb2), 2) * pow(1. / cosh((rLB - Rlb) / p.t6_wlb2), 2) *
           pow(1. / cosh(zz / p.t6_hlb2), 2);
  if (cos_a < cos(M_PI / 4))
    nel1 = 0.;
  return nel1 + nel2;
}

// north polar spur
template <typename T>
T YMW16::nps(const double &xx, const double &yy, const double &zz, const YMW16Parameters<T> &p) const
{
  const T theta_LI = p.t7_thetali / 180. * M_PI;
  // r_LI in ref
  const double rLI = sqrt((xx - x_c) * (xx - x_c) + (yy - y_c) * (yy - y_c) + (zz - z_c) * (zz - z_c));
  const T theta = acos(((xx - x_c) * cos(theta_LI) + (zz - z_c) * sin(theta_LI)) / rLI) * 180. / M_PI;
  if (abs(rLI - p.t7_rli) > cutoff * p.t7_wli or abs(theta) > cutoff * p.t7_detthetali)
    return 0.;
  return p.t7_nli * exp(-pow((rLI - p.t7_rli) / p.t7_wli, 2)) * exp(-pow(theta / p.t7_detthetali, 2));
}


IMAGINE_INSTANTIATE_SCALAR_MODEL(YMW16)

}
