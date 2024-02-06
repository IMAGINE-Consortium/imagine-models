#include <cassert>

#include "hamunits.h"
#include "YMW.h"

number YMW16::_at_position(const double &x, const double &y, const double &z, const YMW16 &p) const
{
  // YMW16 using a different Cartesian frame from our default one
  std::array<double, 3> gc_pos{y, -x, z};
  // cylindrical r
  double r_cyl{sqrt(gc_pos[0] * gc_pos[0] + gc_pos[1] * gc_pos[1])};
  // warp
  if (r_cyl >= p.t0_r_warp) {
    double theta_warp{atan2(gc_pos[1], gc_pos[0])};
    gc_pos[2] -= p.t0_gamma_w * (r_cyl - p.t0_r_warp) * cos(theta_warp - p.t0_theta0 / 180 * M_PI);
  }
  double vec_length = sqrt(pow(gc_pos[0], 2) + pow(gc_pos[1], 2) + pow(gc_pos[2], 2));
  if (vec_length > 25)
  {
    return 0.;
  }
  else
  {
    number ne{0.};
    number ne_comp[8]{0.};
    double weight_localbubble{0.};
    double weight_gum{0.};
    double weight_loop{0.};
    // longitude, in deg
    const double ec_l{atan2(gc_pos[0], p.r0 - gc_pos[1]) * 180 / M_PI};
    // call structure functions
    // since in YMW16, Fermi Bubble is not actually contributing, we ignore FB
    if (do_thick_disc) {
      ne_comp[1] = thick(gc_pos[2], r_cyl, p);
    }
    if (do_thin_disc) {
      ne_comp[2] = thin(gc_pos[2], r_cyl, p);
    }
    if (do_spiral_arms) {
      ne_comp[3] = spiral(gc_pos[0], gc_pos[1], gc_pos[2], r_cyl, p);
    }
    if (do_galactic_center) {
      ne_comp[4] = galcen(gc_pos[0], gc_pos[1], gc_pos[2], p);
    }
    if (do_gum) {
      ne_comp[5] = gum(gc_pos[0], gc_pos[1], gc_pos[2], p);
    }
    if (do_local_bubble) {
      ne_comp[6] = localbubble(gc_pos[0], gc_pos[1], gc_pos[2], ec_l,
                             localbubble_boundary, p);
    }
    if (do_loop) {
      ne_comp[7] = nps(gc_pos[0], gc_pos[1], gc_pos[2], p);
    } 
   
    // adding up rules
    ne_comp[0] = ne_comp[1] + std::max(ne_comp[2], ne_comp[3]);
    // distance to local bubble
    const double rlb{sqrt(pow(((gc_pos[1] - p.r0 - p.t6_offset) * p.t6_zyl1 - p.t6_zyl2 * gc_pos[2]), 2) + gc_pos[0] * gc_pos[0])};
    if (rlb < localbubble_boundary)
    { // inside local bubble
      ne_comp[0] = rlb * ne_comp[1] +
                   std::max(ne_comp[2], ne_comp[3]);
      if (ne_comp[6] > ne_comp[0])
      {
        weight_localbubble = 1;
      }
    }
    else
    { // outside local bubble
      if (ne_comp[6] > ne_comp[0] and ne_comp[6] > ne_comp[5])
      {
        weight_localbubble = 1;
      }
    }
    if (ne_comp[7] > ne_comp[0])
    {
      weight_loop = 1;
    }
    if (ne_comp[5] > ne_comp[0])
    {
      weight_gum = 1;
    }
    // final density
    ne =
        (1 - weight_localbubble) *
            ((1 - weight_gum) * ((1 - weight_loop) * (ne_comp[0] + ne_comp[4]) +
                                 weight_loop * ne_comp[7]) +
             weight_gum * ne_comp[5]) +
        (weight_localbubble) * (ne_comp[6]);
    // assert(std::isfinite(ne));
    #if !autodiff_FOUND
    if (std::isnan(ne)) {
      std::cout << "Found nan at: (x,y,z): ()" << x << ", " << y << ", " << z << ")" << std::endl;
      if (std::isnan(ne_comp[1])) {
        std::cout << "Nan in thick disc" << std::endl;
      }
      if (std::isnan(ne_comp[2])) {
        std::cout << "Nan in thin disc" << std::endl;
      } 
      if (std::isnan(ne_comp[3])) {
        std::cout << "Nan in spiral" << std::endl;
      }
      if (std::isnan(ne_comp[4])) {
        std::cout << "Nan in galcen" << std::endl;
      }
      if (std::isnan(ne_comp[5])) {
        std::cout << "Nan in gum" << std::endl;
      }
      if (std::isnan(ne_comp[6])) {
        std::cout << "Nan in local bubble" << std::endl;
      }
      if (std::isnan(ne_comp[7])) {
        std::cout << "Nan in loop" << std::endl;
      }
    }
    #endif
    //std::cout << "ne: " << ne << std::endl;
    return ne;
  }
}

// convenience function

auto YMW16::_z_scaling(const double &rr, const number &k, const double &h0, const double &h1, const double &h2) const {
double rr_pc = rr * 1000;  // temporarily converting to pc, then back 
return k * (h0  + h1 * rr_pc + h2 * rr_pc * rr_pc) * 0.001;
}

auto YMW16::_cosh_scaling(const double &s, const number &a, const number &b) const {
  return pow(1. / cosh((s - b) / a), 2);
}

// thick disk
number YMW16::thick(const double &zz, const double &rr, const YMW16 &p) const {
  if (zz > 10. * p.t1_h1)
    return 0.; // timesaving
  number gd = 1.;  
  if (rr > p.t1_bd) {
    gd = _cosh_scaling(rr, p.t1_bd, p.t1_ad);
  }
  return p.t1_n1 * gd *_cosh_scaling(zz, p.t1_h1);
}

// thin disk
number YMW16::thin(const double &zz, const double &rr, const YMW16 &p) const
{
  // z scaling, K_2*h0 in ref
  auto k2h = _z_scaling(rr, p.t2_k2, p.h0, p.h1, p.h2); 
  if (zz > 10. * k2h)
    return 0.; // timesaving
  number gd = 1.;  
  if (rr > p.t1_bd) {
    gd = _cosh_scaling(rr, p.t1_bd, p.t1_ad);
  }
  auto gd2 =  _cosh_scaling(rr, p.t2_a2,  p.t2_b2); // pow(1. / cosh((rr -  p.t2_b2) / p.t2_a2), 2);

  return p.t2_n2 * gd * gd2 * _cosh_scaling(zz, k2h);
}

// spiral arms
number YMW16::spiral(const double &xx, const double &yy,
                     const double &zz, const double &rr, const YMW16 &p) const
{
  // structure scaling
  number scaling = 1.;  
  if (rr > p.t1_bd) {
    scaling = _cosh_scaling(rr, p.t1_bd, p.t1_ad);
  }
  // z scaling, K_a*h0 in ref
  auto k3h = _z_scaling(rr, p.t3_ka, p.h0, p.h1, p.h2); 

  if (abs(zz) > 10. * k3h)
    return 0.; // timesaving
  scaling *= _cosh_scaling(zz, k3h);
  if ((rr - p.t3_b2s) > 10. * p.t3_aa)
    return 0.; // timesaving
  // 2nd raidus scaling
  scaling *= _cosh_scaling(rr, p.t3_aa, p.t3_b2s);
  number smin;
  double theta{atan2(yy, xx)};
  if (theta < 0)
    theta += 2 * M_PI;
  number ne3s{0.};
  // looping through arms
  for (int i = 0; i < 5; ++i) { 
    number phimin = p.t3_phimin[i] / 180 * M_PI;
    number tpitch = tan(p.t3_tpitch[i] / 180 * M_PI);
    // get distance to arm center
    if (i != 4) { 
      number d_phi = theta - phimin;
      if (d_phi < 0) {
        d_phi += 2. * M_PI;
      }
      number d = abs(p.t3_rmin[i] * exp(d_phi * tpitch) - rr);
      number d_p = abs(p.t3_rmin[i] * exp((d_phi + 2. * M_PI) * tpitch) - rr);
      // smin = std::min(d, d_p) * tpitch;
      smin = std::min(d, d_p); // * tpitch;
    }
    else if (i == 4 and theta >= phimin and theta < (2 / 180 * M_PI)) { // Local arm
      smin = abs(p.t3_rmin[i] * exp((theta + 2 * M_PI - phimin) * tpitch) - rr);
    }
    else {
      continue;
    }
    if (smin > 10. * p.t3_warm[i])
      continue; // timesaving
    // accumulate density
    if (i != 2) {
      ne3s += p.t3_narm[i] * scaling * pow(1. / cosh(smin / p.t3_warm[i]), 2);
    }
    else if (rr > 6 and
             theta * 180 / M_PI > p.t3_thetacn)
    { // correction for Carina-Sagittarius
      const number ga =
          (1. - (p.t3_nsg) * (exp(-pow((theta  * 180 / M_PI - p.t3_thetasg) / p.t3_wsg, 2)))) *
          (1. + p.t3_ncn) * pow(1. / cosh(smin / p.t3_warm[i]), 2);
      ne3s += p.t3_narm[i] * scaling * ga;
    }
    else
    {
      const number ga =
          (1. - (p.t3_nsg) * (exp(-pow((theta  * 180 / M_PI - p.t3_thetasg) / p.t3_wsg, 2)))) *
          (1. + p.t3_ncn * exp(-pow((theta  * 180 / M_PI - p.t3_thetacn) / p.t3_wcn, 2))) *
          pow(1. / cosh(smin / p.t3_warm[i]), 2);
      ne3s += p.t3_narm[i] * scaling * ga;
    }
  } // end of looping through arms
  return ne3s;
}

// galactic center
number YMW16::galcen(const double &xx, const double &yy, const double &zz, const YMW16 &p) const
{
  // pos of center
  const double R2gc{(xx - p.Xgc) * (xx - p.Xgc) + (yy - p.Ygc) * (yy - p.Ygc)};
  if (R2gc > 10. * p.t4_agc * p.t4_agc)
    return 0.; // timesaving
  const double Ar{exp(-R2gc / (p.t4_agc * p.t4_agc))};
  if (abs(zz - p.Zgc) > 10. * p.t4_hgc)
    return 0.; // timesaving
  const double Az{pow(1. / cosh((zz - p.Zgc) / p.t4_hgc), 2)};
  return p.t4_ngc * Ar * Az;
}

// gum nebula
number YMW16::gum(const double &xx, const double &yy, const double &zz, const YMW16 &p) const
{
  if (yy < 0 or xx > 0)
    return 0.; // timesaving
  // center of Gum Nebula

  const double xc{p.t5_dc * cos(p.t5_bc * M_PI / 180) * sin(p.t5_lc * M_PI / 180)};
  const double yc{p.r0 - p.t5_dc * cos(p.t5_bc * M_PI / 180) * cos(p.t5_lc * M_PI / 180)};
  const double zc{p.t5_dc * sin(p.t5_bc * M_PI / 180)};
  // theta is limited in I quadrant
  const double thetagum{
      atan2(abs(zz - zc),
            sqrt((xx - xc) * (xx - xc) + (yy - yc) * (yy - yc)))};
  const double tantheta = tan(thetagum);
  // zp is positive
  number zp = 0;
  number xyp = 0;
  if (tantheta != 0.) {
    zp +=  (p.t5_agn * p.t5_kgn)/sqrt(1. + p.t5_kgn * p.t5_kgn / (tantheta * tantheta));
    // xyp is positive
    xyp += zp / tantheta;
  } 
  // alpha is positive
  const number xy_dist = {
      sqrt(p.t5_agn * p.t5_agn - xyp * xyp) *
      double(p.t5_agn > xyp)};
  const double alpha{atan2(p.t5_kgn * xyp, xy_dist) +
                     thetagum}; // add theta, timesaving
  const double R2{(xx - xc) * (xx - xc) + (yy - yc) * (yy - yc) + (zz - zc) * (zz - zc)};
  const double r2{zp * zp + xyp * xyp};
  const double D2min{(R2 + r2 - 2. * sqrt(R2 * r2)) * sin(alpha) * sin(alpha)};
  if (D2min > 10. * p.t5_wgn * p.t5_wgn)
    return 0.;
  return p.t5_ngn * exp(-D2min / (p.t5_wgn * p.t5_wgn));
}

// local bubble
number YMW16::localbubble(const double &xx, const double &yy, const double &zz, const double &ll,
                          const double &Rlb, const YMW16 &p) const
{
  if (yy < 0)
    return 0.; // timesaving
  number nel{0.};
  // r_LB in ref
  auto rLB{
      sqrt(pow(((yy - p.r0 - p.t6_offset) * p.t6_zyl1 - p.t6_zyl2 * zz), 2) + pow(xx, 2))};
  // l-l_LB1 in ref

  auto dl1 = std::min(abs(ll + 360. - p.t6_thetalb1), abs(p.t6_thetalb1 - ll));
  if (dl1 < 10. * p.t6_detlb1 or
      (rLB - Rlb) < 10. * p.t6_wlb1 or
      zz < 10. * p.t6_hlb1) // timesaving
    nel += p.t6_nlb1 *
           pow(1. / cosh(dl1 / p.t6_detlb1), 2) *
           pow(1. / cosh((rLB - Rlb) / p.t6_wlb1), 2) *
           pow(1. / cosh(zz / p.t6_hlb1), 2);
  // l-l_LB2 in ref
  auto dl2{
      std::min(abs(ll + 360. - p.t6_thetalb2),
               abs(p.t6_thetalb2 - (ll)))};
  if (dl2 < 10. * p.t6_detlb2 or
      (rLB - Rlb) < 10. * p.t6_wlb2 or
      zz < 10. * p.t6_hlb2) // timesaving
    nel += p.t6_nlb2 *
           pow(1. / cosh(dl2 / p.t6_detlb2), 2) *
           pow(1. / cosh((rLB - Rlb) / p.t6_wlb2), 2) *
           pow(1. / cosh(zz / p.t6_hlb2), 2);
  return nel;
}

// north polar spur
number YMW16::nps(const double &xx, const double &yy, const double &zz, const YMW16 &p) const
{
  if (yy < 0)
    return 0.; // timesaving
  const number theta_LI = p.t7_thetali / 180. * M_PI;
  // r_LI in ref
  const double rLI{sqrt((xx - p.x_c) * (xx - p.x_c) +
                        (yy - p.y_c) * (yy - p.y_c) +
                        (zz - p.z_c) * (zz - p.z_c))};
  const number theta{acos(((xx - p.x_c) * (cos(theta_LI)) +
                           (zz - p.z_c) * (sin(theta_LI))) /
                          rLI)
                     * 180. / M_PI};
  if (theta > 10. * p.t7_detthetali or
      (rLI - p.t7_rli) > 10. * p.t7_wli)
    return 0.; // timesaving
  return (p.t7_nli) *
         exp(-pow((rLI - p.t7_rli) / p.t7_wli, 2)) *
         exp(-pow(theta / p.t7_detthetali, 2));
}

#if autodiff_FOUND

Eigen::VectorXd YMW16::_jac(const double &x, const double &y, const double &z, YMW16 &p) const
{
  ad::real out;
  Eigen::VectorXd _deriv = ad::gradient([&](double _x, double _y, double _z, YMW16 &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(
                                            p.r0,
                                            p.t1_ad, p.t1_bd, p.t1_n1, p.t1_h1,
                                            p.t2_a2, p.t2_b2, p.t2_n2, p.t2_k2,
                                            p.t3_b2s, p.t3_ka, p.t3_aa, p.t3_ncn, p.t3_wcn, p.t3_thetacn, p.t3_nsg, p.t3_wsg, p.t3_thetasg,
                                            p.t4_ngc, p.t4_agc, p.t4_hgc,
                                            p.t5_kgn, p.t5_ngn, p.t5_wgn, p.t5_agn,
                                            p.t6_j_lb, p.t6_nlb1, p.t6_detlb1, p.t6_wlb1, p.t6_hlb1, p.t6_thetalb1, p.t6_nlb2, p.t6_detlb2, p.t6_wlb2, p.t6_hlb2, p.t6_thetalb2,
                                            p.t7_nli, p.t7_rli, p.t7_wli, p.t7_detthetali, p.t7_thetali),
                                        ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif