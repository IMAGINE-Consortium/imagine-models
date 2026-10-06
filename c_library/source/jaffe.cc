#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include "ImagineModels/units.h"
#include "ImagineModels/Jaffe.h"

namespace imagine {

template <typename T>
Vec3<T> JaffeMagneticField::field(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const
{
  if (x == 0. && y == 0. && z == 0.)
  {
    return Vec3<T>{{0., 0., 0.}};
  }
  T inner_b{0};
  if (ring)
  {
    inner_b = p.ring_amp;
  }
  else if (bar)
  {
    inner_b = p.bar_amp;
  }

  Vec3<T> bhat = orientation(x, y, z, p);
  Vec3<T> btot{{0., 0., 0.}};

  auto scaling = radial_scaling(x, y, p) *
                 (p.disk_amp * disk_scaling(z, p) +
                  p.halo_amp * halo_scaling(z, p));

  for (int i = 0; i < bhat.size(); ++i)
  {
    btot[i] = bhat[i] * scaling;
  }

  


  // compress factor for each arm or for ring/bar
  std::vector<T> arm = arm_compress(x, y, z, p);
  // only inner region
  if (arm.size() == 1)
  {
    for (int i = 0; i < bhat.size(); ++i)
    {
      btot[i] += bhat[i] * arm[0] * inner_b;
    }
  }
 // return btot;}

  // spiral arm region
  else
  {
    std::array<T, 4> arm_amp = {p.arm_amp1, p.arm_amp2, p.arm_amp3, p.arm_amp4};
    for (decltype(arm.size()) i = 0; i < arm.size(); ++i)
    {
      for (int j = 0; j < bhat.size(); ++j)
      {
        btot[j] += bhat[j] * arm[i] * arm_amp[i];
      }
    }
  }
  return btot;
}



template <typename T>
Vec3<T> JaffeMagneticField::orientation(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const
{
  if (x == 0. && y == 0.)
  {
    return Vec3<T>{{0., 0., 0.}};
  }

  const double r{
      sqrt(x * x + y * y)}; // cylindrical frame
  const auto r_lim = p.ring_r;
  const auto bar_lim{p.bar_a + 0.5 * p.comp_d};
  auto arm_pitch = p.arm_pitch * M_PI /180;
  const auto cos_p = cos(arm_pitch);
  const auto sin_p = sin(arm_pitch); // pitch angle

  Vec3<T> tmp{{0., 0., 0.}};
  T quadruple{1.};
  if (r < 0.5) // forbidden region
    return tmp;
  if (z > p.disk_z0)
    quadruple = (1 - 2 * this->quadruple);
  // molecular ring
  if (ring)
  {
    // inside spiral arm
    if (r > r_lim)
    {
      tmp[0] = (cos_p * (y / r) - sin_p * (x / r)) * quadruple;  // sin(t-p)
      tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) * quadruple; //-cos(t-p)
    }
    // inside molecular ring
    else
    {
      tmp[0] = (1 - 2 * bss) * y / r; // sin(phi)
      tmp[1] = (2 * bss - 1) * x / r; //-cos(phi)
    }
  }
  // elliptical bar (replace molecular ring)
  else if (bar)
  {
    const auto cos_phi = cos(p.bar_phi0);
    const auto sin_phi = sin(p.bar_phi0);
    auto new_x = cos_phi * x - sin_phi * y;
    auto new_y = sin_phi * x + cos_phi * y;
    double sgn_nx = 1.;
    double sgn_ny = 1.;
    if (new_x < 0)   // manual copysign, as autodiff runs into problems with std::copysign
      sgn_nx = -1.; 
    if (new_y != 0)
      sgn_ny = -1;
    // inside spiral arm
    if (r > bar_lim)
    {
      tmp[0] =
          (cos_p * (y / r) - sin_p * (x / r)) * quadruple; // sin(t-p)
      tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) *
               quadruple; //-cos(t-p)
    }
    // inside elliptical bar
    else
    { if (new_y!= 0) 
      {
      new_x = sgn_ny;
      new_y = -sgn_ny * (new_x / new_y) * p.bar_b * p.bar_b / (p.bar_a * p.bar_a);
      tmp[0] = (cos_phi * new_x + sin_phi * new_y) * (1 - 2 * bss);
      tmp[1] = (-sin_phi * new_x + cos_phi * new_y) * (1 - 2 * bss);
        // versor
      auto tmp_length = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
      if (tmp_length != 0.)
      {
        for (int i = 0; i < tmp.size(); ++i)
        {
          tmp[i] = tmp[i] / tmp_length;
        }
      }
      }
      else
      {
        tmp[0] = (2 * bss - 1) * sgn_nx * sin_phi;
        tmp[1] = (2 * bss - 1) * sgn_nx * cos_phi;
      }
    }
  }
  return tmp;
}

template <typename T>
T JaffeMagneticField::radial_scaling(const double &x, const double &y, const JaffeParameters<T> &p) const
{
  const double r2 = x * x + y * y;
  // separate into 3 parts for better view
  const auto s1{1. - exp(-r2 / (p.r_inner * p.r_inner))};
  const auto s2{exp(-r2 / (p.r_scale * p.r_scale))};
  const auto s3 = p.r_peak == 0 ? 1. : exp(-r2 * r2 / (p.r_peak * p.r_peak * p.r_peak * p.r_peak));
  return s1 * (s2 + s3);
}

template <typename T>
std::vector<T> JaffeMagneticField::arm_compress(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const
{
  const auto r{sqrt(x * x + y * y) / p.comp_r};
  const auto c0{1. / p.comp_c - 1.};
  std::vector<T> a0 = dist2arm(x, y, p);

  const auto r_scaling{radial_scaling(x, y, p)};
  const auto z_scaling{arm_scaling(z, p)};
  // for saving computing time
  const auto d0_inv{(r_scaling * z_scaling) / p.comp_d};
  auto factor{c0 * r_scaling * z_scaling};
  if (r > 1.)
  {
    auto cdrop{pow(r, -p.comp_p)};
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i)
    {
      a0[i] = factor * cdrop * exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
    }
  }
  else
  {
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i)
    {
      a0[i] = factor * exp(-a0[i] * a0[i] * d0_inv * d0_inv);
    }
  }
  return a0;
}

template <typename T>
std::vector<T> JaffeMagneticField::arm_compress_dust(const double &x, const double &y, const double &z, const JaffeParameters<T> &p) const
{
  const auto r{sqrt(x * x + y * y) / p.comp_r};
  const auto c0{1. / p.comp_c - 1.};
  std::vector<T> a0 = dist2arm(x, y, p);
  const auto r_scaling{radial_scaling(x, y, p)};
  const auto z_scaling{arm_scaling(z, p)};
  // only difference from normal arm_compress
  const auto d0_inv{(r_scaling) / p.comp_d};
  auto factor{c0 * r_scaling * z_scaling};
  if (r > 1)
  {
    auto cdrop{pow(r, -p.comp_p)};
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i)
    {
      a0[i] = factor * cdrop * exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
    }
  }
  else
  {
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i)
    {
      a0[i] = factor * exp(-a0[i] * a0[i] * d0_inv * d0_inv);
    }
  }
  return a0;
}

template <typename T>
std::vector<T> JaffeMagneticField::dist2arm(const double &x, const double &y, const JaffeParameters<T> &p) const
{
  const double r{sqrt(x * x + y * y)};
  const auto r_lim{p.ring_r};
  const auto bar_lim{p.bar_a + 0.5 * p.comp_d};
  auto arm_pitch = p.arm_pitch * M_PI /180;
  const auto cos_p = cos(arm_pitch);
  const auto sin_p = sin(arm_pitch); // pitch angle
  const auto beta_inv{-sin_p / cos_p};
  auto theta{atan2(y, x)};

  int arm_num = this->arm_num;

  if (ring)
  {
    if (r < r_lim)
    {
      int arm_num = 1;
    }
  }
  else if (bar)
  {
    if (r < bar_lim)
    {
      int arm_num = 1;
    }
  }

  std::vector<T> d;

  if (theta < 0)
    theta += 2 * M_PI;
  // if molecular ring
  if (ring)
  {
    // in molecular ring, return oly first element of d is used
    if (r < r_lim)
    {
      d.push_back(abs(p.ring_r - r));
    }
    // in spiral arm, return vector with arm_num elements
    else
    {
      // loop through arms
      std::vector<T> arm_phi{p.arm_phi1, p.arm_phi2, p.arm_phi3, p.arm_phi4};
      for (int i = 0; i < this->arm_num; ++i)
      {
        auto d_ang{arm_phi[i]*M_PI/180 - theta};
        auto d_rad{
            abs(p.arm_r0 * exp(d_ang * beta_inv) - r)};
        auto d_rad_p{
            abs(p.arm_r0 * exp((d_ang + 2 * M_PI) * beta_inv) - r)};
        auto d_rad_m{
            abs(p.arm_r0 * exp((d_ang - 2 * M_PI) * beta_inv) - r)};
        d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
      }
    }
  }
  // if elliptical bar
  else if (bar) {
    if (r == 0.) {
      d.push_back(0.);
    }
    else {
      const auto cos_tmp{cos(p.bar_phi0) * x / r - sin(p.bar_phi0) * y / r};
      // cos(phi)cos(phi0) - sin(phi)sin(phi0)
      const auto sin_tmp{cos(p.bar_phi0) * y / r + sin(p.bar_phi0) * x / r};
      // sin(phi)cos(phi0) + cos(phi)sin(phi0)
      // in bar, return single element vector
      if (r < bar_lim)
      {
        d.push_back(abs(p.bar_a * p.bar_b / sqrt(p.bar_a * p.bar_a * sin_tmp * sin_tmp + p.bar_b * p.bar_b * cos_tmp * cos_tmp) - r));
      }
      // in spiral arm, return vector with arm_num elements
      else
      {
        // loop through arms
        std::vector<T> arm_phi{p.arm_phi1, p.arm_phi2, p.arm_phi3, p.arm_phi4};
        for (int i = 0; i < this->arm_num; ++i)
        {
          auto d_ang{arm_phi[i]*M_PI/180 - theta};
          auto d_rad{abs(p.arm_r0 * exp(d_ang * beta_inv) - r)};
          auto d_rad_p{abs(p.arm_r0 * exp((d_ang + 2* M_PI) * beta_inv) - r)};
          auto d_rad_m{abs(p.arm_r0 * exp((d_ang - 2 * M_PI) * beta_inv) - r)};
          d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
        }
      }
    }
  }
  return d;
}

template <typename T>
T JaffeMagneticField::arm_scaling(const double &z, const JaffeParameters<T> &p) const
{
  return 1. / (cosh(z / p.arm_z0) *
               cosh(z / p.arm_z0));
}

template <typename T>
T JaffeMagneticField::disk_scaling(const double &z, const JaffeParameters<T> &p) const
{
  return 1. / (cosh(z / p.disk_z0) *
               cosh(z / p.disk_z0));
}

template <typename T>
T JaffeMagneticField::halo_scaling(const double &z, const JaffeParameters<T> &p) const
{
  return 1. / (cosh(z / p.halo_z0) *
               cosh(z / p.halo_z0));
}


IMAGINE_INSTANTIATE_VECTOR_MODEL(JaffeMagneticField)

}
