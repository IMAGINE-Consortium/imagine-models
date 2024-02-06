#include <cmath>
#include "hamunits.h"
#include "HarariMollerachRoulet.h"

#include "helpers.h"

vector HMRMagneticField::_at_position(const double &x, const double &y, const double &z, const HMRMagneticField &p) const
{

  vector B_vec3{{0, 0, 0}};

  double r = std::sqrt(x * x + y * y);
  const double phi = std::atan2(y, x);

  auto f_z = (1. / (2. * cosh(z / p.b_z1))) + (1. / (2. * cosh(z / p.b_z2)));

  if (r < 0.0000000005)
  {
    r = 0.5;
  }

  auto b_r = (3. * p.b_Rsun / r) * tanh(r / p.b_r1) * tanh(r / p.b_r1) * tanh(r / p.b_r1);

  // BSS model (eq. 2.2 of https://arxiv.org/abs/astro-ph/9906309)
  auto B_r_phi = b_r * cos(-phi - ((1. / tan(p.b_p * (M_PI / 180.))) * log(r / p.b_epsilon0)));

  // B-field in cylindrical coordinates:
  vector B_cyl{{B_r_phi * sin(p.b_p * (M_PI / 180.)) * f_z,
                B_r_phi * cos(p.b_p * (M_PI / 180.)) * f_z,
                0.}};

  B_vec3 = Cyl2Cart<vector>(phi, B_cyl);

  return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd HMRMagneticField::_jac(const double &x, const double &y, const double &z, HMRMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, HMRMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.b_Rsun, p.b_z1, p.b_z2, p.b_r1, p.b_p, p.b_epsilon0), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif