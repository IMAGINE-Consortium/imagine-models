#include <cmath>
#include "hamunits.h"
#include "Sun.h"

#include "helpers.h"

// Sun et al. A&A V.477 2008 ASS+RING model magnetic field
vector SunMagneticField::_at_position(const double &x, const double &y, const double &z, const SunMagneticField &p) const
{

  double r = sqrt(x * x + y * y);
  double phi = atan2(y, x);

  // first we set D2 (ASS+RING model, eq. 8)
  // ------------------------------------------------------------
  double D2;
  if (r > 7.5)
  {
    D2 = 1.;
  }
  else if (r <= 7.5 && r > 6.)
  {
    D2 = -1.;
  }
  else if (r <= 6. && r > 5.)
  {
    D2 = 1.;
  }
  else // if(r <= 5.)
  {
    D2 = -1.;
  }
  // ------------------------------------------------------------

  // now we set D1 (eq. 7)
  // ------------------------------------------------------------
  number D1;
  if (r > p.b_Rc)
  {
    D1 = p.b_B0 * exp(-((r - p.b_Rsun) / p.b_R0) - (std::abs(z) / p.b_z0));
  }
  else // if(r <= b_Rc)
  {
    D1 = p.b_Bc;
  }
  // ------------------------------------------------------------

  auto p_ang = p.b_p * M_PI / 180.;
  vector B_cyl{{D1 * D2 * sin(p_ang),  // eq. 6
                -D1 * D2 * cos(p_ang),
                0.}};

  // [ORIGINAL HAMMURABI COMMENT]  Taking into account the halo field
  number halo_field;

  // [ORIGINAL HAMMURABI COMMENT]  for better overview
  number b3H_z1_actual;
  if (std::abs(z) < p.bH_z0)
  {
    b3H_z1_actual = p.bH_z1a;
  }
  else
  {
    b3H_z1_actual = p.bH_z1b;
  }
  auto hf_piece1 = (b3H_z1_actual * b3H_z1_actual) / (b3H_z1_actual * b3H_z1_actual + (std::abs(z) - p.bH_z0) * (std::abs(z) - p.bH_z0));
  auto hf_piece2 = exp(-(r - p.bH_R0) / (p.bH_R0));

  halo_field = p.bH_B0 * hf_piece1 * (r / p.bH_R0) * hf_piece2;  // eq. 10

  // [ORIGINAL HAMMURABI COMMENT] Flip north.  Not sure how Sun did this. This is his code with no
  // flip though the paper says it's flipped but without this mod,
  // there is no antisymmetry across the disk.  However, it doesn't seem to work.
  if (z > 0 && r >= 5.)
  {
    halo_field *= -1.;
  }

  B_cyl[1] += halo_field;

  vector B_vec3;

  B_vec3 = Cyl2Cart<vector>(phi, B_cyl);

  return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd SunMagneticField::_jac(const double &x, const double &y, const double &z, SunMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, SunMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.b_B0, p.b_Bc, p.b_R0, p.b_Rc, p.b_z0, p.b_Rsun, p.bH_B0, p.bH_R0, p.bH_z0, p.bH_z1a, p.bH_z1b), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif