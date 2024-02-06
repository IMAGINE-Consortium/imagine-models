#include "Han.h"

#include "helpers.h"

// J. L. Han et al 2018 ApJS 234 11
vector HanMagneticField::_at_position(const double &x, const double &y, const double &z, const HanMagneticField &p) const
{
  vector B_cyl{{0., 0., 0.}};
  const double r = sqrt(x * x + y * y);
  const double phi = atan2(y, x);

  if (r < p.R_min || r > p.R_max)
    return B_cyl;

  number B_0 = 0.;

  auto p_ang = p.B_p * M_PI / 180.;
  const double phi_han = -(phi + M_PI); // nneeded to fix different coordinate system convention
  
  number R_0 = r * exp(phi_han * tan(p_ang));  // eq. 4 is wrong, need to change psi and phi!

  std::array<number, 6> B_s = {p.B_s1, p.B_s2, p.B_s3, p.B_s4, p.B_s5, p.B_s6};  // table 5

  if (R_0 < p.R_s[0])
  {
    R_0 = r * exp((phi_han + 2 * M_PI) * tan(p_ang));  // eq. 4
  }
  if (R_0 > p.R_s[6])
  {
    R_0 = r * exp((phi_han - 2 * M_PI) * tan(p_ang));  // eq. 4
  }

  for (int i = 0; i < 6; i++)
  {
    if (p.R_s[i] < R_0)
    { 
      if (R_0 < p.R_s[i+1])
      {
        B_0 = B_s[i];
        break;
      }
    }
  }

  number B_r = B_0 * exp(-r / p.A) * exp(-std::abs(z) / p.H);  // eq. 3

  B_cyl[0] = B_r * sin(p_ang);
  B_cyl[1] = B_r * cos(p_ang);

  vector B_vec3 = Cyl2Cart<vector>(phi, B_cyl);
  return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd HanMagneticField::_jac(const double &x, const double &y, const double &z, HanMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, HanMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.B_p, p.A, p.H, p.B_s1, p.B_s2, p.B_s3, p.B_s4, p.B_s4, p.B_s5, p.B_s6), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif
