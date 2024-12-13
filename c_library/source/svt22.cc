#include <cmath>
#include <cassert>
#include <iostream>
#include "hamunits.h"
#include "SVT22.h"

vector SVT22MagneticField::_at_position(const double &x, const double &y, const double &z, const SVT22MagneticField &p) const
{
  const double r{sqrt(x * x + y * y)};
  const double rho{
      sqrt(x * x + y * y + z * z)};
  const double phi{atan2(y, x)};

    
  number B_cyl[3] = {0, 0, 0}; // the disk field in cylindrical coordinates

  //-------------------------------------------------------------------------
  ////TOROIDAL HALO COMPONENT

  if (do_halo) {
    number b1, rh;
    number B_h = 0.;
    number z_min = 0.1;
    if (z >= 0)
    { // North
      b1 = p.B_val;
    }
    else
    { // South
      b1 = -p.B_val;
    }

    B_h = b1 * (exp(-z_min/std::abs(z)) * exp(-std::abs(r) / p.r_cut) * exp(-(std::abs(z)) / (p.z_cut))); // vertical exponential fall-off
    const number B_cyl_h[3] = {0., B_h * 1, 0.};
    // add fields together
    B_cyl[0] += B_cyl_h[0];
    B_cyl[1] += B_cyl_h[1];
    B_cyl[2] += B_cyl_h[2];
  }


  // convert field to cartesian coordinates
  vector B_cart{{0.0, 0.0, 0.0}};
  B_cart[0] = B_cyl[0] * cos(phi) - B_cyl[1] * sin(phi);
  B_cart[1] = B_cyl[0] * sin(phi) + B_cyl[1] * cos(phi);
  B_cart[2] = B_cyl[2];
  return B_cart;
}

#if autodiff_FOUND

Eigen::MatrixXd SVT22MagneticField::_jac(const double &x, const double &y, const double &z, SVT22MagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, SVT22MagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.B_val, p.r_cut, p.z_cut),
                                        ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif
