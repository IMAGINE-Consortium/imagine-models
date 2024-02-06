#include <cmath>
#include "hamunits.h"
#include "Helix.h"

vector HelixMagneticField::_at_position(const double &x, const double &y, const double &z, const HelixMagneticField &p) const
{

  const double phi = std::atan2(y, x);         // azimuthal angle in cylindrical coordinates
  const double r = std::sqrt(x * x + y * y); // radius in cylindrical coordinates
  vector b{{0.0, 0.0, 0.0}};
  if ((r > rmin) && (r < rmax))
  {
    b[0] = std::cos(phi) * p.ampx;
    b[1] = std::sin(phi) * p.ampy;
    b[2] = p.ampz;
  }
  return b;
}

#if autodiff_FOUND

Eigen::MatrixXd HelixMagneticField::_jac(const double &x, const double &y, const double &z, HelixMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, HelixMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.ampx, p.ampy, p.ampz), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif
