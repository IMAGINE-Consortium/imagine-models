#ifndef SVT22_H
#define SVT22_H

#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "RegularField.h"

class SVT22MagneticField : public RegularVectorField
{
protected:
  vector _at_position(const double &x, const double &y, const double &z, const SVT22MagneticField &p) const;

#if autodiff_FOUND
  Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, SVT22MagneticField &p) const;
#endif

public:
  using RegularVectorField ::RegularVectorField;

  bool do_halo = true;


//// SVT22 model
    number B_val = 3.72;
    number r_cut = 5;
    number z_cut = 6;
 
#if autodiff_FOUND

  const std::set<std::string> all_diff{"B_val", "r_cut", "z_cut"};
  std::set<std::string> active_diff{"B_val", "r_cut", "z_cut"};

  Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
  {
    return _jac(x, y, z, *this);
  }
#endif

  vector at_position(const double &x, const double &y, const double &z) const
  {
    return _at_position(x, y, z, *this);
  }
};

#endif
