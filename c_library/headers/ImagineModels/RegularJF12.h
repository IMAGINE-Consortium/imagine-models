#ifndef REGULARJF12_H
#define REGULARJF12_H

#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "RegularField.h"

class JF12MagneticField : public RegularVectorField
{
protected:
  vector _at_position(const double &x, const double &y, const double &z, const JF12MagneticField &p) const;

#if autodiff_FOUND
  Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, JF12MagneticField &p) const;
#endif

public:
  using RegularVectorField ::RegularVectorField;

  // define fixed parameters
  const double Rmax = 20;   // outer boundary of GMF
  const double rho_GC = 1.; // interior boundary of GMF

  // fixed disk parameters
  const double inc = 11.5; // inclination, in degrees
  const double rmin = 5.;  // outer boundary of the molecular ring region
  const double rcent = 3.; // inner boundary of the molecular ring region (field is
                           // zero within this region)
  const number f[8] = {
      0.130, 0.165, 0.094, 0.122,
      0.13, 0.118, 0.084, 0.156}; // fractions of circumference spanned by each
                                  // spiral, sums to unity
  const double rc_B[8] = {
      5.1, 6.3, 7.1, 8.3, 9.8,
      11.4, 12.7, 15.5}; // the radii where the spiral arm boundaries cross the
                         // negative x-axis

  number b_arm_1 = 0.1;
  number b_arm_2 = 3.0;
  number b_arm_3 = -0.9;
  number b_arm_4 = -0.8;
  number b_arm_5 = -2.0;
  number b_arm_6 = -4.2;
  number b_arm_7 = 0.0;
  number b_ring = 0.1;
  number h_disk = 0.40;
  number w_disk = 0.27;
  // toroidal halo parameters
  bool do_halo = true;
  number Bn = 1.4;
  number Bs = -1.1;
  number rn = 9.22;
  number rs = 16.7;
  number wh = 0.20;
  number z0 = 5.3;
  // X-field parameters
  bool do_X = true;
  number B0_X = 4.6;
  number Xtheta_const = 49.;
  number rpc_X = 4.8;
  number r0_X = 2.9;
  
#if autodiff_FOUND
  const std::set<std::string> all_diff{"b_arm_1", "b_arm_2", "b_arm_3", "b_arm_4", "b_arm_5", "b_arm_6", "b_arm_7", "b_ring", "h_disk", "w_disk", "Bn", "Bs", "rn", "rs", "wh", "z0", "B0_X", "Xtheta_const", "rpc_X", "r0_X"};
  std::set<std::string> active_diff{"b_arm_1", "b_arm_2", "b_arm_3", "b_arm_4", "b_arm_5", "b_arm_6", "b_arm_7", "b_ring", "h_disk", "w_disk", "Bn", "Bs", "rn", "rs", "wh", "z0", "B0_X", "Xtheta_const", "rpc_X", "r0_X"};

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
