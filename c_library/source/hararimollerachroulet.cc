#include <cmath>
#include "../headers/hamunits.h"
#include "../headers/HarariMollerachRoulet.h"

#include "../headers/helpers.h"

std::array<double, 3>  HMRMagneticField::at_position(const double &x, const double &y, const double &z) const { 


   double r = std::sqrt(x*x + y*y);
   double phi = std::atan2(y, x);

   double f_z = (1. / (2.* std::cosh(z / b_z1))) + (1. / (2. * std::cosh(z / b_z2)));

   if (r < 0.0000000005) {
      r = 0.5;
      }

  double b_r = ( 3. * b_Rsun / r) * std::tanh(r / b_r1) * std::tanh(r / b_r1) * std::tanh(r / b_r1);

  double B_r_phi = b_r * std::cos(phi - ((1. / std::tan(b_p)) * std::log(r / b_epsilon0)));

  // B-field in cylindrical coordinates:
  std::array<double, 3> B_cyl{B_r_phi * std::sin(b_p) * f_z, 
                              B_r_phi * std::cos(b_p) * f_z , 
                              0.};

  std::array<double, 3> B_vec3{0, 0, 0};
  if (r <= b_r_max) {
    Cyl2Cart(phi, B_cyl, B_vec3);
  }
  
  return B_vec3;
}