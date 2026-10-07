#include <cmath>
#include "ImagineModels/units.h"
#include "ImagineModels/HarariMollerachRoulet.h"

#include "ImagineModels/helpers.h"

namespace imagine {

template <typename T>
Vec3<T> HMRMagneticField::field(const double &x, const double &y, const double &z, const HMRParameters<T> &p) const
{

  Vec3<T> B_vec3{{0, 0, 0}};

  double r = std::sqrt(x * x + y * y);
  const double phi = std::atan2(y, x);

  auto f_z = (1. / (2. * cosh(z / p.b_z1))) + (1. / (2. * cosh(z / p.b_z2)));

  if (r < 0.0000000005)
  {
    r = 0.5;
  }

  auto b_r = (3. * p.b_Rsun / r) * tanh(r / p.b_r1) * tanh(r / p.b_r1) * tanh(r / p.b_r1);

  // BSS model (eq. 2.2 of https://arxiv.org/abs/astro-ph/9906309)
  const double theta = M_PI - phi;
  auto B_r_phi = b_r * cos(theta - ((1. / tan(p.b_p * (M_PI / 180.))) * log(r / p.b_epsilon0)));

  // B-field in cylindrical coordinates:
  Vec3<T> B_cyl{{B_r_phi * sin(p.b_p * (M_PI / 180.)) * f_z,
                -B_r_phi * cos(p.b_p * (M_PI / 180.)) * f_z,
                0.}};

  B_vec3 = Cyl2Cart<Vec3<T>>(phi, B_cyl);

  return B_vec3;
}


IMAGINE_INSTANTIATE_VECTOR_MODEL(HMRMagneticField)

}
