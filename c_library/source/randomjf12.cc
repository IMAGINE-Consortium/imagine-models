#include <cmath>

#include "ImagineModelsRandom/RandomJF12.h"

namespace imagine {

double JF12RandomField::spectrum(const double &abs_k) const {
  return simple_spectrum(abs_k, spectral_offset, spectral_slope);
}

Vec3<double> JF12RandomField::anisotropy_direction(const double &x, const double &y, const double &z) const {
  return regular_base.at_position(x, y, z);
}

double JF12RandomField::rms(const double &x, const double &y, const double &z) const {

      const double r{sqrt(x * x + y * y)};
      const double rho{
          sqrt(x*x + y*y + z*z)};
      const double phi{atan2(y, x)};

      const double b_arms[8] = {b0_1, b0_2, b0_3, b0_4, b0_5, b0_6, b0_7, b0_8};

      double scaling_disk = 0.0;
      double scaling_halo = 0.0;

      if (r > Rmax || rho < rho_GC) {
        return 0.0;
      }
      if (r < 5.) {
        scaling_disk = b0_int;
      } else {
        double r_negx =
            r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi - M_PI));
        if (r_negx > rc_B[7]) {
          r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + M_PI));
        }
        if (r_negx > rc_B[7]) {
          r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + 3 * M_PI));
        }
        for (int i = 7; i >= 0; i--) {
          if (r_negx < rc_B[i] ) {
            scaling_disk = b_arms[i] * (5.) / r;
          }
        } // "region 8,7,6,..,2"
      }

      scaling_disk = scaling_disk * exp(-0.5 * z * z / (z0_spiral * z0_spiral));
      scaling_halo = b0_halo * exp(-r / r0_halo) * exp(-0.5 * z * z / (z0_halo * z0_halo));

      return std::sqrt(scaling_disk * scaling_disk + scaling_halo * scaling_halo);

}

}
