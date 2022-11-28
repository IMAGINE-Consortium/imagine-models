#include <cmath>
#include <vector>
#include <iostream>


template<typename G>
class JaffeMagneticField : public RegularField<G, std::vector<double>>  {
    protected:
      bool DEBUG = false;
    public:

    using RegularField<G, std::vector<double>> :: RegularField;

      bool quadruple = false; // quadruple pattern in halo
      bool bss = false; // bi-symmetric

      double disk_amp = 0.167;  // disk amplitude, microG
      double disk_z0 = 0.1; // disk height scale, kpc
      double halo_amp = 1.38; // halo amplitude, microG
      double halo_z0 = 3.0; // halo height scale, kpc
      double r_inner = 0.5; // inner R scale, kpc
      double r_scale = 20.; // R scale, kpc
      double r_peak = 0.; // R peak, kpc

      bool ring = false; // molecular ring
      bool bar = true; // elliptical bar
       // either ring or bar!
      double ring_amp = 0.023; // ring field amplitude, microG
      double ring_r = 5.0; // ring radius, kpc
      double bar_amp = 0.023; // bar field amplitude, microG
      double bar_a = 5.0; // major scale, kpc
      double bar_b = 3.0; // minor scale, kpc
      double bar_phi0 = 45.0; // bar major direction

      int arm_num = 4; // # of spiral arms
      double arm_r0 = 7.1; // arm ref radius, kpc
      double arm_z0 = 0.1; // arm heigth scale, kpc
      double arm_phi1 = 70; //arm ref angles, deg
      double arm_phi2 = 160;
      double arm_phi3 = 250;
      double arm_phi4 = 340;
      double arm_amp1 = 2; //  arm field amplitudes, microG
      double arm_amp2 = 0.133;
      double arm_amp3 = -3.78;
      double arm_amp4 = 0.32;
      double arm_pitch = 11.5; // pitch angle, deg

      double comp_c = 0.5; // compress factor
      double comp_d = 0.3; // arm cross-sec scale, kpc
      double comp_r = 12; // radial cutoff scale, kpc
      double comp_p = 3; //cutoff power

      std::vector<double> at_position(const double &x, const double &y, const double &z) const {
          double inner_b{0};
          if (ring) {
            inner_b = ring_amp;}
          else if (bar) {
            inner_b = bar_amp;}
          std::vector<double> bhat = orientation(x, y, z);
          std::vector<double> btot{0., 0., 0.};
          double scaling = radial_scaling(x, y) *
                 (disk_amp * disk_scaling(z) +
                  halo_amp * halo_scaling(z));
          for(int i=0;i<bhat.size();++i) {
              	btot[i] = bhat[i] * scaling;}
          // compress factor for each arm or for ring/bar
          std::vector<double> arm;
          arm = arm_compress(x, y, z);
          // only inner region
          if (arm.size() == 1) {
            for(int i=0;i<bhat.size();++i) {
              btot[i] += bhat[i] * arm[0] * inner_b;}
          }

          // spiral arm region
          else {
            std::vector<double> arm_amp{arm_amp1, arm_amp2, arm_amp3, arm_amp4};
            for (decltype(arm.size()) i = 0; i < arm.size(); ++i) {
              for(int j=0;j<bhat.size();++j) {
                btot[j] += bhat[j] * arm[i] * arm_amp[i];}
            }

          }
          return btot;
        }

    std::vector<double> orientation(const double &x, const double &y, const double &z) const {
          const double r{
              sqrt(x * x + y * y)}; // cylindrical frame
          const double r_lim = ring_r;
          const double bar_lim{bar_a + 0.5 * comp_d};
          const double cos_p = std::cos(arm_pitch);
          const double sin_p = std::sin(arm_pitch); // pitch angle
          std::vector<double> tmp{0., 0., 0.};
          double quadruple{1};
          if (r < 0.5) // forbiden region
            return tmp;
          if (z > disk_z0)
            quadruple = (1 - 2 * quadruple);
          // molecular ring
          if (ring) {
            // inside spiral arm
            if (r > r_lim) {
              tmp[0] =
                  (cos_p * (y / r) - sin_p * (x / r)) * quadruple; // sin(t-p)
              tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) *
                       quadruple; //-cos(t-p)
            }
            // inside molecular ring
            else {
              tmp[0] = (1 - 2 * bss) * y / r; // sin(phi)
              tmp[1] = (2 * bss - 1) * x / r; //-cos(phi)
            }
          }
          // elliptical bar (replace molecular ring)
          else if (bar) {
            const double cos_phi = std::cos(bar_phi0);
            const double sin_phi = std::sin(bar_phi0);
            const double x = cos_phi * x - sin_phi * y;
            const double y = sin_phi * x + cos_phi * y;
            // inside spiral arm
            if (r > bar_lim) {
              tmp[0] =
                  (cos_p * (y / r) - sin_p * (x / r)) * quadruple; // sin(t-p)
              tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) *
                       quadruple; //-cos(t-p)
            }
            // inside elliptical bar
            else {
              if (y != 0) {
                const double new_x = copysign(1, y);
                const double new_y = copysign(1, y) * (x / y) *
                                      bar_b * bar_b /
                                      (bar_a * bar_a);
                tmp[0] =
                    (cos_phi * new_x + sin_phi * new_y) * (1 - 2 * bss);
                tmp[1] = (-sin_phi * new_x + cos_phi * new_y) *
                         (1 - 2 * bss);
                // versor
                double tmp_length = std::sqrt(tmp[0]* tmp[0] + tmp[1]* tmp[1] + tmp[2]* tmp[2]);
                if (tmp_length != 0.) {
                  for(int i=0;i<tmp.size();++i) {
                    tmp[i] = tmp[i]/tmp_length;}
                  }
              } else {
                tmp[0] = (2 * bss - 1) * copysign(1, x) * sin_phi;
                tmp[1] = (2 * bss - 1) * copysign(1, x) * cos_phi;
              }
            }
          }
          return tmp;
        }

    double radial_scaling(const double &x, const double &y) const {
          const double r2 = x * x + y * y;
          // separate into 3 parts for better view
          const double s1{
              1. - std::exp(-r2 / (r_inner * r_inner))};
          const double s2{
              std::exp(-r2 / (r_scale * r_scale))};
          const double s3{
              std::exp(-r2 * r2 /
                       (r_peak * r_peak *
                        r_peak * r_peak))};
          return s1 * (s2 + s3);
        }

    std::vector<double> arm_compress(const double &x, const double &y,  const double &z) const {
          const double r{sqrt(x * x + y * y) /
                            comp_r};
          const double c0{1. / comp_c - 1.};
          std::vector<double> a0 = dist2arm(x, y);
          const double r_scaling{radial_scaling(x, y)};
          const double z_scaling{arm_scaling(z)};
          // for saving computing time
          const double d0_inv{(r_scaling * z_scaling) / comp_d};
          double factor{c0 * r_scaling * z_scaling};
          if (r > 1) {
            double cdrop{std::pow(r, -comp_p)};
            for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
              a0[i] = factor * cdrop *
                      std::exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
            }
          } else {
            for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
              a0[i] = factor * std::exp(-a0[i] * a0[i] * d0_inv * d0_inv);
            }
          }
          return a0;
        }

    std::vector<double> arm_compress_dust(const double &x, const double &y, const double &z) const {
          const double r{sqrt(x * x + y * y) /
                            comp_r};
          const double c0{1. / comp_c - 1.};
          std::vector<double> a0 = dist2arm(x, y);
          const double r_scaling{radial_scaling(z)};
          const double z_scaling{arm_scaling(z)};
          // only difference from normal arm_compress
          const double d0_inv{(r_scaling) / comp_d};
          double factor{c0 * r_scaling * z_scaling};
          if (r > 1) {
            double cdrop{std::pow(r, -comp_p)};
            for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
              a0[i] = factor * cdrop *
                      std::exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
            }
          } else {
            for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
              a0[i] = factor * std::exp(-a0[i] * a0[i] * d0_inv * d0_inv);
            }
          }
          return a0;
        }

    std::vector<double> dist2arm(const double &x, const double &y) const {
          std::vector<double> d;
          const double r{sqrt(x * x + y * y)};
          const double r_lim{ring_r};
          const double bar_lim{bar_a + 0.5 * comp_d};
          const double cos_p{std::cos(arm_pitch)};
          const double sin_p{std::sin(arm_pitch)}; // pitch angle
          const double beta_inv{-sin_p / cos_p};
          double theta{atan2(y, x)};
          if (theta < 0)
            theta += 2 * cgs::pi;
          // if molecular ring
          if (ring) {
            // in molecular ring, return single element vector
            if (r < r_lim) {
              d.push_back(std::fabs(ring_r - r));
            }
            // in spiral arm, return vector with arm_num elements
            else {
              // loop through arms
              std::vector<double> arm_phi{arm_phi1, arm_phi2, arm_phi3, arm_phi4};
              for (int i = 0; i < arm_num; ++i) {
                double d_ang{arm_phi[i] - theta};
                double d_rad{
                    std::fabs(arm_r0 * std::exp(d_ang * beta_inv) - r)};
                double d_rad_p{
                    std::fabs(arm_r0 *
                                  std::exp((d_ang + 2 * cgs::pi) * beta_inv) -
                              r)};
                double d_rad_m{
                    std::fabs(arm_r0 *
                                  std::exp((d_ang - 2 * cgs::pi) * beta_inv) -
                              r)};
                d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
              }
            }
          }
          // if elliptical bar
          else if (bar) {
            const double cos_tmp{std::cos(bar_phi0) * x / r -
                                    std::sin(bar_phi0) * y /
                                        r}; // cos(phi)cos(phi0) - sin(phi)sin(phi0)
            const double sin_tmp{std::cos(bar_phi0) * y / r +
                                    std::sin(bar_phi0) * x /
                                        r}; // sin(phi)cos(phi0) + cos(phi)sin(phi0)
            // in bar, return single element vector
            if (r < bar_lim) {
              d.push_back(
                  std::fabs(bar_a * bar_b /
                                sqrt(bar_a * bar_a *
                                         sin_tmp * sin_tmp +
                                     bar_b * bar_b *
                                         cos_tmp * cos_tmp) -
                            r));
            }
            // in spiral arm, return vector with arm_num elements
            else {
              // loop through arms
              std::vector<double> arm_phi{arm_phi1, arm_phi2, arm_phi3, arm_phi4};
              for (int i = 0; i < arm_num; ++i) {
                double d_ang{arm_phi[i] - theta};
                double d_rad{
                    std::fabs(arm_r0 * std::exp(d_ang * beta_inv) - r)};
                double d_rad_p{
                    std::fabs(arm_r0 *
                                  std::exp((d_ang + 2 * cgs::pi) * beta_inv) -
                              r)};
                double d_rad_m{
                    std::fabs(arm_r0 *
                                  std::exp((d_ang - 2 * cgs::pi) * beta_inv) -
                              r)};
                d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
              }
            }
          }
          return d;
        }

      double arm_scaling(const double &z) const {
          return 1. / (std::cosh(z / arm_z0) *
                       std::cosh(z / arm_z0));
        }

      double disk_scaling(const double &z) const {
          return 1. / (std::cosh(z / disk_z0) *
                       std::cosh(z / disk_z0));
        }

      double halo_scaling(const double &z) const {
          return 1. / (std::cosh(z / halo_z0) *
                       std::cosh(z / halo_z0));
        }
/*
      std::vector<double> orientation(const double &x, const double &y,  const double &z) const;
      double radial_scaling(const double &x, const double &y) const;
      double arm_scaling(const double &z) const;
      double disk_scaling(const double &z) const;
      double halo_scaling(const double &z) const;
      std::vector<double> arm_compress(const double &x, const double &y,  const double &z) const;
      std::vector<double> arm_compress_dust(const double &x, const double &y,  const double &z) const;
      std::vector<double> dist2arm(const double &x, const double &y) const;
*/
};
