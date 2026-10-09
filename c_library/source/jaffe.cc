#include "ImagineModels/Jaffe.h"
#include "ImagineModels/units.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace imagine {

void JaffeMagneticField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown Jaffe model '" + model + "'.");
    active_model = model;
    parameters = JaffeParameters<double>{};
    quadruple = false;
    bss = false;
    ring = false;
    bar = true;
    arm_num = 4;
    hammurabi_v3 = false;
    r_max = 0.;
    if (model == "Jaffe13") {
        ring = true;
        bar = false;
        hammurabi_v3 = true;
        r_max = 20.;
        auto &p = parameters;
        p.disk_amp = 0.;
        p.halo_amp = 1.;
        p.halo_z0 = 6.;
        p.r_inner = 0.;
        p.r_scale = 20.;
        p.ring_amp = -0.8;
        p.ring_r = 5.;
        p.arm_r0 = 7.1;
        p.arm_z0 = 2.;
        p.arm_phi1 = 350.;
        p.arm_phi2 = 260.;
        p.arm_phi3 = 170.;
        p.arm_phi4 = 80.;
        p.arm_amp1 = 3.;
        p.arm_amp2 = 0.5;
        p.arm_amp3 = -4.;
        p.arm_amp4 = 1.2;
        p.arm_pitch = 11.5;
        p.comp_c = 1. / 3.5;
        p.comp_d = 0.3;
        p.comp_r = 12.;
        p.comp_p = 3.;
    }
}

template <typename T>
Vec3<T> JaffeMagneticField::field(const double &x, const double &y, const double &z,
                                  const JaffeParameters<T> &p) const {
    if (x == 0. && y == 0. && z == 0.) {
        return Vec3<T>{{0., 0., 0.}};
    }
    if (r_max > 0. && std::sqrt(x * x + y * y + z * z) > r_max)
        return Vec3<T>{{0., 0., 0.}};
    T inner_b{0};
    if (ring) {
        inner_b = p.ring_amp;
    } else if (bar) {
        inner_b = p.bar_amp;
    }

    Vec3<T> bhat = orientation(x, y, z, p);
    Vec3<T> btot{{0., 0., 0.}};

    auto scaling = radial_scaling(x, y, p) * (p.disk_amp * disk_scaling(z, p) + p.halo_amp * halo_scaling(z, p));
    // reversal inside negative ring
    if (hammurabi_v3 && inner_b < 0.) {
        const double r = std::sqrt(x * x + y * y);
        if ((ring && r < p.ring_r) || (!ring && bar && r < p.bar_a + 0.5 * p.comp_d))
            scaling = -scaling;
    }

    for (int i = 0; i < bhat.size(); ++i) {
        btot[i] = bhat[i] * scaling;
    }

    // compression per arm, ring or bar
    std::vector<T> arm = arm_compress(x, y, z, p);
    if (hammurabi_v3 && !arm.empty()) {
        std::array<T, 4> arm_amp = {p.arm_amp1, p.arm_amp2, p.arm_amp3, p.arm_amp4};
        for (std::size_t i = 0; i < arm.size(); ++i) {
            const T amp = i + 1 < arm.size() ? arm_amp[i] : inner_b;
            for (int j = 0; j < bhat.size(); ++j)
                btot[j] += bhat[j] * arm[i] * amp;
        }
        return btot;
    }
    // only inner region
    if (arm.size() == 1) {
        for (int i = 0; i < bhat.size(); ++i) {
            btot[i] += bhat[i] * arm[0] * inner_b;
        }
    }

    // spiral arm region
    else {
        std::array<T, 4> arm_amp = {p.arm_amp1, p.arm_amp2, p.arm_amp3, p.arm_amp4};
        for (decltype(arm.size()) i = 0; i < arm.size(); ++i) {
            for (int j = 0; j < bhat.size(); ++j) {
                btot[j] += bhat[j] * arm[i] * arm_amp[i];
            }
        }
    }
    return btot;
}

template <typename T>
Vec3<T> JaffeMagneticField::orientation(const double &x, const double &y, const double &z,
                                        const JaffeParameters<T> &p) const {
    if (x == 0. && y == 0.) {
        return Vec3<T>{{0., 0., 0.}};
    }

    const double r{std::sqrt(x * x + y * y)};
    const double r_test{hammurabi_v3 ? std::sqrt(x * x + y * y + z * z) : r};
    const auto r_lim = p.ring_r;
    const auto bar_lim{p.bar_a + 0.5 * p.comp_d};
    auto arm_pitch = p.arm_pitch * units::deg;
    const auto cos_p = cos(arm_pitch);
    const auto sin_p = sin(arm_pitch); // pitch angle

    Vec3<T> tmp{{0., 0., 0.}};
    T quadruple{1.};
    if (r_test < 0.5) // forbidden region
        return tmp;
    if (z > p.disk_z0)
        quadruple = (1 - 2 * this->quadruple);
    // molecular ring
    if (ring) {
        // inside spiral arm
        if (r_test > r_lim) {
            tmp[0] = (cos_p * (y / r) - sin_p * (x / r)) * quadruple;  // sin(t-p)
            tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) * quadruple; //-cos(t-p)
        }
        // inside molecular ring
        else {
            tmp[0] = (1 - 2 * bss) * y / r; // sin(phi)
            tmp[1] = (2 * bss - 1) * x / r; //-cos(phi)
        }
    }
    // elliptical bar, replaces ring
    else if (bar) {
        const auto cos_phi = cos(p.bar_phi0 * units::deg);
        const auto sin_phi = sin(p.bar_phi0 * units::deg);
        const auto x_rot = cos_phi * x - sin_phi * y;
        const auto y_rot = sin_phi * x + cos_phi * y;
        const double sgn_x = x_rot < 0 ? -1. : 1.;
        const double sgn_y = y_rot < 0 ? -1. : 1.;
        // inside spiral arm
        if (r_test > bar_lim) {
            tmp[0] = (cos_p * (y / r) - sin_p * (x / r)) * quadruple;  // sin(t-p)
            tmp[1] = (-cos_p * (x / r) - sin_p * (y / r)) * quadruple; //-cos(t-p)
        }
        // inside elliptical bar
        else {
            if (y_rot != 0) {
                const auto new_x = sgn_y;
                const auto new_y = -sgn_y * (x_rot / y_rot) * p.bar_b * p.bar_b / (p.bar_a * p.bar_a);
                tmp[0] = (cos_phi * new_x + sin_phi * new_y) * (1 - 2 * bss);
                tmp[1] = (-sin_phi * new_x + cos_phi * new_y) * (1 - 2 * bss);
                // versor
                auto tmp_length = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
                if (tmp_length != 0.) {
                    for (int i = 0; i < tmp.size(); ++i) {
                        tmp[i] = tmp[i] / tmp_length;
                    }
                }
            } else {
                tmp[0] = (2 * bss - 1) * sgn_x * sin_phi;
                tmp[1] = (2 * bss - 1) * sgn_x * cos_phi;
            }
        }
    }
    return tmp;
}

template <typename T>
T JaffeMagneticField::radial_scaling(const double &x, const double &y, const JaffeParameters<T> &p) const {
    const double r2 = x * x + y * y;
    const auto s1 = p.r_inner == 0. ? T(1.) : T(1. - exp(-r2 / (p.r_inner * p.r_inner)));
    const auto s2{exp(-r2 / (p.r_scale * p.r_scale))};
    const auto s3 = p.r_peak == 0 ? 0. : exp(-r2 * r2 / (p.r_peak * p.r_peak * p.r_peak * p.r_peak));
    return s1 * (s2 + s3);
}

template <typename T>
std::vector<T> JaffeMagneticField::arm_compress(const double &x, const double &y, const double &z,
                                                const JaffeParameters<T> &p) const {
    const auto r{sqrt(x * x + y * y) / p.comp_r};
    const auto c0{1. / p.comp_c - 1.};
    std::vector<T> a0 = dist2arm(x, y, p);

    const auto r_scaling{radial_scaling(x, y, p)};
    const auto z_scaling{arm_scaling(z, p)};
    // for saving computing time
    const auto d0_inv{(r_scaling * z_scaling) / p.comp_d};
    auto factor{c0 * r_scaling * z_scaling};
    if (r > 1.) {
        auto cdrop{pow(r, -p.comp_p)};
        for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
            a0[i] = factor * cdrop * exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
        }
    } else {
        for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
            a0[i] = factor * exp(-a0[i] * a0[i] * d0_inv * d0_inv);
        }
    }
    return a0;
}

template <typename T>
std::vector<T> JaffeMagneticField::arm_compress_dust(const double &x, const double &y, const double &z,
                                                     const JaffeParameters<T> &p) const {
    const auto r{sqrt(x * x + y * y) / p.comp_r};
    const auto c0{1. / p.comp_c - 1.};
    std::vector<T> a0 = dist2arm(x, y, p);
    const auto r_scaling{radial_scaling(x, y, p)};
    const auto z_scaling{arm_scaling(z, p)};
    // differs from arm_compress
    const auto d0_inv{(r_scaling) / p.comp_d};
    auto factor{c0 * r_scaling * z_scaling};
    if (r > 1) {
        auto cdrop{pow(r, -p.comp_p)};
        for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
            a0[i] = factor * cdrop * exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
        }
    } else {
        for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
            a0[i] = factor * exp(-a0[i] * a0[i] * d0_inv * d0_inv);
        }
    }
    return a0;
}

template <typename T>
std::vector<T> JaffeMagneticField::dist2arm(const double &x, const double &y, const JaffeParameters<T> &p) const {
    const double r{sqrt(x * x + y * y)};
    const auto r_lim{p.ring_r};
    const auto bar_lim{p.bar_a + 0.5 * p.comp_d};
    auto arm_pitch = p.arm_pitch * units::deg;
    const auto cos_p = cos(arm_pitch);
    const auto sin_p = sin(arm_pitch); // pitch angle
    const auto beta_inv{-sin_p / cos_p};
    auto theta{atan2(y, x)};

    if (arm_num < 2 or arm_num > 4)
        throw std::invalid_argument("JaffeMagneticField: arm_num must be between 2 and 4.");

    std::vector<T> d;
    // distance to arm, both conventions
    auto arm_distance = [&](const T &d_ang) -> T {
        if (hammurabi_v3) {
            T best = r;
            for (int k = -4; k <= 4; ++k) {
                const T candidate = abs(p.arm_r0 * exp((d_ang + 2 * k * units::pi) * beta_inv) - r);
                if (candidate < best)
                    best = candidate;
            }
            return best;
        }
        const T d_rad = abs(p.arm_r0 * exp(d_ang * beta_inv) - r);
        const T d_rad_p = abs(p.arm_r0 * exp((d_ang + 2 * units::pi) * beta_inv) - r);
        const T d_rad_m = abs(p.arm_r0 * exp((d_ang - 2 * units::pi) * beta_inv) - r);
        return std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p;
    };

    if (theta < 0)
        theta += 2 * units::pi;
    // if molecular ring
    if (ring) {
        // ring: first element only
        if (r < r_lim) {
            d.push_back(abs(p.ring_r - r));
        }
        // arms: arm_num elements
        else {
            // loop through arms
            std::vector<T> arm_phi{p.arm_phi1, p.arm_phi2, p.arm_phi3, p.arm_phi4};
            for (int i = 0; i < this->arm_num; ++i) {
                d.push_back(arm_distance(T(arm_phi[i] * units::deg - theta)));
            }
        }
    }
    // if elliptical bar
    else if (bar) {
        if (r == 0.) {
            d.push_back(0.);
        } else {
            const auto cos_tmp{cos(p.bar_phi0 * units::deg) * x / r - sin(p.bar_phi0 * units::deg) * y / r};
            // cos(phi)cos(phi0) - sin(phi)sin(phi0)
            const auto sin_tmp{cos(p.bar_phi0 * units::deg) * y / r + sin(p.bar_phi0 * units::deg) * x / r};
            // sin(phi)cos(phi0) + cos(phi)sin(phi0)
            // bar: single element
            if (r < bar_lim) {
                d.push_back(
                    abs(p.bar_a * p.bar_b /
                            sqrt(p.bar_a * p.bar_a * sin_tmp * sin_tmp + p.bar_b * p.bar_b * cos_tmp * cos_tmp) -
                        r));
            }
            // arms: arm_num elements
            else {
                // loop through arms
                std::vector<T> arm_phi{p.arm_phi1, p.arm_phi2, p.arm_phi3, p.arm_phi4};
                for (int i = 0; i < this->arm_num; ++i) {
                    d.push_back(arm_distance(T(arm_phi[i] * units::deg - theta)));
                }
            }
        }
    }
    // inactive components at 100 kpc
    if (hammurabi_v3 && !d.empty()) {
        const T inactive = 100.;
        if (d.size() == 1)
            d.insert(d.begin(), arm_num, inactive);
        else
            d.push_back(inactive);
    }
    return d;
}

template <typename T> T JaffeMagneticField::arm_scaling(const double &z, const JaffeParameters<T> &p) const {
    return 1. / (cosh(z / p.arm_z0) * cosh(z / p.arm_z0));
}

template <typename T> T JaffeMagneticField::disk_scaling(const double &z, const JaffeParameters<T> &p) const {
    return 1. / (cosh(z / p.disk_z0) * cosh(z / p.disk_z0));
}

template <typename T> T JaffeMagneticField::halo_scaling(const double &z, const JaffeParameters<T> &p) const {
    return 1. / (cosh(z / p.halo_z0) * cosh(z / p.halo_z0));
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(JaffeMagneticField)

// used by JaffeRandomField
template Vec3<double> JaffeMagneticField::orientation<double>(const double &, const double &, const double &,
                                                              const JaffeParameters<double> &) const;
template double JaffeMagneticField::radial_scaling<double>(const double &, const double &,
                                                           const JaffeParameters<double> &) const;
template std::vector<double> JaffeMagneticField::arm_compress<double>(const double &, const double &, const double &,
                                                                      const JaffeParameters<double> &) const;

}
