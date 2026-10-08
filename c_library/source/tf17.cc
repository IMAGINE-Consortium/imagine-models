#include "ImagineModels/units.h"
#include <algorithm>
#include <stdexcept>

#include "ImagineModels/TF17.h"
#include "ImagineModels/helpers.h"

namespace imagine {

template <typename T>
Vec3<T> TF17MagneticField::field(const double &x, const double &y, const double &z, const TF17Parameters<T> &p) const {
    const double r = sqrt(x * x + y * y);
    double phi = units::pi - std::atan2(y, x);

    double cosPhi = cos(phi);
    double sinPhi = sin(phi);

    Vec3<T> df = getDiskField(r, z, phi, sinPhi, cosPhi, p);
    Vec3<T> hf = getHaloField(r, z, phi, sinPhi, cosPhi, p);

    Vec3<T> B_cart = addVector<Vec3<T>>({df, hf});
    return B_cart;
}

template <typename T>
Vec3<T> TF17MagneticField::getDiskField(const double &r, const double &z, const double &phi, const double &sinPhi,
                                        const double &cosPhi, const TF17Parameters<T> &p) const {
    Vec3<T> B_cart{{0., 0., 0.}};
    T B_r = 0;
    T B_phi = 0;
    T B_z = 0;

    auto psd = p.phi_star_disk * units::deg;
    auto psh = p.phi_star_halo * units::deg;
    auto p_0 = p.p_0 * units::deg;
    auto cot_p0 = cos(p_0) / sin(p_0);

    if (active_disk_model == "Ad1") {
        if (r > p.r1_disk) {
            auto z1_disk_z = (1. + p.a_disk * p.r1_disk * p.r1_disk) / (1. + p.a_disk * r * r); // z1_disk / z, eq. 4
            // cylindrical components
            auto B_r0 = radialFieldScale(p.B1_disk, psd, z1_disk_z * z, phi, r, z, cot_p0, p);
            B_r = (p.r1_disk / r) * z1_disk_z * B_r0;
            B_z = 2 * p.a_disk * p.r1_disk * z1_disk_z * z / (1 + p.a_disk * r * r) * B_r0;
            B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);
        } else {
            // straight field lines inside r1_disk
            auto phi1_disk = shiftedWindingFunction<T>(p.r1_disk, z, cot_p0, p) + psd;
            auto B_amp = p.B1_disk * exp(-fabs(z) / p.H_disk);
            B_r = cos(phi1_disk - phi) * B_amp;
            B_phi = sin(phi1_disk - phi) * B_amp;
        }
    } else if (active_disk_model == "Bd1") {
        // Bd1: n = 2
        if (r > epsilon) {
            auto r1_disk_r = p.r1_disk / r;
            auto z1_disk_z = 5. / (r1_disk_r * r1_disk_r + 4. / sqrt(r1_disk_r)); // z1_disk / z
            auto B_r0 = radialFieldScale(p.B1_disk, psd, z1_disk_z * z, phi, r, z, cot_p0, p);
            B_r = r1_disk_r * z1_disk_z * B_r0;
            B_z = -0.4 * r1_disk_r / r * z1_disk_z * z1_disk_z * z * (r1_disk_r * r1_disk_r - 1. / sqrt(r1_disk_r)) *
                  B_r0;
        } else {
            auto z1_disk_z = 5. * r * r / (p.r1_disk * p.r1_disk); // z1_disk / z
            auto B_r0 = radialFieldScale(p.B1_disk, psd, z1_disk_z * z, phi, r, z, cot_p0, p);
            B_r = 5. * r / p.r1_disk * B_r0;
            B_z = -10. * z / p.r1_disk * B_r0;
        }
        B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);
    } else if (active_disk_model == "Dd1") {
        // Dd1: n = 0.5
        double z_sign = z >= 0 ? 1. : -1.;
        double z_abs = fabs(z);
        if (z_abs > epsilon) {
            auto z1_disk_z = p.z1_disk / z_abs;
            auto r1_disk_r = 1.5 / (sqrt(z1_disk_z) + 0.5 / z1_disk_z); // r1_disk / r
            auto F_r = r1_disk_r * r <= p.L_disk ? 1. : exp(1. - r1_disk_r * r / p.L_disk);
            // cosine argument
            auto B_z0 = z_sign * p.B1_disk * F_r * cos(phi - shiftedWindingFunction<T>(r, z, cot_p0, p) - psd);
            B_r = -0.5 / 1.5 * r1_disk_r * r1_disk_r * r1_disk_r * r / z_abs * (sqrt(z1_disk_z) - 1 / z1_disk_z) * B_z0;
            B_z = z_sign * r1_disk_r * r1_disk_r * B_z0;
        } else {
            auto z_z1_disk = z_abs / p.z1_disk;
            auto r1_disk_r = 1.5 * sqrt(z_abs / p.z1_disk); // r1_disk / r
            auto F_r = r1_disk_r * r <= p.L_disk ? 1. : exp(1. - r1_disk_r * r / p.L_disk);
            auto B_z0 = z_sign * p.B1_disk * F_r * cos(phi - shiftedWindingFunction<T>(r, z, cot_p0, p) - psd);
            B_r = -1.125 * r / p.z1_disk * (1 - 2.5 * z_z1_disk * sqrt(z_z1_disk)) * B_z0;
            B_z = z_sign * r1_disk_r * r1_disk_r * B_z0;
        }
        B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);
    }

    // Cartesian components
    B_cart[0] = -(B_r * cosPhi - B_phi * sinPhi); // flip x-component at the end
    B_cart[1] = B_r * sinPhi + B_phi * cosPhi;
    B_cart[2] = B_z;
    return B_cart;
}

template <typename T>
Vec3<T> TF17MagneticField::getHaloField(const double &r, const double &z, const double &phi, const double &sinPhi,
                                        const double &cosPhi, const TF17Parameters<T> &p) const {
    int m;
    Vec3<T> B_cart{{0., 0., 0.}};
    auto r1_halo_r = (1. + p.a_halo * p.z1_halo * p.z1_halo) / (1. + p.a_halo * z * z);
    // cylindrical components
    T B_z0;

    auto psd = p.phi_star_disk * units::deg;
    auto psh = p.phi_star_halo * units::deg;
    auto p_0 = p.p_0 * units::deg;
    auto cot_p0 = cos(p_0) / sin(p_0);

    if (active_halo_model == "C0") { // m = 0
        B_z0 = p.B1_halo * exp(-r1_halo_r * r / p.L_halo);
    } else if (active_halo_model == "C1") { // m = 1
        // cosine argument
        auto phi_prime = phi - shiftedWindingFunction<T>(r, z, cot_p0, p) - psh;
        B_z0 = p.B1_halo * exp(-r1_halo_r * r / p.L_halo) * cos(phi_prime);
    }

    // B_r rewritten, equivalent form
    auto B_r = 2 * p.a_halo * r1_halo_r * r1_halo_r * r * z / (1. + p.a_halo * z * z) * B_z0;
    auto B_z = r1_halo_r * r1_halo_r * B_z0;
    auto B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);

    // Cartesian components
    B_cart[0] = -(B_r * cosPhi - B_phi * sinPhi); // flip x-component at the end
    B_cart[1] = B_r * sinPhi + B_phi * cosPhi;
    B_cart[2] = B_z;

    return B_cart;
}

template <typename T>
T TF17MagneticField::azimuthalFieldComponent(const double &r, const double &z, const T &B_r, const T &B_z, const T &cp0,
                                             const TF17Parameters<T> &p) const {
    auto r_ = r / p.L_p;
    auto rscale = r > epsilon ? r_ * exp(-r_) / (1 - exp(-r_)) : 1 - r_ / 2. + r_ * r_ / 12.;
    auto B_phi = cp0 / zscale(z, p) * rscale * B_r;
    B_phi = B_phi - 2 * z * r / (p.H_p * p.H_p) / zscale(z, p) * shiftedWindingFunction<T>(r, z, cp0, p) * B_z;
    return B_phi;
}

template <typename T>
T TF17MagneticField::radialFieldScale(const T &B1, const T &phi_star, const T &z1, const double &phi, const double &r,
                                      const double &z, const T &cp0, const TF17Parameters<T> &p) const {
    // cosine argument
    auto phi_prime = phi - shiftedWindingFunction<T>(r, z, cp0, p) - phi_star;
    // bisymmetric (m = 1) term
    return B1 * exp(-abs(z1) / p.H_disk) * cos(phi_prime);
}

template <typename T>
T TF17MagneticField::shiftedWindingFunction(const T &r, const double &z, const T &cp0,
                                            const TF17Parameters<T> &p) const {
    return cp0 * log(1 - exp(-r / p.L_p) + epsilon) / zscale(z, p);
}

template <typename T> T TF17MagneticField::zscale(const double &z, const TF17Parameters<T> &p) const {
    return 1 + z * z / p.H_p / p.H_p;
}

void TF17MagneticField::set_model(const std::string &dtype, const std::string &htype) {
    if (std::find(available_disk_models.begin(), available_disk_models.end(), dtype) == available_disk_models.end())
        throw std::invalid_argument("Unknown TF17 disk model '" + dtype + "'.");
    if (std::find(available_halo_models.begin(), available_halo_models.end(), htype) == available_halo_models.end())
        throw std::invalid_argument("Unknown TF17 halo model '" + htype + "'.");

    // disk parameters

    bool isAd1 = (dtype == "Ad1");
    bool isBd1 = (dtype == "Bd1");
    bool isDd1 = (dtype == "Dd1");

    bool isC0 = (htype == "C0");
    bool isC1 = (htype == "C1");

    bool isAd1andC0 = (isAd1 && isC0);
    bool isAd1andC1 = (isAd1 && isC1);

    bool isBd1andC0 = (isBd1 && isC0);
    bool isBd1andC1 = (isBd1 && isC1);

    bool isDd1andC0 = (isDd1 && isC0);
    bool isDd1andC1 = (isDd1 && isC1);

    active_disk_model = dtype;
    active_halo_model = htype;

    if (isAd1andC0) {
        active_parameters = {"a_disk",  "r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo",
                             "z1_halo", "B1_halo", "L_halo",  "p_0",           "H_p",    "L_p"};
    }

    if (isAd1andC1) {
        active_parameters = {"a_disk",  "r1_disk", "B1_disk",       "phi_star_disk", "H_disk", "a_halo", "z1_halo",
                             "B1_halo", "L_halo",  "phi_star_halo", "p_0",           "H_p",    "L_p"};
    }

    if (isBd1andC0) {
        active_parameters = {"r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo",
                             "B1_halo", "L_halo",  "p_0",           "H_p",    "L_p"};
    }

    if (isBd1andC1) {
        active_parameters = {"r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo",
                             "B1_halo", "L_halo",  "phi_star_halo", "p_0",    "H_p",    "L_p"};
    }

    if (isDd1andC0) {
        active_parameters = {"z1_disk", "B1_disk", "L_disk", "phi_star_disk", "a_halo", "z1_halo",
                             "B1_halo", "L_halo",  "p_0",    "H_p",           "L_p"};
    }

    if (isDd1andC1) {
        active_parameters = {"z1_disk", "B1_disk", "L_disk",        "phi_star_disk", "a_halo", "z1_halo",
                             "B1_halo", "L_halo",  "phi_star_halo", "p_0",           "H_p",    "L_p"};
    }
    // disk parameters
    parameters.a_disk = isAd1andC0 ? 0.9 : isAd1andC1 ? 0.031 : 0.;
    parameters.r1_disk = isDd1 ? 0. : 3.;
    parameters.z1_disk = isDd1 ? 1.5 : 0.;
    parameters.B1_disk = isAd1andC0   ? 19.
                         : isAd1andC1 ? 32.
                         : isBd1andC0 ? 2.
                         : isBd1andC1 ? 24.
                         : isDd1andC0 ? 0.065
                                      : 0.40;
    parameters.H_disk = isAd1andC0 ? 0.055 : isAd1andC1 ? 0.054 : isBd1andC0 ? .32 : isBd1andC1 ? .09 : 0;

    parameters.phi_star_disk = isAd1andC0   ? -54.
                               : isAd1andC1 ? -31.
                               : isBd1andC0 ? 153.
                               : isBd1andC1 ? -34.
                               : isDd1andC0 ? 14.
                                            : 120.;
    parameters.L_disk = isDd1andC0 ? 9.8 : isDd1andC1 ? 2.9 : 0.;

    // halo parameters

    parameters.z1_halo = 0;
    parameters.B1_halo = isAd1andC0   ? .36
                         : isAd1andC1 ? 9.
                         : isBd1andC0 ? .29
                         : isBd1andC1 ? 8.2
                         : isDd1andC0 ? 0.18
                                      : 9.5;
    parameters.L_halo = isAd1andC0   ? 3.0
                        : isAd1andC1 ? 2.1
                        : isBd1andC0 ? 3.4
                        : isBd1andC1 ? 2.2
                        : isDd1andC0 ? 4.8
                                     : 2.1;
    parameters.a_halo = isAd1andC0   ? 1.17
                        : isAd1andC1 ? 0.33
                        : isBd1andC0 ? 0.88
                        : isBd1andC1 ? 0.38
                        : isDd1andC0 ? 0.61
                                     : 0.45;
    parameters.phi_star_halo = isAd1andC1 ? 198 : isBd1andC1 ? 197 : isDd1andC1 ? 179 : 0;

    // shared parameters
    parameters.p_0 = isAd1andC0   ? -7.9
                     : isAd1andC1 ? -9.1
                     : isBd1andC0 ? -7.2
                     : isBd1andC1 ? -9.0
                     : isDd1andC0 ? -7.4
                                  : -8.4;
    parameters.H_p = isAd1andC0 ? 5. : isBd1andC0 ? 9. : isDd1andC0 ? 4.2 : 1.2;

    parameters.L_p = 50;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(TF17MagneticField)

}
