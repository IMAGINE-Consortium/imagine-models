#include "TF17.h"
#include "helpers.h"

// Terral, Ferriere 2017 - Constraints from Faraday rotation on the magnetic field structure in the galactic halo, DOI: 10.1051/0004-6361/201629572, arXiv:1611.10222, implementation adapted from CRPRopa

vector TFMagneticField::_at_position(const double &x, const double &y, const double &z, const TFMagneticField &p) const
{
    //vector B_cart{{0., 0., 0.}};
    const double r = sqrt(x * x + y * y);
    double phi = M_PI - std::atan2(y, x);

    // double cosPhi = pos.x / r;
    double cosPhi = cos(phi);
    // double sinPhi = pos.y / r;
    double sinPhi = sin(phi);

    vector df = getDiskField(r, z, phi, sinPhi, cosPhi, p);
    vector hf = getHaloField(r, z, phi, sinPhi, cosPhi, p);

    vector B_cart = addVector<vector>({df, hf});
    return B_cart;
}

#if autodiff_FOUND

Eigen::MatrixXd TFMagneticField::_jac(const double &x, const double &y, const double &z, TFMagneticField &p) const
{
    vector out;
    Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, TFMagneticField &_p)
                                          { return _p._at_position(_x, _y, _z, _p); },
                                          ad::wrt(p.a_disk, p.z1_disk, p.r1_disk, p.B1_disk, p.L_disk, p.phi_star_disk, p.H_disk, p.a_halo, p.z1_halo, p.B1_halo, p.L_halo, p.phi_star_halo, p.p_0, p.H_p, p.L_p), ad::at(x, y, z, p), out);
    return _filter_diff(_deriv);
}

#endif

vector TFMagneticField::getDiskField(const double &r, const double &z, const double &phi, const double &sinPhi, const double &cosPhi, const TFMagneticField &p) const
{
    vector B_cart{{0., 0., 0.}};
    number B_r = 0;
    number B_phi = 0;
    number B_z = 0;

    auto psd = p.phi_star_disk * M_PI / 180;
    auto psh = p.phi_star_halo * M_PI / 180;
    auto p_0 = p.p_0 * M_PI / 180;
    auto cot_p0 = cos(p_0) / sin(p_0);

    if (activeDiskModel == "Ad1")
    { // ==========================================================
        if (r > r1_disk)
        {
            auto z1_disk_z = (1. + a_disk * r1_disk * r1_disk) / (1. + a_disk * r * r); // z1_disk / z, eq. 4
            // B components in (r, phi, z)
            auto B_r0 = radialFieldScale(B1_disk, psd, z1_disk_z * z, phi, r, z, cot_p0, p);
            B_r = (r1_disk / r) * z1_disk_z * B_r0;
            B_z = 2 * a_disk * r1_disk * z1_disk_z * z / (1 + a_disk * r * r) * B_r0;
            B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);
        }
        else
        {
            // within r = r1_disk, the field lines are straight in direction g_phi + phi_star_disk
            // and thus z = z1
            auto phi1_disk = shiftedWindingFunction(p.r1_disk, z, cot_p0, p) + psd;
            auto B_amp = B1_disk * exp(-fabs(z) / H_disk);
            B_r = cos(phi1_disk - phi) * B_amp;
            B_phi = sin(phi1_disk - phi) * B_amp;
        }
    }
    else if (activeDiskModel == "Bd1")
    { // ===================================================
        // for model Bd1, best fit for n = 2
        if (r > epsilon)
        {
            auto r1_disk_r = p.r1_disk / r;
            auto z1_disk_z = 5. / (r1_disk_r * r1_disk_r + 4. / sqrt(r1_disk_r)); // z1_disk / z -> remove z dependancy
            auto B_r0 = radialFieldScale(p.B1_disk, psd, z1_disk_z * z, phi, r, z, cot_p0, p);
            B_r = r1_disk_r * z1_disk_z * B_r0;
            B_z = -0.4 * r1_disk_r / r * z1_disk_z * z1_disk_z * z * (r1_disk_r * r1_disk_r - 1. / sqrt(r1_disk_r)) * B_r0;
        }
        else
        {
            auto z1_disk_z = 5. * r * r / (p.r1_disk * p.r1_disk); // z1_disk / z -> remove z dependancy
            auto B_r0 = radialFieldScale(p.B1_disk, psd, z1_disk_z * z, phi, r, z, cot_p0, p);
            B_r = 5. * r / p.r1_disk * B_r0;
            B_z = -10. * z / p.r1_disk * B_r0;
        }
        B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);
    }
    else if (activeDiskModel == "Dd1")
    { // ===================================================
        // for model Dd1, best fit for n = 0.5
        double z_sign = z >= 0 ? 1. : -1.;
        double z_abs = fabs(z);
        if (z_abs > p.epsilon)
        {
            auto z1_disk_z = p.z1_disk / z_abs;
            auto r1_disk_r = 1.5 / (sqrt(z1_disk_z) + 0.5 / z1_disk_z); // r1_disk / r
            auto F_r = r1_disk_r * r <= p.L_disk ? 1. : exp(1. - r1_disk_r * r / p.L_disk);
            // simplication of the equation in the cosinus
            auto B_z0 = z_sign * p.B1_disk * F_r * cos(phi - shiftedWindingFunction(r, z, cot_p0, p) - psd);
            B_r = -0.5 / 1.5 * r1_disk_r * r1_disk_r * r1_disk_r * r / z_abs * (sqrt(z1_disk_z) - 1 / z1_disk_z) * B_z0;
            B_z = z_sign * r1_disk_r * r1_disk_r * B_z0;
        }
        else
        {
            auto z_z1_disk = z_abs / p.z1_disk;
            auto r1_disk_r = 1.5 * sqrt(z_abs / p.z1_disk); // r1_disk / r
            auto F_r = r1_disk_r * r <= p.L_disk ? 1. : exp(1. - r1_disk_r * r / p.L_disk);
            auto B_z0 = z_sign * p.B1_disk * F_r * cos(phi - shiftedWindingFunction(r, z, cot_p0, p) - psd);
            B_r = -1.125 * r / p.z1_disk * (1 - 2.5 * z_z1_disk * sqrt(z_z1_disk)) * B_z0;
            B_z = z_sign * r1_disk_r * r1_disk_r * B_z0;
        }
        B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);
    }

    // Convert to (x, y, z) components
    B_cart[0] = -(B_r * cosPhi - B_phi * sinPhi); // flip x-component at the end
    B_cart[1] = B_r * sinPhi + B_phi * cosPhi;
    B_cart[2] = B_z;
    return B_cart;
}

vector TFMagneticField::getHaloField(const double &r, const double &z, const double &phi, const double &sinPhi, const double &cosPhi, const TFMagneticField &p) const
{
    int m;
    vector B_cart{{0., 0., 0.}};
    auto r1_halo_r = (1. + p.a_halo * p.z1_halo * p.z1_halo) / (1. + p.a_halo * z * z);
    // B components in (r, phi, z)
    number B_z0;

    auto psd = p.phi_star_disk * M_PI / 180;
    auto psh = p.phi_star_halo * M_PI / 180;
    auto p_0 = p.p_0 * M_PI / 180;
    auto cot_p0 = cos(p_0) / sin(p_0);

    if (activeHaloModel == "C0")
    { // m = 0
        B_z0 = B1_halo * exp(-r1_halo_r * r / L_halo);
    }
    else if (activeHaloModel == "C1")
    { // m = 1
        // simplication of the equation in the cosinus
        auto phi_prime = phi - shiftedWindingFunction(r, z, cot_p0, p) - psd;
        B_z0 = B1_halo * exp(-r1_halo_r * r / L_halo) * cos(phi_prime);
    }

    // Contrary to article, Br has been rewriten to a little bit by replacing
    // (2 * a * r1**3 * z) / (r**2) by (2 * a * r1**2 * z) / (r * (1+a*z**2))
    // but that is strictly equivalent except we can reintroduce the z1 in the expression via r1
    auto B_r = 2 * p.a_halo * r1_halo_r * r1_halo_r * r * z / (1. + a_halo * z * z) * B_z0;
    auto B_z = r1_halo_r * r1_halo_r * B_z0;
    auto B_phi = azimuthalFieldComponent(r, z, B_r, B_z, cot_p0, p);

    // Convert to (x, y, z) components
    B_cart[0] = -(B_r * cosPhi - B_phi * sinPhi); // flip x-component at the end
    B_cart[1] = B_r * sinPhi + B_phi * cosPhi;
    B_cart[2] = B_z;

    return B_cart;
}

number TFMagneticField::azimuthalFieldComponent(const double &r, const double &z, const number &B_r, const number &B_z, const number &cp0, const TFMagneticField &p) const
{
    auto r_ = r / p.L_p;
    auto rscale = r > p.epsilon ? r_ * exp(-r_) / (1 - exp(-r_)) : 1 - r_ / 2. - r_ * r_ / 12.;
    auto B_phi = cp0 / zscale(z, p) * rscale * B_r;
    B_phi = B_phi - 2 * z * r / (p.H_p * p.H_p) / zscale(z, p) * shiftedWindingFunction(r, z, cp0, p) * B_z;
    return B_phi;
}

number TFMagneticField::radialFieldScale(const number &B1, const number &phi_star, const number &z1, const double &phi, const double &r, const double &z, const number &cp0, const TFMagneticField &p) const
{
    // simplication of the equation in the cosinus
    auto phi_prime = phi - shiftedWindingFunction(r, z, cp0, p) - phi_star;
    // This term occures is parameterizations of models A and B always bisymmetric (m = 1)
    return B1 * exp(-abs(z1) / p.H_disk) * cos(phi_prime);
}

number TFMagneticField::shiftedWindingFunction(const number &r, const double &z, const number &cp0, const TFMagneticField &p) const
{
    return cp0 * log(1 - exp(-r / p.L_p) + p.epsilon) / zscale(z, p);
}

number TFMagneticField::zscale(const double &z, const TFMagneticField &p) const
{
    return 1 + z * z / p.H_p / p.H_p;
}

void TFMagneticField::set_params(std::string dtype, std::string htype)
{
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

    activeDiskModel = dtype;
    activeHaloModel = htype;

#if autodiff_FOUND
    if (isAd1andC0)
    {
        active_diff = {"a_disk", "r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "p_0", "H_p", "L_p"};
    }

    if (isAd1andC1)
    {
        active_diff = {"a_disk", "r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "phi_star_halo", "p_0", "H_p", "L_p"};
    }

    if (isBd1andC0)
    {
        active_diff = {"r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "p_0", "H_p", "L_p"};
    }

    if (isBd1andC1)
    {
        active_diff = {"r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "phi_star_halo", "p_0", "H_p", "L_p"};
    }

    if (isDd1andC0)
    {
        active_diff = {"z1_disk", "B1_disk", "L_disk", "phi_star_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "p_0", "H_p", "L_p"};
    }

    if (isDd1andC1)
    {
        active_diff = {"z1_disk", "B1_disk", "L_disk", "phi_star_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "phi_star_halo", "p_0", "H_p", "L_p"};
    }
#endif
    // disk parameters
    a_disk = isAd1andC0 ? 0.9 : isAd1andC1 ? 0.031
                                           : 0.;
    r1_disk = isDd1 ? 0. : 3.;
    z1_disk = isDd1 ? 1.5 : 0.;
    B1_disk = isAd1andC0 ? 19. : isAd1andC1 ? 32.
                             : isBd1andC0   ? 2.
                             : isBd1andC1   ? 24.
                             : isDd1andC0   ? 0.065
                                            : 0.40;
    H_disk = isAd1andC0 ? 0.055 : isAd1andC1 ? 0.054
                              : isBd1andC0   ? .32
                              : isBd1andC1   ? .09
                                             : 0;

    phi_star_disk = isAd1andC0 ? -54. : isAd1andC1 ? -31.
                                    : isBd1andC0   ? 153.
                                    : isBd1andC1   ? -34.
                                    : isDd1andC0   ? 14.
                                                   : 120.;
    L_disk = isDd1andC0 ? 9.8 : isDd1andC1 ? 2.9
                                           : 0.;

    // halo parameters

    z1_halo = 0;
    B1_halo = isAd1andC0 ? .36 : isAd1andC1 ? 9.
                             : isBd1andC0   ? .29
                             : isBd1andC1   ? 8.2
                             : isDd1andC0   ? 0.18
                                            : 9.5;
    L_halo = isAd1andC0 ? 3.0 : isAd1andC1 ? 2.1
                            : isBd1andC0   ? 3.4
                            : isBd1andC1   ? 2.2
                            : isDd1andC0   ? 4.8
                                           : 2.1;
    a_halo = isAd1andC0 ? 1.17 : isAd1andC1 ? 0.33
                             : isBd1andC0   ? 0.88
                             : isBd1andC1   ? 0.38
                             : isDd1andC0   ? 0.61
                                            : 0.45;
    phi_star_halo = isAd1andC1 ? 198 : isBd1andC1 ? 197
                                   : isDd1andC1   ? 179
                                                  : 0;

    // shared parameters
    p_0 = isAd1andC0 ? -7.9 : isAd1andC1 ? -9.1
                          : isBd1andC0   ? -7.2
                          : isBd1andC1   ? -9.0
                          : isDd1andC0   ? -7.4
                                         : -8.4;
    H_p = isAd1andC0 ? 5. : isBd1andC0 ? 9.
                        : isDd1andC0   ? 4.2
                                       : 1.2;

    L_p = 50;
}
