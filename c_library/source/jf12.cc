#include "ImagineModels/JF12.h"
#include "ImagineModels/units.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace imagine {

void JF12MagneticField::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown JF12 model '" + model + "'.");
    active_model = model;
    parameters = JF12Parameters<double>{};
    arm_shift = 1.;
    if (model == "JF12")
        return;
    parameters.b_arm_6 = -3.5;
    parameters.B0_X = 1.8;
    if (model == "Planck12c") {
        parameters.Bn = 1.;
        parameters.Bs = -0.8;
        parameters.B0_X = 3.;
        parameters.b_arm_2 = 2.;
        parameters.b_arm_4 = 2.;
        parameters.b_arm_5 = -3.;
        arm_shift = 0.97;
    }
}

template <typename T>
std::array<T, 2> JF12MagneticField::solenoidal_disk(double r, double phi, const JF12Parameters<T> &p) const {
    std::array<T, 2> b{{0., 0.}};
    const double r1 = rmin, r2 = Rmax;
    if (!(r1 < r && r < r2))
        return b;
    const double pitch = inc * units::deg;
    const double sin_p = std::sin(pitch), cos_p = std::cos(pitch), cot_p = 1. / std::tan(pitch);

    // arm sectors at r1
    double arms[11];
    for (int i = 1; i < 9; ++i)
        arms[i] = units::pi - cot_p * std::log(rc_B[i - 1] * arm_shift / r1);
    arms[0] = arms[8] + units::twopi;
    arms[9] = arms[1] - units::twopi;
    arms[10] = arms[2] - units::twopi;
    T bd[11] = {0., p.b_arm_1, p.b_arm_2, p.b_arm_3, p.b_arm_4, p.b_arm_5, p.b_arm_6, p.b_arm_7, 0., 0., 0.};
    T flux = 0.;
    for (int i = 1; i < 8; ++i)
        flux += (arms[i - 1] - arms[i]) * bd[i];
    bd[8] = -flux / (arms[7] - arms[8]);
    bd[0] = bd[8];
    bd[9] = bd[1];
    bd[10] = bd[2];

    // azimuthal flux integral
    T coeff[10];
    coeff[0] = 0.;
    for (int i = 1; i < 10; ++i)
        coeff[i] = coeff[i - 1] + (bd[i - 1] - bd[i]) * arms[i - 1];
    const double phi0 = 0.;
    int idx0 = 1;
    while (phi0 < arms[idx0])
        ++idx0;
    const T corr = coeff[idx0] + bd[idx0] * phi0;
    for (int i = 1; i < 10; ++i)
        coeff[i] -= corr;
    double phi1 = phi - std::log(r / r1) * cot_p;
    phi1 = std::atan2(std::sin(phi1), std::cos(phi1));
    int idx = 1;
    while (phi1 < arms[idx])
        ++idx;
    const T h = phi1 * bd[idx] + coeff[idx];

    // transition polynomial
    const double r1s = r1 + solenoidal_delta;
    const double r2s = solenoidal_outer ? r2 - solenoidal_delta : r2;
    double pd = r1 / r, qd = 0.;
    if (!(r > r1s && r < r2s)) {
        const double ra = r >= r2s ? r2 : r1, rb = r >= r2s ? r2s : r1s;
        const double fakt = (ra / rb - 2.) / ((ra - rb) * (ra - rb));
        pd = (r1 / rb) * (2. - r / rb + fakt * (r - rb) * (r - rb));
        qd = (r1 / rb) * (2. - 2. * r / rb + fakt * (3. * r * r - 4. * r * rb + rb * rb));
    }
    b[0] = pd * bd[idx] * sin_p;
    b[1] = pd * bd[idx] * cos_p - qd * h * sin_p;
    return b;
}

template <typename T>
std::array<T, 2> JF12MagneticField::solenoidal_x(double r, double z, const JF12Parameters<T> &p) const {
    const double zs = solenoidal_zs;
    const T theta0 = p.Xtheta_const * units::deg;
    const T sin0 = sin(theta0), cos0 = cos(theta0), tan0 = tan(theta0);
    const T rxc = p.rpc_X, rx = p.r0_X, bx = p.B0_X;
    const double az = std::abs(z);
    // straight field lines
    if (az > zs || zs <= 0.) {
        if (r == 0.) {
            const T q = 1. + az / tan0 / rxc;
            return {T(0.), T(bx / (q * q))};
        }
        const T rc = rxc + az / tan0;
        T mag, s, c;
        if (r < rc) {
            const T rp = r * rxc / rc;
            mag = bx * exp(-rp / rx) * (rxc / rc) * (rxc / rc);
            const T theta = z == 0. ? T(units::halfpi) : T(atan(az / (r - rp)));
            s = sin(theta);
            c = cos(theta);
        } else {
            const T rp = r - az / tan0;
            mag = bx * exp(-rp / rx) * (rp / r);
            s = sin0;
            c = cos0;
        }
        const double zsign = z < 0. ? -1. : 1.;
        return {T(zsign * mag * c), T(mag * s)};
    }
    // parabolic field lines
    const T r0c = rxc + zs / tan0;
    const double shape = zs - z * z / zs;
    T r0 = r / (1. - 1. / (2. * (zs + rxc * tan0)) * shape);
    T f;
    if (r0 >= r0c) {
        r0 = r + 1. / (2. * tan0) * shape;
        f = 1. + 1. / (2. * r * tan0 / zs) * (1. - (z / zs) * (z / zs));
    } else {
        const T g = 1. - 1. / (2. + 2. * (rxc * tan0 / zs)) * (1. - (z / zs) * (z / zs));
        f = 1. / (g * g);
    }
    T br0, bz0;
    if (r0 < r0c) {
        const T rp = r0 * rxc / r0c;
        const T theta = value(r0 - rp) == 0. ? T(units::halfpi) : T(atan(zs / (r0 - rp)));
        const T mag = bx * exp(-rp / rx) * (rxc / r0c) * (rxc / r0c);
        br0 = mag * cos(theta);
        bz0 = mag * sin(theta);
    } else {
        const T rp = r0 - zs / tan0;
        const T mag = bx * exp(-rp / rx) * (rp / r0);
        br0 = mag * cos0;
        bz0 = mag * sin0;
    }
    return {T(z / zs * f * br0), T(bz0 * f)};
}

template <typename T>
Vec3<T> JF12MagneticField::field(const double &x, const double &y, const double &z, const JF12Parameters<T> &p) const {
    const double r{sqrt(x * x + y * y)};
    const double rho{sqrt(x * x + y * y + z * z)};
    const double phi{atan2(y, x)};

    // zero outside the Galaxy
    if (r > Rmax || rho < rho_GC) {
        return Vec3<T>{{0., 0., 0.}};
    }

    // disk
    const double B0 = (rmin / r); // 1 at r = 5 kpc
    // logistic disk-halo transition
    const auto zprofile{1. / (1 + exp(-2. / p.w_disk * (std::abs(z) - p.h_disk)))};

    T B_cyl[3] = {0, 0, 0}; // cylindrical components

    if ((r > rcent)) // disk field zero elsewhere
    {
        if (r < rmin) { // circular field in molecular ring
            B_cyl[1] = B0 * p.b_ring * (1 - zprofile);
        } else if (solenoidal) {
            const auto b = solenoidal_disk(r, phi, p);
            B_cyl[0] = b[0] * (1 - zprofile);
            B_cyl[1] = b[1] * (1 - zprofile);
        } else {
            // b8 from flux conservation
            T bv_B[8] = {p.b_arm_1, p.b_arm_2, p.b_arm_3, p.b_arm_4, p.b_arm_5, p.b_arm_6, p.b_arm_7, 0.};
            T b8 = 0.;

            for (int i = 0; i < 7; i++) {
                b8 -= f[i] * bv_B[i] / f[7];
            }
            bv_B[7] = b8;

            // find spiral region
            T b_disk = 0.;
            double r_negx = r * exp(-1 / tan(units::deg * (90 - inc)) * (phi - units::pi));

            if (r_negx > rc_B[7] * arm_shift) {
                r_negx = r * exp(-1 / tan(units::deg * (90 - inc)) * (phi + units::pi));
            }
            if (r_negx > rc_B[7] * arm_shift) {
                r_negx = r * exp(-1 / tan(units::deg * (90 - inc)) * (phi + 3 * units::pi));
            }
            for (int i = 7; i >= 0; i--) {
                if (r_negx < rc_B[i] * arm_shift) {
                    b_disk = bv_B[i];
                }
            } // "region 8,7,6,..,2"

            B_cyl[0] = b_disk * B0 * sin(units::deg * inc) * (1 - zprofile);
            B_cyl[1] = b_disk * B0 * cos(units::deg * inc) * (1 - zprofile);
        }
    }

    // toroidal halo

    if (do_halo) {
        T b1, rh;
        T B_h = 0.;

        if (z >= 0) { // North
            b1 = p.Bn;
            rh = p.rn; // transition radius
        } else {       // South
            b1 = p.Bs;
            rh = p.rs;
        }

        B_h = b1 * (1. - 1. / (1. + exp(-2. / p.wh * (r - rh)))) *
              exp(-(std::abs(z)) / (p.z0)); // vertical exponential fall-off
        const T B_cyl_h[3] = {0., B_h * zprofile, 0.};
        // add fields together
        B_cyl[0] += B_cyl_h[0];
        B_cyl[1] += B_cyl_h[1];
        B_cyl[2] += B_cyl_h[2];
    }

    // X-field

    if (do_X && solenoidal) {
        const auto b = solenoidal_x(r, z, p);
        B_cyl[0] += b[0];
        B_cyl[2] += b[1];
    } else if (do_X) {
        T Xtheta = 0.;
        T rp_X = 0.; // mid-plane radius of field line
        T B_X = 0.;
        double r_sign = 1.; // +1 for north, -1 for south
        if (z < 0) {
            r_sign = -1.;
        }

        // interior/exterior boundary
        T rc_X = p.rpc_X + std::abs(z) / tan(p.Xtheta_const * units::deg);
        if (r < rc_X) { // interior, varying elevation
            rp_X = r * p.rpc_X / rc_X;
            B_X = p.B0_X * pow(p.rpc_X / rc_X, 2.) * exp(-rp_X / p.r0_X);
            Xtheta = atan(std::abs(z) / (r - rp_X)); // interior elevation angle
            if (z == 0. or r == 0.) {
                Xtheta = units::pi / 2.;
            } // to avoid some NaN
        } else { // exterior, constant elevation
            Xtheta = p.Xtheta_const * units::deg;
            rp_X = r - std::abs(z) / tan(Xtheta);
            B_X = p.B0_X * rp_X / r * exp(-rp_X / p.r0_X);
        }

        // cylindrical components
        T B_cyl_X[3] = {B_X * cos(Xtheta) * r_sign, 0., B_X * sin(Xtheta)};
        // add fields together
        B_cyl[0] += B_cyl_X[0];
        B_cyl[1] += B_cyl_X[1];
        B_cyl[2] += B_cyl_X[2];
    }

    // convert field to cartesian coordinates
    Vec3<T> B_cart{{0.0, 0.0, 0.0}};
    B_cart[0] = B_cyl[0] * cos(phi) - B_cyl[1] * sin(phi);
    B_cart[1] = B_cyl[0] * sin(phi) + B_cyl[1] * cos(phi);
    B_cart[2] = B_cyl[2];
    return B_cart;
}

IMAGINE_INSTANTIATE_VECTOR_MODEL(JF12MagneticField)

}
