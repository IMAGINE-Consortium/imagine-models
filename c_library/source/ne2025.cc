/*
This file contains code adapted from mwprop (NE2025 Fortran code, density.NE2025.f, neLISM.NE2025.f,
neclumpN.NE2025.f, nevoidN.NE2025.f).

The original copyright statement is reproduced below:

Copyright (C) 2026, J.M. Cordes & S.K. Ocker

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModels/NE2025.h"

namespace imagine {

namespace {

constexpr double rad = 57.2957795130823;
constexpr double radian = 57.29577951;
constexpr double pihalf = 1.5707963267948966;
constexpr int narms = 5;
constexpr int armmap[narms] = {1, 3, 4, 2, 5};
constexpr double arm_a[narms] = {4.25, 4.25, 4.89, 4.89, 4.57};
constexpr double arm_rmin[narms] = {3.48, 3.48, 4.90, 3.76, 8.10};
constexpr double arm_thmin[narms] = {0., 3.141, 2.525, 4.24, 5.847};
constexpr double arm_extent[narms] = {6., 6., 6., 6., 0.55};

template <typename T> T sech2(const T &z) {
    if (abs(z) >= 20.)
        return T(0.);
    const T s = 2. / (exp(z) + exp(-z));
    return s * s;
}

struct Spline {
    std::vector<double> x, y, y2;

    Spline(const std::vector<double> &x, const std::vector<double> &y) : x(x), y(y), y2(x.size(), 0.) {
        const std::size_t n = x.size();
        std::vector<double> u(n, 0.);
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
            const double p = sig * y2[i - 1] + 2.;
            y2[i] = (sig - 1.) / p;
            u[i] = (6. * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
                        (x[i + 1] - x[i - 1]) -
                    sig * u[i - 1]) /
                   p;
        }
        y2[n - 1] = 0.;
        for (std::size_t k = n - 1; k-- > 0;)
            y2[k] = y2[k] * y2[k + 1] + u[k];
    }

    double operator()(double xout) const {
        std::size_t lo = 0, hi = x.size() - 1;
        while (hi - lo > 1) {
            const std::size_t k = (hi + lo) / 2;
            if (x[k] > xout)
                hi = k;
            else
                lo = k;
        }
        const double h = x[hi] - x[lo];
        const double a = (x[hi] - xout) / h;
        const double b = (xout - x[lo]) / h;
        return a * y[lo] + b * y[hi] + ((a * a * a - a) * y2[lo] + (b * b * b - b) * y2[hi]) * h * h / 6.;
    }
};

}

NE2025::NE2025(const std::string &model) {
    setup_arms();
    set_model(model);
}

void NE2025::setup_arms() {
    constexpr int nn = 20;
    constexpr int narmpoints = 500;
    for (int j = 0; j < narms; ++j) {
        std::vector<double> th1(nn), r1(nn);
        for (int n = 0; n < nn; ++n) {
            th1[n] = arm_thmin[j] + n * arm_extent[j] / (nn - 1.);
            r1[n] = arm_rmin[j] * std::exp((th1[n] - arm_thmin[j]) / arm_a[j]);
            th1[n] *= rad;
            // arm shape modifications
            if (armmap[j] == 3) {
                if (th1[n] > 370. && th1[n] <= 410.)
                    r1[n] *= 1. + 0.04 * std::cos((th1[n] - 390.) * 180. / (40. * rad));
                if (th1[n] > 315. && th1[n] <= 370.)
                    r1[n] *= 1. - 0.07 * std::cos((th1[n] - 345.) * 180. / (55. * rad));
                if (th1[n] > 180. && th1[n] <= 315.)
                    r1[n] *= 1 + 0.16 * std::cos((th1[n] - 260.) * 180. / (135. * rad));
            }
            if (armmap[j] == 2) {
                if (th1[n] > 290. && th1[n] <= 395.)
                    r1[n] *= 1. - 0.11 * std::cos((th1[n] - 350.) * 180. / (105. * rad));
            }
        }
        const Spline spline(th1, r1);
        const double dth = 5. / r1[0];
        double th = th1[0] - 0.999 * dth;
        auto &arm = arms[j];
        arm.clear();
        int k = 1;
        for (; k <= narmpoints - 1; ++k) {
            th += dth;
            if (th > th1[nn - 1])
                break;
            const double r = spline(th);
            arm.push_back({-r * std::sin(th / rad), r * std::cos(th / rad)});
        }
        // unfilled Fortran entry kmax
        arm.push_back({0., 0.});
    }
}

void NE2025::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown NE2025 model '" + model + "'.");
    active_model = model;
    parameters = NE2025Parameters<double>{};
    const bool ne2001 = model == "NE2001";
    if (ne2001) {
        parameters.n1h1 = 0.033;
        parameters.h1 = 0.97;
        parameters.A2 = 3.8;
        parameters.narm2 = 1.2;
        parameters.narm3 = 1.3;
        parameters.narm4 = 1.;
        parameters.warm4 = 0.8;
        parameters.hgc = 0.026;
    }
    clumps.clear();
    for (const auto &c : ne2001 ? clumps_ne2001 : clumps_ne2025) {
        const double rgalc = c.d * std::cos(c.b / radian);
        clumps.push_back({rgalc * std::sin(c.l / radian), r_sun - rgalc * std::cos(c.l / radian),
                          c.d * std::sin(c.b / radian), c.ne, c.r, c.edge});
    }
    voids.clear();
    for (const auto &v : ne2001 ? voids_ne2001 : voids_ne2025) {
        const double rgalc = v.d * std::cos(v.b / radian);
        const double s1 = std::sin(v.theta_y / radian), c1 = std::cos(v.theta_y / radian);
        const double s2 = std::sin(v.theta_z / radian), c2 = std::cos(v.theta_z / radian);
        voids.push_back({rgalc * std::sin(v.l / radian), r_sun - rgalc * std::cos(v.l / radian),
                         v.d * std::sin(v.b / radian), v.ne, v.a, v.b_axis, v.c, c1 * c2, s1 * s2, c2 * s1, c1 * s2, s1,
                         c1, s2, c2, v.edge});
    }
}

template <typename T> T NE2025::arms_density(double x, double y, double z, const NE2025Parameters<T> &p) const {
    const std::array<T, narms> narm{p.narm1, p.narm2, p.narm3, p.narm4, p.narm5};
    const std::array<T, narms> warm{p.warm1, p.warm2, p.warm3, p.warm4, p.warm5};
    const std::array<T, narms> harm{p.harm1, p.harm2, p.harm3, p.harm4, p.harm5};
    constexpr int ks = 3;
    T nea = 0.;
    if (!(abs(z / p.ha) < 10.))
        return nea;
    const double rr = std::sqrt(x * x + y * y);
    double thxy = std::atan2(-x, y) * rad;
    if (thxy < 0.)
        thxy += 360.;
    for (int j = 0; j < narms; ++j) {
        const int jj = armmap[j];
        const auto &arm = arms[j];
        const int kmax = int(arm.size());
        auto sq = [&](int k) {
            const double dx = x - arm[k - 1][0], dy = y - arm[k - 1][1];
            return dx * dx + dy * dy;
        };
        // nearest arm point
        double sqmin = 1e10;
        int kk = 1;
        for (int k = 1 + ks; k <= kmax - ks; k += 2 * ks + 1) {
            if (sq(k) < sqmin) {
                sqmin = sq(k);
                kk = k;
            }
        }
        const int kmi = std::max(kk - 2 * ks, 1), kma = std::min(kk + 2 * ks, kmax);
        for (int k = kmi; k <= kma; ++k) {
            if (sq(k) < sqmin) {
                sqmin = sq(k);
                kk = k;
            }
        }
        double exx, eyy;
        if (kk > 1 && kk < kmax) {
            const int kl = sq(kk - 1) < sq(kk + 1) ? kk - 1 : kk + 1;
            const auto &a = arm[kk - 1], &b = arm[kl - 1];
            const double emm = (a[1] - b[1]) / (a[0] - b[0]);
            const double ebb = a[1] - emm * a[0];
            exx = (x + emm * y - emm * ebb) / (1. + emm * emm);
            const double test = (exx - a[0]) / (b[0] - a[0]);
            if (test < 0. || test > 1.)
                exx = a[0];
            eyy = emm * exx + ebb;
        } else {
            exx = arm[kk - 1][0];
            eyy = arm[kk - 1][1];
        }
        const double smin = std::sqrt((x - exx) * (x - exx) + (y - eyy) * (y - eyy));
        if (!(smin < 3. * p.wa))
            continue;
        const T width = warm[jj - 1] * p.wa;
        T ga = exp(-(smin / width) * (smin / width));
        if (rr > p.Aa)
            ga *= sech2(T((rr - p.Aa) / 2.));
        // arm 3 and arm 2 tapers
        const double th3a = 290., th3b = 363.;
        double test3 = thxy - th3a;
        if (test3 < 0.)
            test3 += 360.;
        if (jj == 3 && 0. <= test3 && test3 < th3b - th3a) {
            const double fac = std::pow((1. + std::cos(6.2831853 * (thxy - th3a) / (th3b - th3a))) / 2., 4.);
            ga *= fac;
        }
        const double th2a = 340., th2b = 370., fac2min = 0.1;
        double test2 = thxy - th2a;
        if (test2 < 0.)
            test2 += 360.;
        if (jj == 2 && 0. <= test2 && test2 < th2b - th2a)
            ga *= (1. + fac2min + (1. - fac2min) * std::cos(6.2831853 * (thxy - th2a) / (th2b - th2a))) / 2.;
        nea += narm[jj - 1] * p.na * ga * sech2(T(z / (harm[jj - 1] * p.ha)));
    }
    return nea;
}

template <typename T> T NE2025::lism(double x, double y, double z, const NE2025Parameters<T> &p, bool &inside) const {
    auto ellipsoid = [&](const T &a, const T &b, const T &c, const T &x0, const T &y0, const T &z0, const T &theta) {
        const T s = sin(theta / radian), co = cos(theta / radian);
        const T ap = (co / a) * (co / a) + (s / b) * (s / b);
        const T bp = (s / a) * (s / a) + (co / b) * (co / b);
        const T cp = 1. / (c * c);
        const T dp = 2. * co * s * (1. / (a * a) - 1. / (b * b));
        const T q =
            (x - x0) * (x - x0) * ap + (y - y0) * (y - y0) * bp + (z - z0) * (z - z0) * cp + (x - x0) * (y - y0) * dp;
        return q <= 1.;
    };
    const bool ldr = ellipsoid(p.aldr, p.bldr, p.cldr, p.xldr, p.yldr, p.zldr, p.thetaldr);
    const bool lsb = ellipsoid(p.alsb, p.blsb, p.clsb, p.xlsb, p.ylsb, p.zlsb, p.thetalsb);

    // Local Hot Bubble
    const T yaxis = p.ylhb + tan(p.thetalhb / radian) * z;
    T aa = p.alhb;
    if (z <= 0. && z >= p.zlhb - p.clhb)
        aa = 0.001 + (p.alhb - 0.001) * (1. - (1. / (p.zlhb - p.clhb)) * z);
    const T qxy = ((x - p.xlhb) / aa) * ((x - p.xlhb) / aa) + ((y - yaxis) / p.blhb) * ((y - yaxis) / p.blhb);
    const T qz = abs(z - p.zlhb) / p.clhb;
    const bool lhb = qxy <= 1. && qz <= 1.;

    // Loop I
    bool loop = false;
    T neloop = 0.;
    if (z >= 0.) {
        const T r = sqrt((x - p.xlpI) * (x - p.xlpI) + (y - p.ylpI) * (y - p.ylpI) + (z - p.zlpI) * (z - p.zlpI));
        if (r <= p.rlpI) {
            loop = true;
            neloop = p.nelpI;
        } else if (r <= p.rlpI + p.drlpI) {
            loop = true;
            neloop = p.dnelpI;
        }
    }

    inside = ldr || lsb || lhb || loop;
    if (lhb)
        return p.nelhb;
    if (loop)
        return neloop;
    if (lsb)
        return p.nelsb;
    return ldr ? p.neldr : T(0.);
}

template <typename T>
T NE2025::field(const double &x_in, const double &y_in, const double &z, const NE2025Parameters<T> &p) const {
    // NE2001 frame
    const double x = y_in;
    const double y = -x_in;
    const double rr = std::sqrt(x * x + y * y);

    const T g1 = rr > p.A1 ? T(0.) : T(cos(pihalf * rr / p.A1) / cos(pihalf * r_sun / p.A1));
    const T ne1 = (p.n1h1 / p.h1) * g1 * sech2(T(z / p.h1));

    const T rrarg = ((rr - p.A2) / 1.8) * ((rr - p.A2) / 1.8);
    const T ne2 = rrarg < 10. ? T(p.n2 * exp(-rrarg) * sech2(T(z / p.h2))) : T(0.);

    const T nea = arms_density(x, y, z, p);

    T negc = 0.;
    const T rgc = sqrt((x - p.xgc) * (x - p.xgc) + (y - p.ygc) * (y - p.ygc));
    const T zz = abs(z - p.zgc);
    if (rgc <= p.rgc && zz <= p.hgc && (rgc / p.rgc) * (rgc / p.rgc) + (zz / p.hgc) * (zz / p.hgc) <= 1.)
        negc = p.negc0;

    bool in_lism = false;
    const T nelism = lism(x, y, z, p, in_lism);
    T ne = in_lism ? nelism : T(ne1 + ne2 + nea + negc);

    // clumps and voids
    double necN = 0.;
    for (const auto &c : clumps) {
        const double arg = ((x - c.x) * (x - c.x) + (y - c.y) * (y - c.y) + (z - c.z) * (z - c.z)) / (c.r * c.r);
        if (c.edge == 0 && arg < 5.)
            necN += c.ne * std::exp(-arg);
        if (c.edge == 1 && arg <= 1.)
            necN += c.ne;
    }
    bool in_void = false;
    double nevN = 0.;
    for (const auto &v : voids) {
        const double dx = x - v.x, dy = y - v.y, dz = z - v.z;
        const double q1 = v.cc12 * dx + v.s2 * dy + v.cs21 * dz;
        const double q2 = -v.cs12 * dx + v.c2 * dy - v.ss12 * dz;
        const double q3 = -v.s1 * dx + v.c1 * dz;
        const double q = q1 * q1 / (v.a * v.a) + q2 * q2 / (v.b * v.b) + q3 * q3 / (v.c * v.c);
        if (v.edge == 0 && q < 3.) {
            nevN = v.ne * std::exp(-q);
            in_void = true;
        }
        if (v.edge == 1 && q <= 1.) {
            nevN = v.ne;
            in_void = true;
        }
    }
    if (in_void)
        ne = T(nevN);
    return ne + necN;
}

IMAGINE_INSTANTIATE_SCALAR_MODEL(NE2025)

}
