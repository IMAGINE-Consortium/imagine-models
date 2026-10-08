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

// Reference: Ocker & Cordes 2026, arXiv:2602.11838 (NE2025); Cordes & Lazio 2002, arXiv:astro-ph/0207156 (NE2001)
// Based on: authors' Fortran code in mwprop (github.com/stella-ocker/mwprop), GPL-3.0-or-later
// Deviations:
// - electron density only; the fluctuation parameters F and scattering are not included
// - double instead of single precision
// - position parameters (Galactic Centre, local ISM) in the NE2001 frame: x towards l = 90 deg, Sun at (0, 8.5, 0) kpc

#pragma once

#include <array>
#include <string>
#include <vector>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define NE2025_PARAMETERS(X)     \
    X(n1h1, 0.0275)              \
    X(h1, 1.589)                 \
    X(A1, 17.5)                  \
    X(n2, 0.08)                  \
    X(h2, 0.15)                  \
    X(A2, 4.3)                   \
    X(na, 0.028)                 \
    X(ha, 0.23)                  \
    X(wa, 0.65)                  \
    X(Aa, 10.5)                  \
    X(narm1, 0.5)                \
    X(narm2, 1.5)                \
    X(narm3, 2.7)                \
    X(narm4, 3.7)                \
    X(narm5, 0.25)               \
    X(warm1, 1.)                 \
    X(warm2, 1.5)                \
    X(warm3, 1.)                 \
    X(warm4, 0.96)               \
    X(warm5, 1.)                 \
    X(harm1, 1.)                 \
    X(harm2, 0.8)                \
    X(harm3, 1.3)                \
    X(harm4, 1.5)                \
    X(harm5, 1.)                 \
    X(xgc, -0.01)                \
    X(ygc, 0.)                   \
    X(zgc, -0.02)                \
    X(rgc, 0.145)                \
    X(hgc, 0.05)                 \
    X(negc0, 10.)                \
    X(aldr, 1.5)                 \
    X(bldr, 0.75)                \
    X(cldr, 0.5)                 \
    X(xldr, 1.36)                \
    X(yldr, 8.06)                \
    X(zldr, 0.)                  \
    X(thetaldr, -24.2) /* deg */ \
    X(neldr, 0.012)              \
    X(alsb, 1.05)                \
    X(blsb, 0.425)               \
    X(clsb, 0.325)               \
    X(xlsb, -0.75)               \
    X(ylsb, 9.)                  \
    X(zlsb, -0.05)               \
    X(thetalsb, 139.) /* deg */  \
    X(nelsb, 0.016)              \
    X(alhb, 0.085)               \
    X(blhb, 0.1)                 \
    X(clhb, 0.33)                \
    X(xlhb, 0.01)                \
    X(ylhb, 8.45)                \
    X(zlhb, 0.17)                \
    X(thetalhb, 15.) /* deg */   \
    X(nelhb, 0.005)              \
    X(xlpI, -0.045)              \
    X(ylpI, 8.4)                 \
    X(zlpI, 0.07)                \
    X(rlpI, 0.12)                \
    X(drlpI, 0.06)               \
    X(nelpI, 0.0125)             \
    X(dnelpI, 0.0125)

IMAGINE_PARAMETERS(NE2025Parameters, NE2025_PARAMETERS)

class NE2025 : public RegularScalarModel<NE2025, NE2025Parameters> {
public:
    struct Clump {
        double l, b, ne, d, r;
        int edge;
    };
    struct Void {
        double l, b, d, ne, a, b_axis, c, theta_y, theta_z;
        int edge;
    };

    const std::array<std::string, 2> available_models{"NE2025", "NE2001"};
    explicit NE2025(const std::string &model = "NE2025");
    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    template <typename T>
    T field(const double &x, const double &y, const double &z, const NE2025Parameters<T> &p) const;

private:
    static constexpr double r_sun = 8.5;
    static const std::vector<Clump> clumps_ne2025, clumps_ne2001;
    static const std::vector<Void> voids_ne2025, voids_ne2001;

    struct ClumpXYZ {
        double x, y, z, ne, r;
        int edge;
    };
    struct VoidXYZ {
        double x, y, z, ne, a, b, c, cc12, ss12, cs21, cs12, s1, c1, s2, c2;
        int edge;
    };

    std::string active_model = "NE2025";
    std::vector<ClumpXYZ> clumps;
    std::vector<VoidXYZ> voids;
    std::array<std::vector<std::array<double, 2>>, 5> arms;

    void setup_arms();
    template <typename T> T arms_density(double x, double y, double z, const NE2025Parameters<T> &p) const;
    template <typename T> T lism(double x, double y, double z, const NE2025Parameters<T> &p, bool &inside) const;
};

}
