#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define JF12_PARAMETERS(X) \
    X(b_arm_1, 0.1)        \
    X(b_arm_2, 3.0)        \
    X(b_arm_3, -0.9)       \
    X(b_arm_4, -0.8)       \
    X(b_arm_5, -2.0)       \
    X(b_arm_6, -4.2)       \
    X(b_arm_7, 0.0)        \
    X(b_ring, 0.1)         \
    X(h_disk, 0.40)        \
    X(w_disk, 0.27)        \
    X(Bn, 1.4)             \
    X(Bs, -1.1)            \
    X(rn, 9.22)            \
    X(rs, 16.7)            \
    X(wh, 0.20)            \
    X(z0, 5.3)             \
    X(B0_X, 4.6)           \
    X(Xtheta_const, 49.)   \
    X(rpc_X, 4.8)          \
    X(r0_X, 2.9)

IMAGINE_PARAMETERS(JF12Parameters, JF12_PARAMETERS)

class JF12MagneticField : public RegularVectorModel<JF12MagneticField, JF12Parameters> {
public:
    // define fixed parameters
    const double Rmax = 20;   // outer boundary of GMF
    const double rho_GC = 1.; // interior boundary of GMF

    // fixed disk parameters
    const double inc = 11.5; // inclination, in degrees
    const double rmin = 5.;  // outer boundary of the molecular ring region
    const double rcent = 3.; // inner boundary of the molecular ring region (field is
                             // zero within this region)
    const double f[8] = {0.130, 0.165, 0.094, 0.122,
                         0.13,  0.118, 0.084, 0.156}; // fractions of circumference spanned by each
                                                      // spiral, sums to unity
    const double rc_B[8] = {5.1, 6.3,  7.1,  8.3,
                            9.8, 11.4, 12.7, 15.5}; // the radii where the spiral arm boundaries cross the
                                                    // negative x-axis

    // toroidal halo parameters
    bool do_halo = true;
    // X-field parameters
    bool do_X = true;

    const std::array<std::string, 3> available_models{"JF12", "Planck12b", "Planck12c"};
    explicit JF12MagneticField(const std::string &model = "JF12") { set_model(model); }
    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }
    double arm_shift = 1.;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const JF12Parameters<T> &p) const;

private:
    std::string active_model = "JF12";
};

}
