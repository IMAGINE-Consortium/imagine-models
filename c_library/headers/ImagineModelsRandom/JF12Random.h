// Reference: Jansson & Farrar 2012, arXiv:1210.7820; Planck variants: Planck XLII 2016, arXiv:1601.00546 (Table C.1)
// Based on: hammurabiX, GPL-3.0; compared with CRPropa (JF12Field, PlanckJF12bField), GPL-3.0
// Deviations:
// - Planck variants without the striation factor beta

#pragma once

#include <array>
#include <cmath>
#include <string>

#include "ImagineModels/JF12.h"
#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

class JF12RandomField : public RandomVectorField {
public:
    double b0_1 = 10.81;     // muG
    double b0_2 = 6.96;      // muG
    double b0_3 = 9.59;      // muG
    double b0_4 = 6.96;      // muG
    double b0_5 = 1.96;      // muG
    double b0_6 = 16.34;     // muG
    double b0_7 = 37.29;     // muG
    double b0_8 = 10.35;     // muG
    double b0_int = 7.63;    // muG
    double z0_spiral = 0.61; // kpc
    double b0_halo = 4.68;   // muG
    double r0_halo = 10.97;  // kpc
    double z0_halo = 2.84;   // kpc
    double Rmax = 20.;
    double rho_GC = 1.;
    const double rc_B[8] = {5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5}; // arm boundaries on negative x-axis
    const double inc = 11.5;                                            // inclination, in degrees
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    JF12MagneticField regular_base = JF12MagneticField();

    const std::array<std::string, 3> available_models{"JF12", "Planck12b", "Planck12c"};
    explicit JF12RandomField(const std::string &model = "JF12") { set_model(model); }
    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }
    double arm_shift = 1.;

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;
    Vec3<double> anisotropy_direction(const double &x, const double &y, const double &z) const override;

private:
    std::string active_model = "JF12";
};

}
