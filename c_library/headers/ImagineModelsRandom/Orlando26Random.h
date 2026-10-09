// Reference: Orlando 2026, arXiv:2608.07679 (Sect. 2.2.2, Table 1)
// Deviations:
// - ordered random halo amplitude sqrt(B_OH^2 - B_H^2): the paper fits only the total ordered field B_OH
// - ordered random part as a Gaussian projection (g2.e)e of an independent random field along the XH24 toroid
// - the paper gives only the rms; the power spectrum is the library default
// - r_sun = 8.5 kpc (not stated in the paper)

#pragma once

#include <array>
#include <string>

#include "ImagineModels/XH24.h"
#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

class Orlando26RandomField : public RandomVectorField {
public:
    const std::array<std::string, 2> available_models{"halo4kpc", "halo10kpc"};

    explicit Orlando26RandomField(const std::string &model = "halo4kpc") { set_model(model); }

    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    double b_ran = 4.9;     // muG
    double r0_ran = 30.;    // kpc
    double z0_ran = 4.;     // kpc
    double r_sun = 8.5;     // kpc
    double b_ordered = 4.3; // muG, total ordered
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    XH24MagneticField regular_base = XH24MagneticField();

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;
    double isotropic_rms(const double &x, const double &y, const double &z) const override;
    double ordered_amplitude(const double &x, const double &y, const double &z) const override;
    Vec3<double> anisotropy_direction(const double &x, const double &y, const double &z) const override;

private:
    std::string active_model = "halo4kpc";
};

}
