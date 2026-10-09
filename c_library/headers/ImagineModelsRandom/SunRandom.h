// Reference: Sun et al. 2008, arXiv:0711.1572; Sun10b variant: Planck XLII 2016, arXiv:1601.00546 (Sect. 3.3.1, Table C.1)
// Deviations:
// - the papers give only the rms; the power spectrum is the library default
// - Sun10b: ordered random component (beta = 3) not modelled

#pragma once

#include <array>
#include <string>

#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

class SunRandomField : public RandomVectorField {
public:
    const std::array<std::string, 2> available_models{"Sun10", "Sun10b"};

    explicit SunRandomField(const std::string &model = "Sun10") { set_model(model); }

    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    double b_iso = 3.;  // muG
    double r_sun = 8.5; // kpc
    double r0 = 30.;    // kpc
    double h_disk = 1.; // kpc
    double h_halo = 3.; // kpc
    double f_disk = 0.5;
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;

private:
    std::string active_model = "Sun10";
};

}
