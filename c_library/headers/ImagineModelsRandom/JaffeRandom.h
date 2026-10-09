// Reference: Jaffe et al. 2013, arXiv:1302.0143 (Table A1)
// Based on: conventions of hammurabi v3.01 (field_bb, bran_arms_compress), GPL-3.0
// Deviations:
// - arm and ring geometry, field directions and arm widths of JaffeMagneticField("Jaffe13") (see its deviations)
// - ring included like an arm, weighted by |ring_amp| (as in hammurabi v3.01; Table A1 sums over arms)
// - field set to zero beyond 20 kpc (sphere), as the coherent field

#pragma once

#include <array>
#include <string>

#include "ImagineModels/Jaffe.h"
#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

class JaffeRandomField : public RandomVectorField {
public:
    const std::array<std::string, 1> available_models{"Jaffe13"};

    explicit JaffeRandomField(const std::string &model = "Jaffe13") { set_model(model); }

    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    double b_rms = 3.5; // muG
    double h_rms = 2.;  // kpc
    double r_grf = 20.; // kpc
    double f_ord = 0.15;
    double spectral_offset = 0.;
    double spectral_slope = 0.37;

    JaffeMagneticField regular_base = JaffeMagneticField("Jaffe13");

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;
    double isotropic_rms(const double &x, const double &y, const double &z) const override;
    double ordered_amplitude(const double &x, const double &y, const double &z) const override;
    Vec3<double> anisotropy_direction(const double &x, const double &y, const double &z) const override;

private:
    std::string active_model = "Jaffe13";
    bool outside(const double &x, const double &y, const double &z) const;
    double compression(const double &x, const double &y, const double &z) const;
};

}
