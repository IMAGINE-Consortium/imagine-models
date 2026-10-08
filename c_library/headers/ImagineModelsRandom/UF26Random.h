#ifndef UF26RANDOM_H
#define UF26RANDOM_H

#include <array>
#include <string>

#include "ImagineModelsRandom/RandomVectorField.h"

namespace imagine {

class UF26RandomField : public RandomVectorField {
public:
    const std::array<std::string, 2> available_models{"expDisk", "ringDisk"};

    explicit UF26RandomField(const std::string &model = "expDisk") { set_model(model); }

    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    double b_disk = 4.4;  // muG
    double z_disk = 1.0;  // kpc
    double l_r = 15.;     // kpc
    double r_c = 1.;      // kpc
    double b_ring = 0.;   // muG
    double r_ring = 4.5;  // kpc
    double w_ring = 2.0;  // kpc
    double z_ring = 1.9;  // kpc
    double r_max = 18.;   // kpc
    double w_max = 2.;    // kpc
    double r_sun = 8.178; // kpc
    double spectral_offset = 1.;
    double spectral_slope = 2.;

    double spectrum(const double &abs_k) const override;
    double rms(const double &x, const double &y, const double &z) const override;

private:
    std::string active_model = "expDisk";
};

}

#endif
