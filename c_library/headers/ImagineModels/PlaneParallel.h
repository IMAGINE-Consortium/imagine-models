// Reference: Ocker, Cordes & Chatterjee 2020, doi:10.3847/1538-4357/ab98f9, arXiv:2004.11921 (eq. 1, Ocker20); Gaensler, Madsen, Chatterjee & Mao 2008, doi:10.1071/AS08004, arXiv:0808.2550 (eq. 5, Gaensler08)
// Deviations:
// - Ocker20: smooth plane-parallel component only; the paper's clumps and voids are specific to single lines of sight

#pragma once

#include <array>
#include <string>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define PLANE_PARALLEL_PARAMETERS(X) \
    X(n0, 0.015)                     \
    X(z0, 1.57)

IMAGINE_PARAMETERS(PlaneParallelParameters, PLANE_PARALLEL_PARAMETERS)

class PlaneParallelDensity : public RegularScalarModel<PlaneParallelDensity, PlaneParallelParameters> {
public:
    const std::array<std::string, 2> available_models{"Ocker20", "Gaensler08"};
    explicit PlaneParallelDensity(const std::string &model = "Ocker20") { set_model(model); }
    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

    template <typename T>
    T field(const double &x, const double &y, const double &z, const PlaneParallelParameters<T> &p) const;

private:
    std::string active_model = "Ocker20";
};

}
