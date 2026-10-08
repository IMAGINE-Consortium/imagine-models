#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "ImagineModels/PlaneParallel.h"

namespace imagine {

void PlaneParallelDensity::set_model(const std::string &model) {
    if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
        throw std::invalid_argument("Unknown plane-parallel model '" + model + "'.");
    active_model = model;
    parameters = PlaneParallelParameters<double>{};
    if (model == "Gaensler08") {
        parameters.n0 = 0.014;
        parameters.z0 = 1.83;
    }
}

template <typename T>
T PlaneParallelDensity::field(const double &x, const double &y, const double &z,
                              const PlaneParallelParameters<T> &p) const {
    return p.n0 * exp(-std::abs(z) / p.z0);
}

IMAGINE_INSTANTIATE_SCALAR_MODEL(PlaneParallelDensity)

}
