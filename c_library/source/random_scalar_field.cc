#include <cmath>

#include "ImagineModelsRandom/RandomScalarField.h"

namespace imagine {

ScalarGridData RandomScalarField::random_numbers(const RegularGrid &grid, const int seed) const {
    FFTWWorkspace ws(grid.shape);
    ScalarGridData out(grid.shape);
    seed_complex_random_numbers(ws.complex(), grid.shape, grid.increment, seed);
    ws.backward();
    ws.copy_unpadded(out.component(0));
    const double norm = 1. / std::sqrt(double(ws.size()));
    for (double &v : out.data)
        v *= norm;
    return out;
}

ScalarGridData RandomScalarField::sample(const RegularGrid &grid, const int seed) const {
    ScalarGridData out = random_numbers(grid, seed);
    for_each_point(
        grid, [&](std::size_t idx, double x, double y, double z) { out(0, idx) = transform(out(0, idx), x, y, z); });
    return out;
}

}
