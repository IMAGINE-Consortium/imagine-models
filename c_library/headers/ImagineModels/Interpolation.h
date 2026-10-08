#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#include "ImagineModels/Grid.h"

namespace imagine {

enum class Interpolation { linear, nearest };

template <int N>
GridData<N> interpolate(const GridData<N> &data, const RegularGrid &grid, const PointCloud &points,
                        Interpolation method = Interpolation::linear, bool nan_outside = false) {
    if (data.shape != grid.shape)
        throw GridException("interpolate: data shape does not match the grid shape.");
    for (int a = 0; a < 3; ++a)
        if (grid.shape[a] > 1 && grid.increment[a] == 0.)
            throw GridException("interpolate: grid increment must be nonzero.");

    const auto &n = grid.shape;
    const double tolerance = 1e-9;
    GridData<N> out(points.shape());
    for (std::size_t p = 0; p < points.size(); ++p) {
        const double position[3] = {points.x[p], points.y[p], points.z[p]};
        int lower[3];
        double weight[3];
        bool inside = true;
        for (int a = 0; a < 3; ++a) {
            double u = n[a] > 1 ? (position[a] - grid.reference_point[a]) / grid.increment[a] : 0.;
            if (n[a] == 1 && std::abs(position[a] - grid.reference_point[a]) > tolerance * std::abs(grid.increment[a]))
                inside = false;
            if (!std::isfinite(position[a]) || u < -tolerance || u > n[a] - 1 + tolerance)
                inside = false;
            u = std::fmin(std::fmax(u, 0.), double(n[a] - 1));
            if (method == Interpolation::nearest) {
                lower[a] = int(std::lround(u));
                weight[a] = 0.;
            } else {
                lower[a] = std::min(int(std::floor(u)), std::max(n[a] - 2, 0));
                weight[a] = u - lower[a];
            }
        }
        if (!inside) {
            if (!nan_outside)
                throw GridException("interpolate: position outside the grid.");
            for (int c = 0; c < N; ++c)
                out(c, p) = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        for (int c = 0; c < N; ++c) {
            const double *values = data.component(c);
            double sum = 0.;
            for (int corner = 0; corner < 8; ++corner) {
                int index[3];
                double w = 1.;
                for (int a = 0; a < 3; ++a) {
                    const int upper = (corner >> a) & 1;
                    index[a] = lower[a] + upper;
                    w *= upper ? weight[a] : 1. - weight[a];
                }
                if (w == 0.)
                    continue;
                sum += w * values[(std::size_t(index[0]) * n[1] + index[1]) * n[2] + index[2]];
            }
            out(c, p) = sum;
        }
    }
    return out;
}

}
