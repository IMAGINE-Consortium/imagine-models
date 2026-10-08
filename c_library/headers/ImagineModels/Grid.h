#pragma once

#include <array>
#include <cstddef>
#include <variant>
#include <vector>

#include "ImagineModels/exceptions.h"

namespace imagine {

struct RegularGrid {
    std::array<int, 3> shape;
    std::array<double, 3> reference_point;
    std::array<double, 3> increment;

    RegularGrid(const std::array<int, 3> &shape, const std::array<double, 3> &reference_point,
                const std::array<double, 3> &increment)
        : shape(shape), reference_point(reference_point), increment(increment) {
        for (int s : shape)
            if (s < 1)
                throw GridException("RegularGrid: all shape entries must be positive.");
    }

    std::size_t size() const { return std::size_t(shape[0]) * shape[1] * shape[2]; }
};

struct IrregularGrid {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    IrregularGrid(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z)
        : x(x), y(y), z(z) {
        if (x.empty() || y.empty() || z.empty())
            throw GridException("IrregularGrid: coordinate vectors must not be empty.");
    }

    std::array<int, 3> shape() const { return {int(x.size()), int(y.size()), int(z.size())}; }
    std::size_t size() const { return x.size() * y.size() * z.size(); }
};

struct PointCloud {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    PointCloud(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z)
        : x(x), y(y), z(z) {
        if (x.empty())
            throw GridException("PointCloud: coordinate vectors must not be empty.");
        if (y.size() != x.size() || z.size() != x.size())
            throw GridException("PointCloud: x, y and z must have the same length.");
    }

    std::array<int, 3> shape() const { return {int(x.size()), 1, 1}; }
    std::size_t size() const { return x.size(); }
};

using Grid = std::variant<RegularGrid, IrregularGrid, PointCloud>;

inline std::array<int, 3> grid_shape(const Grid &grid) {
    if (auto g = std::get_if<RegularGrid>(&grid))
        return g->shape;
    if (auto g = std::get_if<IrregularGrid>(&grid))
        return g->shape();
    return std::get<PointCloud>(grid).shape();
}

template <typename F> void for_each_point(const RegularGrid &grid, F &&f) {
    const auto &size = grid.shape;
    const auto &rpt = grid.reference_point;
    const auto &inc = grid.increment;
    for (int i = 0; i < size[0]; i++) {
        int m = i * size[1] * size[2];
        for (int j = 0; j < size[1]; j++) {
            int n = j * size[2];
            for (int k = 0; k < size[2]; k++)
                f(m + n + k, rpt[0] + i * inc[0], rpt[1] + j * inc[1], rpt[2] + k * inc[2]);
        }
    }
}

template <typename F> void for_each_point(const IrregularGrid &grid, F &&f) {
    const int sx = grid.x.size(), sy = grid.y.size(), sz = grid.z.size();
    for (int i = 0; i < sx; i++) {
        int m = i * sy * sz;
        for (int j = 0; j < sy; j++) {
            int n = j * sz;
            for (int k = 0; k < sz; k++)
                f(m + n + k, grid.x[i], grid.y[j], grid.z[k]);
        }
    }
}

template <typename F> void for_each_point(const PointCloud &grid, F &&f) {
    for (std::size_t i = 0; i < grid.x.size(); i++)
        f(i, grid.x[i], grid.y[i], grid.z[i]);
}

template <typename F> void for_each_point(const Grid &grid, F &&f) {
    std::visit([&f](const auto &g) { for_each_point(g, f); }, grid);
}

template <int N> struct GridData {
    std::array<int, 3> shape;
    std::vector<double> data;

    explicit GridData(const std::array<int, 3> &shape)
        : shape(shape), data(N * std::size_t(shape[0]) * shape[1] * shape[2]) {}

    static constexpr int components = N;
    std::size_t size() const { return std::size_t(shape[0]) * shape[1] * shape[2]; }
    double *component(int c) { return data.data() + c * size(); }
    const double *component(int c) const { return data.data() + c * size(); }
    double &operator()(int c, std::size_t idx) { return data[c * size() + idx]; }
    double operator()(int c, std::size_t idx) const { return data[c * size() + idx]; }
};

using ScalarGridData = GridData<1>;
using VectorGridData = GridData<3>;

template <int N> struct JacobianData {
    std::array<int, 3> shape;
    std::size_t columns;
    std::vector<double> data;

    JacobianData(const std::array<int, 3> &shape, std::size_t columns)
        : shape(shape), columns(columns), data(N * std::size_t(shape[0]) * shape[1] * shape[2] * columns) {}

    static constexpr int components = N;
    std::size_t size() const { return std::size_t(shape[0]) * shape[1] * shape[2]; }
    double &operator()(int c, std::size_t idx, std::size_t col) { return data[(c * size() + idx) * columns + col]; }
    double operator()(int c, std::size_t idx, std::size_t col) const {
        return data[(c * size() + idx) * columns + col];
    }
};

}
