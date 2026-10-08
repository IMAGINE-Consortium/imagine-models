#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "ImagineModels/ImagineModels.h"

using namespace imagine;
using Catch::Approx;

namespace {

const RegularGrid grid({5, 4, 3}, {-2., 1., -0.5}, {0.5, 0.25, 0.4});

double affine(int c, double x, double y, double z) {
    return 1. + c + 2. * x - 3. * y + 0.5 * z;
}

VectorGridData affine_data(const RegularGrid &g) {
    VectorGridData data(g.shape);
    for_each_point(g, [&](std::size_t idx, double x, double y, double z) {
        for (int c = 0; c < 3; ++c)
            data(c, idx) = affine(c, x, y, z);
    });
    return data;
}

}

TEST_CASE("linear interpolation reproduces affine fields", "[interpolation]") {
    const auto data = affine_data(grid);
    const PointCloud points({-2., -1.3, 0., -0.01, 0.}, {1., 1.6, 1.75, 1.2, 1.3}, {-0.5, 0.1, 0.3, 0.25, -0.5});
    const auto result = interpolate(data, grid, points);
    REQUIRE(result.shape == points.shape());
    for (std::size_t p = 0; p < points.size(); ++p)
        for (int c = 0; c < 3; ++c)
            CHECK(result(c, p) == Approx(affine(c, points.x[p], points.y[p], points.z[p])).epsilon(1e-12));
}

TEST_CASE("interpolation at grid points returns the grid values", "[interpolation]") {
    const auto data = affine_data(grid);
    std::vector<double> x, y, z;
    for_each_point(grid, [&](std::size_t, double a, double b, double c) {
        x.push_back(a);
        y.push_back(b);
        z.push_back(c);
    });
    for (auto method : {Interpolation::linear, Interpolation::nearest}) {
        const auto result = interpolate(data, grid, PointCloud(x, y, z), method);
        for (std::size_t i = 0; i < grid.size(); ++i)
            for (int c = 0; c < 3; ++c)
                CHECK(result(c, i) == Approx(data(c, i)).epsilon(1e-12));
    }
}

TEST_CASE("nearest interpolation returns the closest grid value", "[interpolation]") {
    ScalarGridData data(grid.shape);
    for (std::size_t i = 0; i < grid.size(); ++i)
        data(0, i) = double(i);
    const auto result = interpolate(data, grid, PointCloud({-1.4}, {1.3}, {-0.25}), Interpolation::nearest);
    CHECK(result(0, 0) == double((1 * 4 + 1) * 3 + 1));
}

TEST_CASE("interpolation outside the grid throws or returns NaN", "[interpolation]") {
    const auto data = affine_data(grid);
    const PointCloud outside({0.1, -1.}, {1.2, 1.2}, {0., std::nan("")});
    CHECK_THROWS_AS(interpolate(data, grid, outside), GridException);
    const auto result = interpolate(data, grid, outside, Interpolation::linear, true);
    for (std::size_t p = 0; p < 2; ++p)
        for (int c = 0; c < 3; ++c)
            CHECK(std::isnan(result(c, p)));
}

TEST_CASE("interpolation on a grid with a single layer", "[interpolation]") {
    const RegularGrid flat({3, 2, 1}, {0., 0., 1.}, {1., 1., 1.});
    const auto data = affine_data(flat);
    const auto result = interpolate(data, flat, PointCloud({0.5, 1.5}, {0.5, 1.}, {1., 1.}));
    CHECK(result(0, 0) == Approx(affine(0, 0.5, 0.5, 1.)).epsilon(1e-12));
    CHECK(result(2, 1) == Approx(affine(2, 1.5, 1., 1.)).epsilon(1e-12));
    CHECK_THROWS_AS(interpolate(data, flat, PointCloud({0.5}, {0.5}, {1.2})), GridException);
}

TEST_CASE("interpolation rejects data of another shape", "[interpolation]") {
    const VectorGridData data({4, 4, 3});
    CHECK_THROWS_AS(interpolate(data, grid, PointCloud({0.}, {1.}, {0.})), GridException);
}
