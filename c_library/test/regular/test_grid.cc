#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "test_helpers.h"

using namespace imagine;
using namespace imagine::test;

namespace {

const RegularGrid regular({4, 3, 2}, {-4., 0.1, -0.3}, {2.1, 0.3, 1.});
const IrregularGrid irregular({2., 4., 0., 1., .4, -12.}, {4., 6., 0.1, 0., .2}, {-0.2, 0.8, 0.2, 0., 1.});

IrregularGrid as_irregular(const RegularGrid &g) {
    std::vector<double> x, y, z;
    for (int i = 0; i < g.shape[0]; ++i)
        x.push_back(g.reference_point[0] + i * g.increment[0]);
    for (int j = 0; j < g.shape[1]; ++j)
        y.push_back(g.reference_point[1] + j * g.increment[1]);
    for (int k = 0; k < g.shape[2]; ++k)
        z.push_back(g.reference_point[2] + k * g.increment[2]);
    return IrregularGrid(x, y, z);
}

}

TEMPLATE_LIST_TEST_CASE("evaluate on a regular grid equals the same irregular grid", "[grid]", AllModels) {
    TestType model;
    auto on_regular = model.evaluate(regular);
    auto on_irregular = model.evaluate(as_irregular(regular));
    CHECK(on_regular.shape == regular.shape);
    CHECK(on_regular.data.size() == decltype(on_regular)::components * regular.size());
    CHECK(on_regular.data == on_irregular.data);
}

TEMPLATE_LIST_TEST_CASE("evaluate on an irregular grid equals at_position", "[grid]", AllModels) {
    TestType model;
    auto result = model.evaluate(irregular);
    REQUIRE(result.shape == irregular.shape());
    std::size_t idx = 0;
    for (double x : irregular.x)
        for (double y : irregular.y)
            for (double z : irregular.z) {
                auto expected = value_at(model, {x, y, z});
                for (std::size_t c = 0; c < expected.size(); ++c)
                    CHECK(result(int(c), idx) == expected[c]);
                ++idx;
            }
}

TEMPLATE_LIST_TEST_CASE("evaluate on a point cloud equals at_position", "[grid]", AllModels) {
    TestType model;
    std::vector<double> x, y, z;
    for (const auto &p : positions) {
        x.push_back(p[0]);
        y.push_back(p[1]);
        z.push_back(p[2]);
    }
    const PointCloud cloud(x, y, z);
    auto result = model.evaluate(cloud);
    REQUIRE(result.shape == std::array<int, 3>{int(positions.size()), 1, 1});
    for (std::size_t i = 0; i < positions.size(); ++i) {
        auto expected = value_at(model, positions[i]);
        for (std::size_t c = 0; c < expected.size(); ++c)
            CHECK(result(int(c), i) == expected[c]);
    }
}

TEST_CASE("invalid grids throw", "[grid]") {
    CHECK_THROWS_AS(RegularGrid({4, 0, 2}, {0., 0., 0.}, {1., 1., 1.}), GridException);
    CHECK_THROWS_AS(RegularGrid({-1, 2, 2}, {0., 0., 0.}, {1., 1., 1.}), GridException);
    CHECK_THROWS_AS(IrregularGrid({1., 2.}, {}, {0.}), GridException);
    CHECK_THROWS_AS(PointCloud({}, {}, {}), GridException);
    CHECK_THROWS_AS(PointCloud({1., 2.}, {1.}, {1., 2.}), GridException);
}

TEST_CASE("grid shapes and sizes", "[grid]") {
    CHECK(regular.size() == 24);
    CHECK(irregular.shape() == std::array<int, 3>{6, 5, 5});
    CHECK(grid_shape(Grid(regular)) == regular.shape);
    CHECK(grid_shape(Grid(irregular)) == irregular.shape());
    const PointCloud cloud({1., 2., 3.}, {0., 0., 0.}, {-1., 0., 1.});
    CHECK(cloud.size() == 3);
    CHECK(grid_shape(Grid(cloud)) == std::array<int, 3>{3, 1, 1});
}
