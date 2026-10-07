#include "test_helpers.h"

#if IMAGINE_HAS_AUTODIFF

#include <algorithm>
#include <type_traits>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace imagine;
using namespace imagine::test;

namespace {

std::vector<Position> off_axis_positions() {
  std::vector<Position> out;
  for (const auto &p : positions)
    if (p[0] != 0. || p[1] != 0.)
      out.push_back(p);
  return out;
}

template <typename Model>
std::vector<double> finite_difference(Model &model, const std::string &name, const Position &p, double h) {
  const double p0 = model.get_parameter(name);
  model.set_parameter(name, p0 + h);
  auto up = value_at(model, p);
  model.set_parameter(name, p0 - h);
  auto down = value_at(model, p);
  model.set_parameter(name, p0);
  for (std::size_t k = 0; k < up.size(); ++k)
    up[k] = (up[k] - down[k]) / (2 * h);
  return up;
}

bool all_close(const std::vector<double> &a, const std::vector<double> &b, double rtol, double atol) {
  for (std::size_t k = 0; k < a.size(); ++k)
    if (std::abs(a[k] - b[k]) > atol + rtol * std::abs(b[k]))
      return false;
  return true;
}

bool all_finite(const Eigen::MatrixXd &m) { return m.allFinite(); }

}

TEMPLATE_LIST_TEST_CASE("Jacobian has one column per active parameter", "[derivatives]", AllModels) {
  TestType model;
  const int rows = std::is_same_v<decltype(model.at_position(0., 0., 0.)), double> ? 1 : 3;
  auto jac = model.derivative(-8.5, 1., .2);
  CHECK(jac.rows() == rows);
  CHECK(jac.cols() == Eigen::Index(model.active_parameters.size()));
}

TEMPLATE_LIST_TEST_CASE("Jacobian matches finite differences", "[derivatives]", AllModels) {
  TestType model;
  for (const auto &p : off_axis_positions()) {
    CAPTURE(to_string(p));
    auto jac = model.derivative(p[0], p[1], p[2]);
    REQUIRE(all_finite(jac));
    for (std::size_t c = 0; c < model.active_parameters.size(); ++c) {
      const auto &name = model.active_parameters[c];
      CAPTURE(name);
      const double scale = std::max(1., std::abs(model.get_parameter(name)));
      auto fine = finite_difference(model, name, p, 1e-6 * scale);
      auto coarse = finite_difference(model, name, p, 1e-4 * scale);
      if (!all_close(fine, coarse, 1e-2, 1e-6))
        continue;
      double largest = 1.;
      for (double v : fine)
        largest = std::max(largest, std::abs(v));
      std::vector<double> column(jac.rows());
      for (Eigen::Index k = 0; k < jac.rows(); ++k)
        column[k] = jac(k, c);
      CAPTURE(column, fine);
      CHECK(all_close(column, fine, 1e-4, 1e-6 * largest));
    }
  }
}

TEMPLATE_LIST_TEST_CASE("Jacobian columns follow active_parameters", "[derivatives]", AllModels) {
  TestType model;
  const auto all = model.active_parameters;
  const Position p{-8.5, 1., .2};
  const auto full = model.derivative(p[0], p[1], p[2]);
  std::vector<std::string> subset{all.back()};
  if (all.size() > 1)
    subset.push_back(all.front());
  model.active_parameters = subset;
  const auto partial = model.derivative(p[0], p[1], p[2]);
  REQUIRE(partial.cols() == Eigen::Index(subset.size()));
  for (std::size_t c = 0; c < subset.size(); ++c) {
    const auto full_column = std::find(all.begin(), all.end(), subset[c]) - all.begin();
    CHECK(partial.col(c) == full.col(full_column));
  }
  model.active_parameters = {"no_such_parameter"};
  CHECK_THROWS_AS(model.derivative(p[0], p[1], p[2]), std::invalid_argument);
}

TEMPLATE_LIST_TEST_CASE("Jacobian is finite on the z-axis", "[derivatives]", AllModels) {
  TestType model;
  for (const auto &p : z_axis) {
    CAPTURE(to_string(p));
    CHECK(all_finite(model.derivative(p[0], p[1], p[2])));
  }
}

#endif
