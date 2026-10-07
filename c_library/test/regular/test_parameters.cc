#include <set>
#include <stdexcept>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "test_helpers.h"

using namespace imagine;
using namespace imagine::test;

TEST_CASE("UniformMagneticField parameter update", "[parameters]") {
  UniformMagneticField umf;
  CHECK(umf.parameters.bx == 0.);
  CHECK(umf.parameters.by == 0.);
  CHECK(umf.parameters.bz == 0.);
  CHECK(umf.at_position(2.4, 2.1, -.2) == Vec3<double>{0., 0., 0.});

  umf.parameters.bx = -3.2;
  CHECK(umf.get_parameter("bx") == -3.2);
  CHECK(umf.at_position(2.4, 2.1, -.2) == Vec3<double>{-3.2, 0., 0.});

  umf.set_parameter_map({{"by", 1.5}, {"bz", 0.25}});
  CHECK(umf.parameters.by == 1.5);
  CHECK(umf.parameter_map().at("bz") == 0.25);
  CHECK(umf.at_position(2.4, 2.1, -.2) == Vec3<double>{-3.2, 1.5, 0.25});
}

TEST_CASE("Han parameter names follow the declaration order", "[parameters]") {
  auto names = HanMagneticField::parameter_names();
  REQUIRE(names.size() == 9);
  CHECK(names.front() == "B_p");
}

TEMPLATE_LIST_TEST_CASE("parameter registry is consistent", "[parameters]", AllModels) {
  TestType model;
  auto names = TestType::parameter_names();
  REQUIRE(!names.empty());
  CHECK(std::set<std::string>(names.begin(), names.end()).size() == names.size());
  const std::set<std::string> active(model.active_parameters.begin(), model.active_parameters.end());
  CHECK(active.size() == model.active_parameters.size());
  for (const auto &name : model.active_parameters)
    CHECK_NOTHROW(TestType::parameter_index(name));

  auto map = model.parameter_map();
  CHECK(map.size() == names.size());
  for (std::size_t i = 0; i < names.size(); ++i) {
    CAPTURE(names[i]);
    CHECK(TestType::parameter_index(names[i]) == i);
    CHECK(map.at(names[i]) == model.get_parameter(names[i]));
  }
}

TEMPLATE_LIST_TEST_CASE("set_parameter round-trips and only changes that parameter", "[parameters]", AllModels) {
  TestType model;
  const auto before = model.parameter_map();
  for (const auto &name : TestType::parameter_names()) {
    CAPTURE(name);
    const double value = before.at(name) * 1.5 + 0.25;
    model.set_parameter(name, value);
    auto after = model.parameter_map();
    CHECK(after.at(name) == value);
    after[name] = before.at(name);
    CHECK(after == before);
    model.set_parameter(name, before.at(name));
  }
}

TEMPLATE_LIST_TEST_CASE("unknown parameter names throw", "[parameters]", AllModels) {
  TestType model;
  const auto before = model.parameter_map();
  const auto first = TestType::parameter_names().front();
  CHECK_THROWS_AS(model.get_parameter("no_such_parameter"), std::invalid_argument);
  CHECK_THROWS_AS(model.set_parameter("no_such_parameter", 1.), std::invalid_argument);
  CHECK_THROWS_AS(model.set_parameter_map({{first, before.at(first) + 1.}, {"no_such_parameter", 1.}}), std::invalid_argument);
  CHECK(model.parameter_map() == before);
}
