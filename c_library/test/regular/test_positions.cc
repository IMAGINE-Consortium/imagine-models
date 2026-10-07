#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "test_helpers.h"

using namespace imagine;
using namespace imagine::test;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

const Vec3<double> zero{0., 0., 0.};

}

TEST_CASE("JF12 is zero outside its boundaries", "[positions]") {
  JF12MagneticField jf12;
  CHECK(jf12.at_position(0., 0., 0.) == zero);
  CHECK(jf12.at_position(.1, .3, .4) == zero);
  CHECK(jf12.at_position(20.5, 0., 0.) == zero);
  CHECK(jf12.at_position(-15., -15., 2.) == zero);
  CHECK(jf12.at_position(-8.5, 0., 0.) != zero);
}

TEST_CASE("Jaffe and Pshirkov at the Galactic centre", "[positions]") {
  CHECK(JaffeMagneticField().at_position(0., 0., 0.) == zero);
  CHECK(JaffeMagneticField().at_position(0., 0., 1.) == zero);
  CHECK(PshirkovMagneticField().at_position(0., 0., 0.) == zero);
}

TEST_CASE("Helix", "[positions]") {
  HelixMagneticField helix;
  CHECK(helix.at_position(.3, .4, 0.) == zero);
  CHECK(helix.at_position(20., 1., 0.) == zero);

  auto b = helix.at_position(3., 4., .2);
  CHECK_THAT(b[0], WithinRel(.6, 1e-15));
  CHECK_THAT(b[1], WithinRel(.8, 1e-15));
  CHECK(b[2] == 1.);

  helix.parameters = {2., 3., -1.};
  b = helix.at_position(0., -5., 7.);
  CHECK_THAT(b[0], WithinAbs(0., 1e-15));
  CHECK_THAT(b[1], WithinRel(-3., 1e-15));
  CHECK(b[2] == -1.);
}

TEST_CASE("uniform fields are constant", "[positions]") {
  UniformMagneticField b;
  b.parameters = {1., -2., .5};
  UniformDensityField n;
  n.parameters.n0 = .03;
  for (const auto &p : positions) {
    CAPTURE(to_string(p));
    CHECK(b.at_position(p[0], p[1], p[2]) == Vec3<double>{1., -2., .5});
    CHECK(n.at_position(p[0], p[1], p[2]) == .03);
  }
}
