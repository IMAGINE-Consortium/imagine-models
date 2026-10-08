#include <stdexcept>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include "test_helpers.h"

using namespace imagine;
using namespace imagine::test;

TEMPLATE_LIST_TEST_CASE("values are finite", "[models]", AllModels) {
    TestType model;
    for (const auto &p : positions) {
        CAPTURE(to_string(p));
        CHECK(all_finite(value_at(model, p)));
    }
}

TEMPLATE_LIST_TEST_CASE("values are finite on the z-axis", "[models]", AllModels) {
    TestType model;
    for (const auto &p : z_axis) {
        CAPTURE(to_string(p));
        CHECK(all_finite(value_at(model, p)));
    }
}

TEST_CASE("UF24 variants load their published parameters", "[models][variants]") {
    UFMagneticField uf;
    for (const auto &variant : uf.available_models) {
        CAPTURE(variant);
        uf.set_parameter("fDiskB1", 123.);
        uf.set_model(variant);
        CHECK(uf.model() == variant);
        CHECK(uf.parameter_map() == UFMagneticField(variant).parameter_map());
        for (const auto &[name, value] : uf.all_parameters.at(variant))
            CHECK(uf.get_parameter(name) == value);
        for (const auto &p : positions)
            CHECK(all_finite(value_at(uf, p)));
    }
    CHECK_THROWS_AS(uf.set_model("no_such_model"), std::invalid_argument);
    CHECK_THROWS_AS(UFMagneticField("no_such_model"), std::invalid_argument);
}

TEST_CASE("TF17 variants load their published parameters", "[models][variants]") {
    const auto disk = GENERATE(as<std::string>{}, "Ad1", "Bd1", "Dd1");
    const auto halo = GENERATE(as<std::string>{}, "C0", "C1");
    CAPTURE(disk, halo);
    TFMagneticField tf("Ad1", "C0");
    tf.set_parameter("B1_disk", 123.);
    tf.set_model(disk, halo);
    const TFMagneticField fresh(disk, halo);
    CHECK(tf.disk_model() == disk);
    CHECK(tf.halo_model() == halo);
    CHECK(tf.parameter_map() == fresh.parameter_map());
    CHECK(tf.active_parameters == fresh.active_parameters);
    for (const auto &p : positions)
        CHECK(all_finite(value_at(tf, p)));
}

TEST_CASE("TF17 rejects unknown variants", "[models][variants]") {
    TFMagneticField tf;
    CHECK_THROWS_AS(tf.set_model("Xd1", "C0"), std::invalid_argument);
    CHECK_THROWS_AS(tf.set_model("Ad1", "C9"), std::invalid_argument);
    CHECK(tf.disk_model() == "Ad1");
    CHECK(tf.halo_model() == "C0");
}

TEST_CASE("switches change the field", "[models]") {
    const Position p{-8.5, 1., 1.5};
    JF12MagneticField jf12;
    const auto full = value_at(jf12, p);
    jf12.do_halo = false;
    CHECK(value_at(jf12, p) != full);
    jf12.do_halo = true;
    jf12.do_X = false;
    CHECK(value_at(jf12, p) != full);

    SVT22MagneticField svt;
    const auto with_halo = value_at(svt, p);
    svt.do_halo = false;
    CHECK(value_at(svt, p) != with_halo);
}

TEST_CASE("Jaffe rejects unsupported arm numbers", "[models]") {
    JaffeMagneticField jaffe;
    jaffe.arm_num = 3;
    CHECK(all_finite(value_at(jaffe, {-8.5, 1., .2})));
    jaffe.arm_num = 5;
    CHECK_THROWS_AS(jaffe.at_position(-8.5, 1., .2), std::invalid_argument);
    jaffe.arm_num = 1;
    CHECK_THROWS_AS(jaffe.at_position(-8.5, 1., .2), std::invalid_argument);
}

TEST_CASE("Pshirkov variants load their published parameters", "[models][variants]") {
    PshirkovMagneticField ps;
    CHECK(ps.model() == "BSS");
    CHECK(ps.parameters.pitch == -6.);
    CHECK(ps.parameters.B0_Hs == 4.);
    ps.set_parameter("B0_D", 123.);
    ps.set_model("ASS");
    CHECK(ps.model() == "ASS");
    CHECK(ps.parameters.pitch == -5.);
    CHECK(ps.parameters.B0_Hs == 2.);
    CHECK(ps.parameter_map() == PshirkovMagneticField("ASS").parameter_map());
    ps.set_model("BSS");
    CHECK(ps.parameter_map() == PshirkovMagneticField().parameter_map());
    CHECK_THROWS_AS(ps.set_model("no_such_model"), std::invalid_argument);
    ps.useDisk = false;
    ps.useHalo = false;
    CHECK(ps.at_position(-8.5, 1., .2) == Vec3<double>{0., 0., 0.});
}
