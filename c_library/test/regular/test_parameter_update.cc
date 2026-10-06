#include <iostream>
#include <stdexcept>

#include "ImagineModels/RegularModels.h"

using namespace imagine;

#define check(exp) do { if (!(exp)) { std::cerr << "check failed: " #exp " (" __FILE__ ":" << __LINE__ << ")" << std::endl; std::exit(1); } } while (0)

void test_parameter_update() {
    UniformMagneticField umf;
    check(umf.parameters.bx == 0.);
    check(umf.parameters.by == 0.);
    check(umf.parameters.bz == 0.);

    Vec3<double> zeros{{0., 0., 0.}};
    check(umf.at_position(2.4, 2.1, -.2) == zeros);

    umf.parameters.bx = -3.2;
    check(umf.get_parameter("bx") == -3.2);
    Vec3<double> updated{{-3.2, 0., 0.}};
    check(umf.at_position(2.4, 2.1, -.2) == updated);

    umf.set_parameter_map({{"by", 1.5}, {"bz", 0.25}});
    check(umf.parameters.by == 1.5);
    check(umf.parameter_map().at("bz") == 0.25);
}

void test_parameter_names() {
    check(HanMagneticField::parameter_names().size() == 9);
    check(HanMagneticField::parameter_names()[0] == "B_p");
    bool thrown = false;
    try { UniformMagneticField().set_parameter("bw", 1.); } catch (const std::invalid_argument &) { thrown = true; }
    check(thrown);
    thrown = false;
    try { UniformMagneticField().set_parameter_map({{"bx", 1.}, {"bw", 1.}}); } catch (const std::invalid_argument &) { thrown = true; }
    check(thrown);
}

int main() {
    test_parameter_update();
    test_parameter_names();
}
