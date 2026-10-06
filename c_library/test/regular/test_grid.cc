#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "ImagineModels/RegularModels.h"

using namespace imagine;

#define check(exp) do { if (!(exp)) { std::cerr << "check failed: " #exp " (" __FILE__ ":" << __LINE__ << ")" << std::endl; std::exit(1); } } while (0)


IrregularGrid as_irregular(const RegularGrid &g) {
    std::vector<double> x, y, z;
    for (int i = 0; i < g.shape[0]; ++i) x.push_back(g.reference_point[0] + i * g.increment[0]);
    for (int j = 0; j < g.shape[1]; ++j) y.push_back(g.reference_point[1] + j * g.increment[1]);
    for (int k = 0; k < g.shape[2]; ++k) z.push_back(g.reference_point[2] + k * g.increment[2]);
    return IrregularGrid(x, y, z);
}

void test_grid(const std::map<std::string, std::shared_ptr<RegularVectorField>> &models, const RegularGrid &regular, const IrregularGrid &irregular) {
    for (const auto &[name, model] : models) {
        VectorGridData eval_regular = model->evaluate(regular);
        VectorGridData eval_as_irregular = model->evaluate(as_irregular(regular));
        check(eval_regular.shape == regular.shape);
        check(eval_regular.data == eval_as_irregular.data);

        VectorGridData eval_irregular = model->evaluate(irregular);
        check(eval_irregular.shape == irregular.shape());
        size_t idx = 0;
        for (double x : irregular.x)
            for (double y : irregular.y)
                for (double z : irregular.z) {
                    Vec3<double> v = model->at_position(x, y, z);
                    for (int c = 0; c < 3; ++c)
                        check(eval_irregular(c, idx) == static_cast<double>(v[c]));
                    ++idx;
                }
    }
}

void test_scalar_grid(const RegularGrid &regular) {
    YMW16 ymw;
    ScalarGridData eval_regular = ymw.evaluate(regular);
    ScalarGridData eval_as_irregular = ymw.evaluate(as_irregular(regular));
    check(eval_regular.data == eval_as_irregular.data);
}

void test_invalid_grids() {
    bool thrown = false;
    try { RegularGrid({4, 0, 2}, {0., 0., 0.}, {1., 1., 1.}); } catch (const GridException &) { thrown = true; }
    check(thrown);
    thrown = false;
    try { IrregularGrid({1., 2.}, {}, {0.}); } catch (const GridException &) { thrown = true; }
    check(thrown);
}


int main() {
    const IrregularGrid irregular({2., 4., 0., 1., .4}, {4., 6., 0.1, 0., .2}, {-0.2, 0.8, 0.2, 0., 1.});
    const RegularGrid regular({4, 3, 2}, {-4., 0.1, -0.3}, {2.1, 0.3, 1.});

    std::map<std::string, std::shared_ptr<RegularVectorField>> models;
    models["JF12"] = std::make_shared<JF12MagneticField>();
    models["Jaffe"] = std::make_shared<JaffeMagneticField>();
    models["Helix"] = std::make_shared<HelixMagneticField>();
    models["UF24"] = std::make_shared<UFMagneticField>();

    test_grid(models, regular, irregular);
    test_scalar_grid(regular);
    test_invalid_grids();
}
