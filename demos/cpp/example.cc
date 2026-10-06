#include <iostream>
#include <map>
#include <memory>
#include <string>

#include "ImagineModels/ImagineModels.h"

using namespace imagine;

void print_positions(const std::map<std::string, std::shared_ptr<RegularVectorField>> &models,
                     const std::map<std::string, std::array<double, 3>> &positions) {
  for (const auto &[name, model] : models) {
    std::cout << "The model " << name << " is evaluated:\n";
    for (const auto &[label, p] : positions) {
      Vec3<double> b = model->at_position(p[0], p[1], p[2]);
      std::cout << "  " << label << " (" << p[0] << ", " << p[1] << ", " << p[2] << ") kpc: "
                << b[0] << " " << b[1] << " " << b[2] << " muG\n";
    }
  }
}

template <int N>
void print_grid(const std::string &name, const GridData<N> &data) {
  std::cout << name << " on a " << data.shape[0] << "x" << data.shape[1] << "x" << data.shape[2] << " grid:\n";
  for (std::size_t idx = 0; idx < data.size(); ++idx) {
    std::cout << "  " << idx << ":";
    for (int c = 0; c < N; ++c)
      std::cout << " " << data(c, idx);
    std::cout << "\n";
  }
}

int main() {
  std::map<std::string, std::shared_ptr<RegularVectorField>> models;
  models["Jansson Farrar regular"] = std::make_shared<JF12MagneticField>();
  models["Jaffe"] = std::make_shared<JaffeMagneticField>();
  models["Helix"] = std::make_shared<HelixMagneticField>();

  std::map<std::string, std::array<double, 3>> positions;
  positions["origin"] = {0., 0., 0.};
  positions["sun"] = {-8.5, 0., 0.};
  positions["above plane"] = {3., -2., 1.5};
  print_positions(models, positions);

  RegularGrid regular({4, 3, 2}, {-4., 0.1, -0.3}, {2.1, 0.3, 1.});
  IrregularGrid irregular({2., 4., 0.}, {4., 6., 0.1}, {-0.2, 0.8});

  print_grid("Jansson Farrar regular", models["Jansson Farrar regular"]->evaluate(regular));
  print_grid("Helix", models["Helix"]->evaluate(irregular));
  print_grid("YMW16", YMW16().evaluate(irregular));

#if IMAGINE_HAS_FFTW
  print_grid("Gaussian random field", GaussianScalarField().sample(RegularGrid({2, 2, 2}, {0., 0., 0.}, {1., 1., 1.}), 13));
  print_grid("Jansson Farrar random", JF12RandomField().sample(RegularGrid({2, 2, 2}, {-8.5, 0., 0.}, {1., 1., 1.}), 13));
#endif
}
