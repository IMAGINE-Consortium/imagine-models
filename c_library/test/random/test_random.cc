#include <cmath>
#include <complex>
#include <memory>
#include <numeric>
#include <vector>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "ImagineModelsRandom/RandomModels.h"
#include "test_helpers.h"

using namespace imagine;
using namespace imagine::test;
using Catch::Matchers::WithinRel;

namespace {

const RegularGrid small_grid({16, 16, 16}, {-1., -1., -1.}, {.125, .125, .125});
const RegularGrid stat_grid({48, 32, 24}, {-8., -4., -2.}, {.25, .25, .25});
const int n_seeds = 12;

class ConstantRandomField : public RandomVectorField {
public:
  double value = 2.;
  double slope = 0.;

  explicit ConstantRandomField(double slope = 0.) : slope(slope) { apply_spectrum = slope != 0.; }
  double rms(const double &, const double &, const double &) const override { return value; }
  double spectrum(const double &k) const override { return std::pow(k, -slope); }
};

class VerticalAnisotropy : public ConstantRandomField {
public:
  Vec3<double> anisotropy_direction(const double &, const double &, const double &) const override { return {0., 0., 3.}; }
};

void check_within(const std::vector<double> &values, double expected, double n_sigma = 5.) {
  const double n = values.size();
  const double mean = std::accumulate(values.begin(), values.end(), 0.) / n;
  double var = 0.;
  for (double v : values)
    var += (v - mean) * (v - mean);
  const double error = std::sqrt(var / (n - 1) / n);
  CAPTURE(mean, expected, error);
  CHECK(std::abs(mean - expected) < n_sigma * std::max(error, 1e-12));
}

template <int N>
double mean_square(const GridData<N> &g) {
  double s = 0.;
  for (double v : g.data)
    s += v * v;
  return s / g.size();
}

double mean_square(const VectorGridData &g, int c) {
  double s = 0.;
  for (std::size_t i = 0; i < g.size(); ++i)
    s += g(c, i) * g(c, i);
  return s / g.size();
}

double frequency(int i, int n, double d) { return (i < (n + 1) / 2 ? i : i - n) / (n * d); }

double relative_divergence(const VectorGridData &b, const RegularGrid &grid) {
  const auto [n0, n1, n2] = grid.shape;
  const int m2 = n2 / 2 + 1;
  std::vector<std::vector<std::complex<double>>> spectra;
  for (int c = 0; c < 3; ++c) {
    std::vector<double> in(b.component(c), b.component(c) + b.size());
    std::vector<std::complex<double>> out(std::size_t(n0) * n1 * m2);
    fftw_plan plan = fftw_plan_dft_r2c_3d(n0, n1, n2, in.data(), reinterpret_cast<fftw_complex *>(out.data()), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    spectra.push_back(std::move(out));
  }
  double divergence = 0., scale = 0.;
  for (int i = 0; i < n0; ++i)
    for (int j = 0; j < n1; ++j)
      for (int l = 0; l < m2; ++l) {
        const std::size_t idx = (std::size_t(i) * n1 + j) * m2 + l;
        const double k[3] = {frequency(i, n0, grid.increment[0]), frequency(j, n1, grid.increment[1]), l / (n2 * grid.increment[2])};
        std::complex<double> kb = 0.;
        double b2 = 0.;
        for (int c = 0; c < 3; ++c) {
          kb += k[c] * spectra[c][idx];
          b2 += std::norm(spectra[c][idx]);
        }
        divergence += std::abs(kb);
        scale += std::sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]) * std::sqrt(b2);
      }
  return divergence / scale;
}

}

using RandomModels = std::tuple<JF12RandomField, ESRandomField, GaussianScalarField, LogNormalScalarField>;
using RandomVectorModels = std::tuple<JF12RandomField, ESRandomField>;

TEMPLATE_LIST_TEST_CASE("samples have the grid shape and are finite", "[random]", RandomModels) {
  TestType model;
  auto field = model.sample(small_grid, 3);
  CHECK(field.shape == small_grid.shape);
  CHECK(all_finite(field.data));
  CHECK(mean_square(field) > 0.);
}

TEMPLATE_LIST_TEST_CASE("samples are reproducible per seed", "[random]", RandomModels) {
  TestType model;
  CHECK(model.sample(small_grid, 7).data == model.sample(small_grid, 7).data);
  CHECK(model.sample(small_grid, 7).data != model.sample(small_grid, 8).data);
}

TEMPLATE_LIST_TEST_CASE("random numbers have zero mean and unit variance", "[random]", RandomModels) {
  TestType model;
  model.apply_spectrum = GENERATE(false, true);
  CAPTURE(model.apply_spectrum);
  std::vector<double> variances;
  for (int seed = 0; seed < n_seeds; ++seed) {
    auto g = model.random_numbers(stat_grid, seed);
    for (int c = 0; c < decltype(g)::components; ++c) {
      double mean = 0.;
      for (std::size_t i = 0; i < g.size(); ++i)
        mean += g(c, i);
      CHECK(std::abs(mean / g.size()) < 1e-10);
    }
    variances.push_back(mean_square(g));
  }
  check_within(variances, 1.);
}

TEST_CASE("GaussianScalarField realises mu and sigma", "[random]") {
  GaussianScalarField model;
  model.mu = 2.5;
  model.sigma = .3;
  model.apply_spectrum = false;
  CHECK(model.mean(0., 0., 0.) == 2.5);
  CHECK(model.rms(0., 0., 0.) == .3);
  std::vector<double> means, variances;
  for (int seed = 0; seed < n_seeds; ++seed) {
    auto s = model.sample(stat_grid, seed);
    const double mean = std::accumulate(s.data.begin(), s.data.end(), 0.) / s.size();
    double var = 0.;
    for (double v : s.data)
      var += (v - mean) * (v - mean);
    means.push_back(mean);
    variances.push_back(var / s.size());
  }
  check_within(means, 2.5);
  check_within(variances, .09);
}

TEST_CASE("LogNormalScalarField realises its mean and rms", "[random]") {
  LogNormalScalarField model;
  model.log_mu = .2;
  model.log_sigma = .5;
  model.apply_spectrum = false;
  CHECK_THAT(model.mean(0., 0., 0.), WithinRel(std::exp(.2 + .125), 1e-14));
  CHECK_THAT(model.rms(0., 0., 0.), WithinRel(std::sqrt(std::expm1(.25)) * std::exp(.2 + .125), 1e-14));
  std::vector<double> means, variances;
  for (int seed = 0; seed < n_seeds; ++seed) {
    auto s = model.sample(stat_grid, seed);
    const double mean = std::accumulate(s.data.begin(), s.data.end(), 0.) / s.size();
    double var = 0.;
    for (double v : s.data)
      var += (v - mean) * (v - mean);
    means.push_back(mean);
    variances.push_back(var / s.size());
  }
  check_within(means, model.mean(0., 0., 0.));
  check_within(variances, model.variance(0., 0., 0.));
}

TEST_CASE("vector amplitude matches rms", "[random]") {
  ConstantRandomField model(GENERATE(0., 2.));
  model.clean_divergence = GENERATE(false, true);
  CAPTURE(model.slope, model.clean_divergence);
  std::vector<double> energies;
  for (int seed = 0; seed < n_seeds; ++seed)
    energies.push_back(mean_square(model.sample(stat_grid, seed)));
  check_within(energies, 4.);
}

TEMPLATE_LIST_TEST_CASE("model amplitude follows the rms profile", "[random]", RandomVectorModels) {
  TestType model;
  model.clean_divergence = false;
  model.apply_spectrum = false;
  const auto rms = model.evaluate_rms(stat_grid);
  std::vector<double> ratios;
  for (int seed = 0; seed < n_seeds; ++seed) {
    auto b = model.sample(stat_grid, seed);
    double sum = 0.;
    int count = 0;
    for (std::size_t i = 0; i < b.size(); ++i)
      if (rms(0, i) > 1e-3) {
        sum += (b(0, i) * b(0, i) + b(1, i) * b(1, i) + b(2, i) * b(2, i)) / (rms(0, i) * rms(0, i));
        ++count;
      }
    ratios.push_back(sum / count);
  }
  check_within(ratios, 1.);
}

TEMPLATE_LIST_TEST_CASE("divergence cleaning preserves total power", "[random]", RandomVectorModels) {
  TestType model;
  const auto rms = model.evaluate_rms(stat_grid);
  const double rms2 = mean_square(rms);
  std::vector<double> ratios;
  for (int seed = 0; seed < n_seeds; ++seed)
    ratios.push_back(mean_square(model.sample(stat_grid, seed)) / rms2);
  check_within(ratios, 1.);
}

TEST_CASE("divergence cleaning removes the divergence", "[random]") {
  std::unique_ptr<RandomVectorField> model;
  SECTION("constant rms") { model = std::make_unique<ConstantRandomField>(); }
  SECTION("constant rms with spectrum") { model = std::make_unique<ConstantRandomField>(2.); }
  SECTION("JF12") { model = std::make_unique<JF12RandomField>(); }
  SECTION("ES") { model = std::make_unique<ESRandomField>(); }
  model->clean_divergence = true;
  CHECK(relative_divergence(model->sample(stat_grid, 4), stat_grid) < 1e-12);
  model->clean_divergence = false;
  CHECK(relative_divergence(model->sample(stat_grid, 4), stat_grid) > .3);
}

TEST_CASE("anisotropy sets the parallel power share", "[random]") {
  VerticalAnisotropy model;
  model.clean_divergence = false;
  model.anisotropy_rho = GENERATE(1., 2., .5);
  model.apply_anisotropy = GENERATE(true, false);
  CAPTURE(model.anisotropy_rho, model.apply_anisotropy);
  const double rho4 = std::pow(model.anisotropy_rho, 4);
  const double expected_share = model.apply_anisotropy ? rho4 / (rho4 + 2.) : 1. / 3.;
  std::vector<double> energies, shares;
  for (int seed = 0; seed < n_seeds; ++seed) {
    auto b = model.sample(stat_grid, seed);
    const double energy = mean_square(b);
    energies.push_back(energy);
    shares.push_back(mean_square(b, 2) / energy);
  }
  check_within(energies, 4.);
  check_within(shares, expected_share);
}

TEMPLATE_LIST_TEST_CASE("evaluate_rms equals rms at each point", "[random]", RandomModels) {
  TestType model;
  const IrregularGrid grid({-8.5, 3.}, {0., 4.}, {.1, -.4});
  const auto on_grid = model.evaluate_rms(grid);
  REQUIRE(on_grid.shape == grid.shape());
  for_each_point(Grid(grid), [&](std::size_t idx, double x, double y, double z) {
    CHECK(on_grid(0, idx) == model.rms(x, y, z));
    CHECK(model.variance(x, y, z) == model.rms(x, y, z) * model.rms(x, y, z));
  });
}

TEST_CASE("JF12 anisotropy follows the regular JF12 field", "[random]") {
  const JF12RandomField random;
  const JF12MagneticField regular;
  for (const auto &p : positions) {
    CAPTURE(to_string(p));
    const auto a = random.anisotropy_direction(p[0], p[1], p[2]);
    const auto b = regular.at_position(p[0], p[1], p[2]);
    for (int c = 0; c < 3; ++c)
      CHECK_THAT(a[c], WithinRel(b[c], 1e-12) || Catch::Matchers::WithinAbs(b[c], 1e-14));
  }
}
