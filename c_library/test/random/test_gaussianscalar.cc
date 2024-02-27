#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>


#include "RandomModels.h"
#include "../test_helpers.h"

TEST_CASE("GaussianScalarField") {

    const int n_seeds = 10;

    std::array<int,  n_seeds> seeds{23, 35, 645, 1, 75, 968876 , 323424, 89987, 86786, 2342};
    
    UNSCOPED_INFO("Start testing gaussian scalar test case");

    //const std::array<int, 3> shape {{150, 200, 600}};
    const std::array<int, 3> shape {{4, 4, 3}};
    const std::array<double, 3> refpoint {{-4., -0.1, -3.2}};
    const std::array<double, 3> increment {{0.5, 0.01, 1.6}};

    int size = shape[0]*shape[1]*shape[2];

    const std::vector<double> grid_x {{2.  , -4. , 0. , 1., 0.4}};
    const std::vector<double> grid_y {{4.  , 6. , -0.1, 0., 0.2}};
    const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};
    

    GaussianScalarField gauss_plain = GaussianScalarField();

    gauss_plain.apply_spectrum = false;

    double model_var = gauss_plain.rms * gauss_plain.rms;

    double variance_of_the_sample_variance = 2 * model_var / (size - 2);
    double variance_of_the_sample_mean = model_var/ (size - 1);

    std::cout << "size: " << size << std::endl;
    std::cout << "sqrt size: " << std::sqrt(size) << std::endl;

    std::cout << "variance_of_the_sample_mean: " << variance_of_the_sample_mean << std::endl;
    std::cout << "variance_of_the_sample_variance: " << variance_of_the_sample_variance << std::endl;

    for (int i = 0; i < n_seeds; i++) {
        double sample_mean = 0.;
        double sample_variance = 0.;

        fftw_complex* ev_random_c = gauss_plain.test_random_numbers(shape, increment, seeds[i]);
        double* ev_random = gauss_plain.random_numbers_on_grid(shape, increment, seeds[i]);

        for (int j = 0; j < size - 1; j++) {
            sample_mean = sample_mean + ev_random[j];
            //std::cout << j << ": " << ev_random[j] << ", " << sample_mean << std::endl;
        }

        sample_mean = sample_mean / (size - 1);

        CHECK_THAT(sample_mean, Catch::Matchers::WithinAbs(gauss_plain.mean, 2.* std::sqrt(variance_of_the_sample_mean)) );    

        for (int j = 0; j < size - 1; j++) {
            double diff = ev_random[j] - gauss_plain.mean;
            sample_variance = sample_variance + diff*diff;
        }
        sample_variance = sample_variance / (size - 1) ; /// M_PI;

        std::cout << "sample_mean: " << sample_mean << std::endl;
        std::cout << "sample_variance: " << sample_variance << std::endl;

        CHECK_THAT(sample_variance, Catch::Matchers::WithinAbs(model_var, 2.*std::sqrt(variance_of_the_sample_variance)) ); 
    }
  

}