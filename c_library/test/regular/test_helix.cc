#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>


#include "RegularModels.h"
#include "../test_helpers.h"


vector model(double x, double y, number ampx, number ampy, number ampz, double rmin, double rmax) {
    const double phi = std::atan2(y, x);         // azimuthal angle in cylindrical coordinates
    const double r = std::sqrt(x * x + y * y); // radius in cylindrical coordinates
    vector b{{0.0, 0.0, 0.0}};
    if ((r > rmin) && (r < rmax))
    {
        b[0] = std::cos(phi) * ampx;
        b[1] = std::sin(phi) * ampy;
        b[2] = ampz;
    } 
    return b;
}

#if autodiff_FOUND

Eigen::MatrixXd model_derivative(double x, double y, number ampx, number ampy, number ampz, double rmin, double rmax) {
    const double phi = std::atan2(y, x);         // azimuthal angle in cylindrical coordinates
    const double r = std::sqrt(x * x + y * y); // radius in cylindrical coordinates
    if ((r < rmin) || (r > rmax)) return Eigen::MatrixXd::Zero(3, 3);
    Eigen::MatrixXd b_deriv;

    b_deriv(0, 0) = std::cos(phi);
    b_deriv(0, 1) = 0.;
    b_deriv(0, 2) = 0.;
    
    b_deriv(1, 0) = 0.;
    b_deriv(1, 1) = std::sin(phi);
    b_deriv(1, 2) = 0.;
    
    b_deriv(2, 0) = 0;
    b_deriv(2, 1) = 0.;
    b_deriv(2, 2) = 1.;
    
    return b_deriv;
}
# endif

TEST_CASE("HelixMagneticField") {

    // constructor s
    HelixMagneticField helix = HelixMagneticField();
    // DEFINITIONS

    // define positions in Galactic cartesian coordinates (units are kpc)

    std::array<double, 3> origin{{0., 0., 0.}};
    std::array<double, 3> smaller_rmin{{.09, -0.47, 0.65}};
    std::array<double, 3> bigger_rmin_smaller_rmax{{-17., 13.4, 9.5}};
    std::array<double, 3> bigger_rmax{{24., -23.13, 19.5}};

    // define parameters to be updated
    
    number ampx = -3.;
    number ampy = 3.;
    number ampz = 13.4;
    double rmax = 2.;
    double rmin = .5;

    // define expected results of at_position
    vector zeros{{0., 0., 0.}};

    // define parameters for a regular grid in Galactic cartesian coordinates (units are kpc)
    const std::array<int, 3> shape {{4, 3, 2}};
    const std::array<double, 3> refpoint {{-4., 0.1, -0.3}};
    const std::array<double, 3> increment {{2.1, 0.3, 1.}};

    size_t arr_sz = 4*3*2;

    // define positions of a irregular grid in Galactic cartesian coordinates (units are kpc)
    const std::vector<double> grid_x {{0., 0. , 0.   , 0.  , 0.  , 0.3, -17.2, 
    0.     , 0.    , 0.   , 0.    , -35.032, 2.2, -6. , .234,  -0.00424, -2.2, -2.32 , .234,  
    0.43 , 8.23    , 3.   , 1., -0.001   , -5.342, -5.23, -2.432  }};
    const std::vector<double> grid_y {{0., 0. , 0.   , 11.5, -7.2, 0. , 0.   , 
    -0.0032, 23.353, -6.1 , 4.    , 0.     , 0. , 0.  , 0.  ,  -0.00132, 2.2 ,  -6.42, .234,   
    6.123, 8.1     , -3.  , -1., 0.0033  , 0.45  , -30.2, -1.543  }};
    const std::vector<double> grid_z {{0., 3.1, -0.01, 0.  , 0.  , 0. , 0.   , 
    -5.32  , 0.012 , 34.02, -8.431, -5.132, 9.3 , -6.1, 4.  , 0.       , 0.  , 0.    , 0.  ,
    10.3 , -0.02   , 3.   , -1., 0.000432, -15.43, 0.001, -7.123  }};

    const int size_pos = grid_x.size();

    std::array<double *, 3>  default_irregular_grid;   
    default_irregular_grid[0] = new double[size_pos];
    default_irregular_grid[1] = new double[size_pos];
    default_irregular_grid[2] = new double[size_pos];

    std::array<double *, 3>  updated_irregular_grid;   
    updated_irregular_grid[0] = new double[size_pos];
    updated_irregular_grid[1] = new double[size_pos];
    updated_irregular_grid[2] = new double[size_pos];

    for (int i=0; i < size_pos; i++) {
        vector def_mod_at_grid_point = model(grid_x[0], grid_y[0], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax);
        default_irregular_grid[0][i] = def_mod_at_grid_point[0].val(); 
        default_irregular_grid[1][i] = def_mod_at_grid_point[1].val(); 
        default_irregular_grid[2][i] = def_mod_at_grid_point[2].val(); 
        vector up_mod_at_grid_point = model(grid_x[0], grid_y[0], ampx, ampy, ampz, rmin, rmax);
        updated_irregular_grid[0][i] = up_mod_at_grid_point[0].val(); 
        updated_irregular_grid[1][i] = up_mod_at_grid_point[1].val(); 
        updated_irregular_grid[2][i] = up_mod_at_grid_point[2].val(); 
    }

    // test at_position -- default parameters
    REQUIRE_THAT(helix.at_position(origin[0], origin[1], origin[2]), EqualsVector(zeros));
    REQUIRE_THAT(helix.at_position(smaller_rmin[0], smaller_rmin[1], smaller_rmin[2]), EqualsVector(zeros));
    REQUIRE_THAT(helix.at_position(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], bigger_rmin_smaller_rmax[2]), EqualsVector(model(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax)));
    REQUIRE_THAT(helix.at_position(bigger_rmax[0], bigger_rmax[1], bigger_rmax[2]), EqualsVector(zeros));

    // test on_grid -- default parameters
    REQUIRE_THROWS_WITH(helix.on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
    REQUIRE_NOTHROW(helix.on_grid(shape, refpoint, increment));
    REQUIRE_THAT(helix.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(default_irregular_grid, size_pos));

    #if autodiff_FOUND
    // test derivative
    CAPTURE(helix.derivative(origin[0], origin[1], origin[2]), model_derivative(origin[0], origin[1], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax));
    REQUIRE_THAT(helix.derivative(origin[0], origin[1], origin[2]), EqualsMatrix(model_derivative(origin[0], origin[1], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax)));

    REQUIRE_THAT(helix.derivative(smaller_rmin[0], smaller_rmin[1], smaller_rmin[2]), EqualsMatrix(model_derivative(smaller_rmin[0], smaller_rmin[1], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax)));
    REQUIRE_THAT(helix.derivative(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], bigger_rmin_smaller_rmax[2]), EqualsMatrix(model_derivative(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax)));
    REQUIRE_THAT(helix.derivative(bigger_rmax[0], bigger_rmax[1], bigger_rmax[2]), EqualsMatrix(model_derivative(bigger_rmax[0], bigger_rmax[1], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax)));
    #endif
   

    // parameter update
    helix.ampx = ampx; 
    helix.ampy = ampy;
    helix.ampz = ampz;  
    helix.rmin = rmin;
    helix.rmax = rmax;  
    REQUIRE(helix.ampx == ampx); 
    REQUIRE(helix.ampy == ampy); 
    REQUIRE(helix.ampz == ampz); 
    REQUIRE(helix.rmin == rmin); 
    REQUIRE(helix.rmax == rmax); 

    
    // test at_position -- default parameters
    REQUIRE_THAT(helix.at_position(origin[0], origin[1], origin[2]), EqualsVector(zeros));
    REQUIRE_THAT(helix.at_position(smaller_rmin[0], smaller_rmin[1], smaller_rmin[2]), EqualsVector(model(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], ampx, ampy, ampz, rmin, rmax)));
    REQUIRE_THAT(helix.at_position(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], bigger_rmin_smaller_rmax[2]), EqualsVector(zeros));
    REQUIRE_THAT(helix.at_position(bigger_rmax[0], bigger_rmax[1], bigger_rmax[2]), EqualsVector(zeros));
    // test on_grid -- default parameters
    REQUIRE_THROWS_WITH(helix.on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
    REQUIRE_NOTHROW(helix.on_grid(shape, refpoint, increment));
    REQUIRE_THAT(helix.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(updated_irregular_grid, size_pos));
    
    
    #if autodiff_FOUND
    // test derivative

    REQUIRE_THAT(helix.derivative(origin[0], origin[1], origin[2]), EqualsMatrix(model_derivative(origin[0], origin[1], ampx, ampy, ampz, rmin, rmax)));
    REQUIRE_THAT(helix.derivative(smaller_rmin[0], smaller_rmin[1], smaller_rmin[2]), EqualsMatrix(model_derivative(smaller_rmin[0], smaller_rmin[1], ampx, ampy, ampz, rmin, rmax)));
    REQUIRE_THAT(helix.derivative(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], bigger_rmin_smaller_rmax[2]), EqualsMatrix(model_derivative(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], ampx, ampy, ampz, rmin, rmax)));
    REQUIRE_THAT(helix.derivative(bigger_rmax[0], bigger_rmax[1], bigger_rmax[2]), EqualsMatrix(model_derivative(bigger_rmax[0], bigger_rmax[1], ampx, ampy, ampz, rmin, rmax)));
    #endif
    
}