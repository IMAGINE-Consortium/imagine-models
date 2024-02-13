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

Eigen::Matrix3d model_derivative(double x, double y, number ampx, number ampy, number ampz, double rmin, double rmax) {
    const double phi = std::atan2(y, x);         // azimuthal angle in cylindrical coordinates
    const double r = std::sqrt(x * x + y * y); // radius in cylindrical coordinates
    Eigen::Matrix3d b_deriv;
    if ((r > rmin) && (r < rmax))
    {
        b(0, 0) = std::cos(phi);
        b(0, 1) = 0.;
        b(0, 2) = 0.;
        
        b(1, 0) = 0.;
        b(1, 1) = std::sin(phi);
        b(1, 2) = 0.;
        
        b(2, 0) = 0;
        b(2, 1) = 0.;
        b(2, 2) = 1.;
    } 
    return b_deriv;
}
# endif

TEST_CASE("HelixMagneticField") {
    // DEFINITIONS

    // define positions in Galactic cartesian coordinates (units are kpc)
    double posx = -3.5:
    double posy = 2.6;
    double posz = 0.2;

    std::array<double, 3> origin{{0., 0., 0.}};
    std::array<double, 3> smaller_rmin{{.9, -0.87, 0.65}};
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
    const std::vector<double> grid_x {{2.  , -4. , 0. , 1., -32.4}};
    const std::vector<double> grid_y {{4.  , 6. , -0.1, 0., 0.2}};
    const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 12.}};
    
    // constructor s
    HelixMagneticField helix = HelixMagneticField();
    // test at_position -- default parameters
    REQUIRE_THAT(helix.at_position(origin[0], origin[1], origin[2]), EqualsVector(zeros));
    REQUIRE_THAT(helix.at_position(smaller_rmin[0], smaller_rmin[1], smaller_rmin[2]), EqualsVector(zeros));
    REQUIRE_THAT(helix.at_position(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], bigger_rmin_smaller_rmax[2]), EqualsVector(model(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], helix.ampx, helix.ampy, helix.ampz, helix.rmin, helix.rmax)));
    REQUIRE_THAT(helix.at_position(bigger_rmax[0], bigger_rmax[1], bigger_rmax[2]), EqualsVector(zeros));

    // test on_grid -- default parameters
    REQUIRE_THROWS_WITH(helix.on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
    REQUIRE_NO_THROW(helix.on_grid(shape, refpoint, increment));
    REQUIRE_NO_THROW(helix.on_grid(grid_x, grid_y, grid_z));
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
    REQUIRE_THAT(helix.at_position(smaller_rmin[0], smaller_rmin[1], smaller_rmin[2]), EqualsVector(zeros));
    REQUIRE_THAT(helix.at_position(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], bigger_rmin_smaller_rmax[2]), EqualsVector(model(bigger_rmin_smaller_rmax[0], bigger_rmin_smaller_rmax[1], ampx, ampy, ampz, rmin, rmax)));
    REQUIRE_THAT(helix.at_position(bigger_rmax[0], bigger_rmax[1], bigger_rmax[2]), EqualsVector(zeros));
    // test on_grid -- default parameters
    REQUIRE_THROWS_WITH(helix.on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
    REQUIRE_THAT(helix.on_grid(shape, refpoint, increment), EqualsPointerArray(updated_regular_grid, arr_sz));
    REQUIRE_THAT(helix.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(updated_irregular_grid, 5));
    #if autodiff_FOUND
    // test derivative

    REQUIRE_THAT(helix.derivative(posx, posy, posz), EqualsVector(zeros));
    #endif
    
}
