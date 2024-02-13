#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>


#include "RegularModels.h"
#include "../test_helpers.h"

TEST_CASE("UniformMagneticField") {
    // DEFINITIONS

    // define positions in Galactic cartesian coordinates (units are kpc)

    double posx = -3.5;
    double posy = 2.6;
    double posz = 0.2;

    // define parameters to be updated

    double bx = -3.2;
    double by = 2.279;
    double bz = -0.1;

    // define expected results of at_position
    vector zeros{{0., 0., 0.}};
    vector updated{{bx, by, bz}};

    // define parameters for a regular grid in Galactic cartesian coordinates (units are kpc)
    const std::array<int, 3> shape {{4, 3, 2}};
    const std::array<double, 3> refpoint {{-4., 0.1, -0.3}};
    const std::array<double, 3> increment {{2.1, 0.3, 1.}};

    size_t arr_sz = 4*3*2;

    std::array<double *, 3>  default_regular_grid;  
    default_regular_grid[0] = new double[arr_sz];
    default_regular_grid[1] = new double[arr_sz];
    default_regular_grid[2] = new double[arr_sz];

    std::array<double *, 3>  updated_regular_grid;  
    updated_regular_grid[0] = new double[arr_sz];
    updated_regular_grid[1] = new double[arr_sz];
    updated_regular_grid[2] = new double[arr_sz];

    for (int i=0; i < arr_sz; i++) {
        default_regular_grid[0][i] = 0.; 
        default_regular_grid[1][i] = 0.; 
        default_regular_grid[2][i] = 0.; 
        updated_regular_grid[0][i] = bx; 
        updated_regular_grid[1][i] = by; 
        updated_regular_grid[2][i] = bz; 
    }

    // define positions of a irregular grid in Galactic cartesian coordinates (units are kpc)
    const std::vector<double> grid_x {{2.  , -4. , 0. , 1., 0.4}};
    const std::vector<double> grid_y {{4.  , 6. , -0.1, 0., 0.2}};
    const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};
    
    std::array<double *, 3>  default_irregular_grid;   
    default_irregular_grid[0] = new double[5];
    default_irregular_grid[1] = new double[5];
    default_irregular_grid[2] = new double[5];

    std::array<double *, 3>  updated_irregular_grid;   
    updated_irregular_grid[0] = new double[5];
    updated_irregular_grid[1] = new double[5];
    updated_irregular_grid[2] = new double[5];

    for (int i=0; i < 5; i++) {
        default_regular_grid[0][i] = 0.; 
        default_regular_grid[1][i] = 0.; 
        default_regular_grid[2][i] = 0.; 
        updated_irregular_grid[0][i] = bx; 
        updated_irregular_grid[1][i] = by; 
        updated_irregular_grid[2][i] = bz; 
    }

    SECTION("Empty constructor") {
        // constructor s
        UniformMagneticField umf_plain = UniformMagneticField();
        // test at_position -- default parameters
        REQUIRE_THAT(umf_plain.at_position(0., 0., 0.), EqualsVector(zeros));
        REQUIRE_THAT(umf_plain.at_position(posx, posy, posz), EqualsVector(zeros));
        // test on_grid -- default parameters
        REQUIRE_THROWS_WITH(umf_plain.on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
        REQUIRE_THAT(umf_plain.on_grid(shape, refpoint, increment), EqualsPointerArray(default_regular_grid, arr_sz));
        REQUIRE_THAT(umf_plain.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(default_irregular_grid, 5));
        // parameter update
        umf_plain.bx = bx; 
        umf_plain.by = by;
        umf_plain.bz = bz;  
        REQUIRE(umf_plain.bx == bx); 
        REQUIRE(umf_plain.by == by); 
        REQUIRE(umf_plain.bz == bz); 
        // test at_position -- default parameters
        REQUIRE_THAT(umf_plain.at_position(0., 0., 0.), EqualsVector(updated));
        REQUIRE_THAT(umf_plain.at_position(posx, posy, posz), EqualsVector(updated));
        // test on_grid -- default parameters
        REQUIRE_THROWS_WITH(umf_plain.on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
        REQUIRE_THAT(umf_plain.on_grid(shape, refpoint, increment), EqualsPointerArray(updated_regular_grid, arr_sz));
        REQUIRE_THAT(umf_plain.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(updated_irregular_grid, 5));
        #if autodiff_FOUND
        // test derivative 
        REQUIRE_THAT(umf_plain.derivative(posx, posy, posz), EqualsVector(zeros));
        #endif
    }

    SECTION("Regular grid constructor") {
        // constructor s
        UniformMagneticField umf_regular_grid = UniformMagneticField(shape, refpoint, increment);
        // test at_position -- default parameters
        REQUIRE_THAT(umf_regular_grid.at_position(0., 0., 0.), EqualsVector(zeros));
        REQUIRE_THAT(umf_regular_grid.at_position(posx, posy, posz), EqualsVector(zeros));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_regular_grid.on_grid(), EqualsPointerArray(default_regular_grid, arr_sz));       
        REQUIRE_THAT(umf_regular_grid.on_grid(shape, refpoint, increment), EqualsPointerArray(default_regular_grid, arr_sz));
        REQUIRE_THAT(umf_regular_grid.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(default_irregular_grid, 5));
        // parameter update
        umf_regular_grid.bx = bx; 
        umf_regular_grid.by = by;
        umf_regular_grid.bz = bz;  
        REQUIRE(umf_regular_grid.bx == bx); 
        REQUIRE(umf_regular_grid.by == by); 
        REQUIRE(umf_regular_grid.bz == bz); 
        // test at_position -- default parameters
        REQUIRE_THAT(umf_regular_grid.at_position(0., 0., 0.), EqualsVector(updated));
        REQUIRE_THAT(umf_regular_grid.at_position(posx, posy, posz), EqualsVector(updated));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_regular_grid.on_grid(), EqualsPointerArray(updated_regular_grid, arr_sz));       
        REQUIRE_THAT(umf_regular_grid.on_grid(shape, refpoint, increment), EqualsPointerArray(updated_regular_grid, arr_sz));
        REQUIRE_THAT(umf_regular_grid.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(updated_irregular_grid, 5));
        #if autodiff_FOUND
        // test derivative
        REQUIRE_THAT(umf_regular_grid.derivative(posx, posy, posz), EqualsVector(zeros));
        #endif
    }

    SECTION("Irregular grid constructor") {
        // constructor s
        UniformMagneticField umf_irregular_grid = UniformMagneticField(grid_x, grid_y, grid_z);
        // test at_position -- default parameters
        REQUIRE_THAT(umf_irregular_grid.at_position(0., 0., 0.), EqualsVector(zeros));
        REQUIRE_THAT(umf_irregular_grid.at_position(posx, posy, posz), EqualsVector(zeros));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_irregular_grid.on_grid(), EqualsPointerArray(default_irregular_grid, 5));       
        REQUIRE_THAT(umf_irregular_grid.on_grid(shape, refpoint, increment), EqualsPointerArray(default_regular_grid, arr_sz));
        REQUIRE_THAT(umf_irregular_grid.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(default_irregular_grid, 5));
        // parameter update
        umf_irregular_grid.bx = bx; 
        umf_irregular_grid.by = by;
        umf_irregular_grid.bz = bz;  
        REQUIRE(umf_irregular_grid.bx == bx); 
        REQUIRE(umf_irregular_grid.by == by); 
        REQUIRE(umf_irregular_grid.bz == bz); 
        // test at_position -- default parameters
        REQUIRE_THAT(umf_irregular_grid.at_position(0., 0., 0.), EqualsVector(updated));
        REQUIRE_THAT(umf_irregular_grid.at_position(posx, posy, posz), EqualsVector(updated));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_irregular_grid.on_grid(), EqualsPointerArray(updated_irregular_grid, 5));       
        REQUIRE_THAT(umf_irregular_grid.on_grid(shape, refpoint, increment), EqualsPointerArray(updated_regular_grid, arr_sz));
        REQUIRE_THAT(umf_irregular_grid.on_grid(grid_x, grid_y, grid_z), EqualsPointerArray(updated_irregular_grid, 5));
        #if autodiff_FOUND
        // test derivative
        REQUIRE_THAT(umf_irregular_grid.derivative(posx, posy, posz), EqualsVector(zeros));
        #endif
    }
}
