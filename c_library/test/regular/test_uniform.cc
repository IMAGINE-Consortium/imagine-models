#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "RegularModels.h"

#define assertm(exp, msg) assert(((void)msg, exp))


void _check_array_equality_from_pointer(std::array<double*, 3> a, std::array<double*, 3> b , size_t &n) {
    for (int d = 0; d==3; ++d) {
        std::vector<double>  arr_a(a[d], a[d] + n);
        std::vector<double>  arr_b(b[d], b[d] + n);
        assert (arr_a == arr_b); 
    }


void test_uniform() {

    // DEFINITIONS

    // define positions in Galactic cartesian coordinates (units are kpc)
    vector zeros{{0., 0., 0.}};
    vector position{{-3.5, 2.6, 0.2}};
    
    // define parameters to be updated

    double bx = -3.2;
    double by = 2.279;
    double bz = -0.1;

    vector updated{{bx, by, bz}};

    // define parameters for a regular grid in Galactic cartesian coordinates (units are kpc)
    const std::array<int, 3> shape {{4, 3, 2}};
    const std::array<double, 3> refpoint {{-4., 0.1, -0.3}};
    const std::array<double, 3> increment {{2.1, 0.3, 1.}};

    size_t arr_sz = 4*3*2;

    const std::array<double *, 3>  default_regular_grid;  
    default_regular_grid[0] = new double[arr_sz];
    default_regular_grid[1] = new double[arr_sz];
    default_regular_grid[2] = new double[arr_sz];

    const std::array<double *, 3>  updated_regular_grid;  
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
    
    const std::array<double *, 3>  default_irregular_grid;   
    default_irregular_grid[0] = new double[5];
    default_irregular_grid[1] = new double[5];
    default_irregular_grid[2] = new double[5];

    const std::array<double *, 3>  updated_irregular_grid;   
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

    // constructor 
    
    UniformMagneticField umf_plain = UniformMagneticField();
    UniformMagneticField umf_regular_grid = UniformMagneticField(shape, refpoint, increment);
    UniformMagneticField umf_irregular_grid = UniformMagneticField(grid_x, grid_y, grid_z);

    // TESTING
    
    // test at_position -- default parameters 

    assert (umf_plain.at_position(zeros[0], zeros[1], zeros[2]) == zeros);
    assert (umf_plain.at_position(position[0], position[1], position[2]) == zeros);
    assert (umf_regular_grid.at_position(zeros[0], zeros[1], zeros[2]) == zeros);
    assert (umf_regular_grid.at_position(position[0], position[1], position[2]) == zeros);
    assert (umf_irregular_grid.at_position(zeros[0], zeros[1], zeros[2]) == zeros);
    assert (umf_irregular_grid.at_position(position[0], position[1], position[2]) == zeros);
    
    // test on_grid -- default parameters 

    // assert (umf_plain.on_grid() == zeros); need to catch excpetion here

    std::array<double * , 3> plain_regular = umf_plain.on_grid(shape, refpoint, increment);
    _check_array_equality_from_pointer(plain_regular, default_regular_grid);

    std::array<double * , 3> plain_irregular = umf_plain.on_grid(grid_x, grid_y, grid_z);
    _check_array_equality_from_pointer(plain_irregular, default_irregular_grid);

    std::array<double * , 3> regular_noinput = umf_regular_grid.on_grid();
    _check_array_equality_from_pointer(regular_noinput, default_regular_grid);

    std::array<double * , 3> irregular_noinput = umf_irregular_grid.on_grid();
    _check_array_equality_from_pointer(irregular_noinput, default_irregular_grid);

    // test parameter update
    
    umf_plain.bx = bx; 
    umf_plain.by = by;
    umf_plain.bz = bz;  
    
    assert (umf_plain.bx == bx); 
    assert (umf_plain.by == by); 
    assert (umf_plain.bz == bz); 

    // test at_position -- updated parameters 
    
    assert (umf_plain.at_position(zeros[0], zeros[1], zeros[2]) == updated);
    assert (umf_plain.at_position(position[0], position[1], position[2]) == updated);


    // test on_grid -- updated parameters 

    std::array<double * , 3> plain_regular_updated = umf_plain.on_grid(shape, refpoint, increment);
    _check_array_equality_from_pointer(plain_regular_updated, updated_regular_grid);

    std::array<double * , 3> plain_irregular_updated = umf_plain.on_grid(grid_x, grid_y, grid_z);
    _check_array_equality_from_pointer(plain_irregular_updated, updated_irregular_grid);


    // test derivative
}