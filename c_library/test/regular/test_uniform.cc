#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "RegularModels.h"

#define assertm(exp, msg) assert(((void)msg, exp))


// Custom matcher to check array equality
template <typename T, size_t N>
struct ArrayEqualsMatcher : Catch::Matchers::Impl::MatcherBase<std::array<T, N>> {
    ArrayEqualsMatcher(const std::array<T, N>& expected) : m_expected(expected) {}

    bool match(const std::array<T, N>& arr) const override {
        for (size_t i = 0; i < N; ++i) {
            if (arr[i] != m_expected[i]) {
                return false;
            }
        }
        return true;
    }

    std::string describe() const override {
        std::ostringstream oss;
        oss << "equals: [";
        for (size_t i = 0; i < N; ++i) {
            if (i > 0) oss << ", ";
            oss << m_expected[i];
        }
        oss << "]";
        return oss.str();
    }

    const std::array<T, N>& m_expected;
};

// Helper function to create the matcher
template <typename T, size_t N>
ArrayEqualsMatcher<T, N> EqualsArray(const std::array<T, N>& expected) {
    return ArrayEqualsMatcher<T, N>(expected);
}


void _check_array_equality_from_pointer(std::array<double*, 3> a, std::array<double*, 3> b , size_t &n) {
    for (int d = 0; d==3; ++d) {
        std::vector<double>  arr_a(a[d], a[d] + n);
        std::vector<double>  arr_b(b[d], b[d] + n);
        assert (arr_a == arr_b); 
    }
}



TEST_CASE("UniformMagneticField") {
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

    SECTION("Empty constructor") {
        // constructor s
        UniformMagneticField umf_plain = UniformMagneticField();
        // test at_position -- default parameters
        REQUIRE_THAT(umf_plain.at_position(zeros[0], zeros[1], zeros[2]), EqualsArray(zeros));
        REQUIRE_THAT(umf_plain.at_position(position[0], position[1], position[2]), EqualsArray(zeros));
        // test on_grid -- default parameters
        REQUIRE_THROWS(umf_plain.on_grid(), GridException);
        REQUIRE_THAT(umf_plain.on_grid(shape, refpoint, increment), EqualsArray(default_regular_grid));
        REQUIRE_THAT( umf_plain.on_grid(grid_x, grid_y, grid_z), EqualsArray(default_irregular_grid));
        // parameter update
        umf_plain.bx = bx; 
        umf_plain.by = by;
        umf_plain.bz = bz;  
        REQUIRE(umf_plain.bx == bx); 
        REQUIRE(umf_plain.by == by); 
        REQUIRE(umf_plain.bz == bz); 
        // test at_position -- default parameters
        REQUIRE_THAT(umf_plain.at_position(zeros[0], zeros[1], zeros[2]), EqualsArray(updated));
        REQUIRE_THAT(umf_plain.at_position(position[0], position[1], position[2]), EqualsArray(updated));
        // test on_grid -- default parameters
        REQUIRE_THROWS(umf_plain.on_grid(), GridException);
        REQUIRE_THAT(umf_plain.on_grid(shape, refpoint, increment), EqualsArray(updated_regular_grid));
        REQUIRE_THAT(umf_plain.on_grid(grid_x, grid_y, grid_z), EqualsArray(updated_irregular_grid));
        #if autodiff_FOUND
        // test derivative 
        REQUIRE_THAT(umf_plain.derivative(position[0], position[1], position[2]), EqualsArray(zeros));
        #endif
    }

    SECTION("Regular grid constructor") {
        // constructor s
        UniformMagneticField umf_regular_grid = UniformMagneticField(shape, refpoint, increment);
        // test at_position -- default parameters
        REQUIRE_THAT(umf_regular_grid.at_position(zeros[0], zeros[1], zeros[2]), EqualsArray(zeros));
        REQUIRE_THAT(umf_regular_grid.at_position(position[0], position[1], position[2]), EqualsArray(zeros));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_regular_grid.on_grid(), EqualsArray(default_regular_grid));       
        REQUIRE_THAT(umf_regular_grid.on_grid(shape, refpoint, increment), EqualsArray(default_regular_grid));
        REQUIRE_THAT(umf_regular_grid.on_grid(grid_x, grid_y, grid_z), EqualsArray(default_irregular_grid));
        // parameter update
        umf_regular_grid.bx = bx; 
        umf_regular_grid.by = by;
        umf_regular_grid.bz = bz;  
        REQUIRE(umf_regular_grid.bx == bx); 
        REQUIRE(umf_regular_grid.by == by); 
        REQUIRE(umf_regular_grid.bz == bz); 
        // test at_position -- default parameters
        REQUIRE_THAT(umf_regular_grid.at_position(zeros[0], zeros[1], zeros[2]), EqualsArray(updated));
        REQUIRE_THAT(umf_regular_grid.at_position(position[0], position[1], position[2]), EqualsArray(updated));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_regular_grid.on_grid(), EqualsArray(updated_regular_grid));       
        REQUIRE_THAT(umf_regular_grid.on_grid(shape, refpoint, increment), EqualsArray(updated_regular_grid));
        REQUIRE_THAT(umf_regular_grid.on_grid(grid_x, grid_y, grid_z), EqualsArray(updated_irregular_grid));
        #if autodiff_FOUND
        // test derivative
        REQUIRE_THAT(umf_regular_grid.derivative(position[0], position[1], position[2]), EqualsArray(zeros));
        #endif
    }

    SECTION("Irregular grid constructor") {
        // constructor s
        UniformMagneticField umf_irregular_grid = UniformMagneticField(grid_x, grid_y, grid_z);
        // test at_position -- default parameters
        REQUIRE_THAT(umf_irregular_grid.at_position(zeros[0], zeros[1], zeros[2]), EqualsArray(zeros));
        REQUIRE_THAT(umf_irregular_grid.at_position(position[0], position[1], position[2]), EqualsArray(zeros));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_irregular_grid.on_grid(), EqualsArray(default_irregular_grid));       
        REQUIRE_THAT(umf_irregular_grid.on_grid(shape, refpoint, increment), EqualsArray(default_regular_grid));
        REQUIRE_THAT(umf_irregular_grid.on_grid(grid_x, grid_y, grid_z), EqualsArray(default_irregular_grid));
        // parameter update
        umf_irregular_grid.bx = bx; 
        umf_irregular_grid.by = by;
        umf_irregular_grid.bz = bz;  
        REQUIRE(umf_irregular_grid.bx == bx); 
        REQUIRE(umf_irregular_grid.by == by); 
        REQUIRE(umf_irregular_grid.bz == bz); 
        // test at_position -- default parameters
        REQUIRE_THAT(umf_irregular_grid.at_position(zeros[0], zeros[1], zeros[2]), EqualsArray(updated));
        REQUIRE_THAT(umf_irregular_grid.at_position(position[0], position[1], position[2]), EqualsArray(updated));
        // test on_grid -- default parameters
        REQUIRE_THAT(umf_irregular_grid.on_grid(), EqualsArray(updated_irregular_grid));       
        REQUIRE_THAT(umf_irregular_grid.on_grid(shape, refpoint, increment), EqualsArray(updated_regular_grid));
        REQUIRE_THAT(umf_irregular_grid.on_grid(grid_x, grid_y, grid_z), EqualsArray(updated_irregular_grid));
        #if autodiff_FOUND
        // test derivative
        REQUIRE_THAT(umf_irregular_grid.derivative(position[0], position[1], position[2]), EqualsArray(zeros));
        #endif
    }
}
