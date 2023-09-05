#ifndef REGULAR_TRAMPOLINE_H
#define REGULAR_TRAMPOLINE_H

#include "hamunits.h"
#include "Field.h"

#include "RegularField.h"

#include "../../c_library/headers/autodiff.hh"

#if defined autodiff_FOUND
    #include <autodiff/forward/real.hpp>
    #include <autodiff/forward/dual.hpp>
    #include <autodiff/forward/real/eigen.hpp>
    namespace ad = autodiff;
    typedef ad::real number;
    typedef ad::VectorXreal vector;
#else
    typedef double number;  // only used for differentiable numbers! 
    typedef std::array<double, 3> vector;
#endif


#include <iostream>

//using Array3Type = std::array<double, 3>;
using Array3PointerType = std::array<double*, 3>; // Only for PYBIND11_OVERRIDE_PURE macro, else gets confused by commas 

// These classes are necessary to override virtual functions when binding abstract c++ classes

class PyScalarFieldBase: public Field<double, double*> {
public:
    using Field<double, double*>:: Field; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, Field, at_position, x, y, z); }

    double* on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, seed); }
    
    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(double*, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }


    double* allocate_memory(std::array<int, 3> shp) override {PYBIND11_OVERRIDE_PURE(double*, Field, allocate_memory, shp); }
    
    void free_memory(double* grid_eval) override {PYBIND11_OVERRIDE_PURE(void, Field, free_memory, grid_eval); }
    
};


class PyVectorFieldBase: public Field<vector, std::array<double*, 3>> {
public:
    using Field<vector, std::array<double*, 3>>:: Field; // Inherit constructors
    vector at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(vector, Field, at_position, x, y, z); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }


    std::array<double*, 3> allocate_memory(std::array<int, 3> shp) override {PYBIND11_OVERRIDE_PURE(Array3PointerType, Field, allocate_memory, shp); }
    
    void free_memory(std::array<double*, 3> grid_eval) override {PYBIND11_OVERRIDE_PURE(void, Field, free_memory, grid_eval); }
    
};


class PyRegularVectorField: public RegularVectorField {
public:
    using RegularVectorField:: RegularVectorField; // Inherit constructors
    vector at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(vector, RegularVectorField, at_position, x, y, z); }

    std::array<double*, 3> on_grid(int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RegularVectorField, on_grid, seed); }
    
    std::array<double*, 3> on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RegularVectorField, on_grid, grid_x, grid_y, grid_z, seed); }

    std::array<double*, 3> on_grid(const std::array<int, 3>& grid_shape, const std::array<double, 3>& grid_zeropoint, const std::array<double, 3>& grid_increment, int seed) override {PYBIND11_OVERRIDE(Array3PointerType, RegularVectorField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }
    
};


class PyRegularScalarField : public RegularScalarField {
public:
    using RegularScalarField:: RegularScalarField; // Inherit constructors
    double at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(double, RegularScalarField, at_position, x, y, z); }
    
    double* on_grid(const int seed=0) override {PYBIND11_OVERRIDE_PURE(double*, RegularScalarField, on_grid,  seed); }

    double* on_grid(const std::vector<double>& grid_x, const std::vector<double>& grid_y, const std::vector<double>& grid_z, const int seed=0) override {PYBIND11_OVERRIDE_PURE(double*, RegularScalarField, on_grid, grid_x, grid_y, grid_z, seed); }

    double* on_grid(const std::array<int, 3> &grid_shape, const std::array<double, 3> &grid_zeropoint, const std::array<double, 3> &grid_increment, const int seed = 0) override {PYBIND11_OVERRIDE_PURE(double*, RegularScalarField, on_grid, grid_shape, grid_zeropoint, grid_increment, seed); }
    
};


#endif