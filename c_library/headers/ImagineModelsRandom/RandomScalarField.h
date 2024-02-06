#ifndef RANDOMSCALARFIELD_H
#define RANDOMSCALARFIELD_H

#include <vector>
#include <cassert>
#include <functional>
#include <iostream>
#include <algorithm>
#include <complex> 
#include <fftw3.h>
#include <random>
#include <memory>
#include <initializer_list>

#include "exceptions.h"
#include "RandomField.h"

class RandomScalarField : public RandomField<number, double*>  {
protected:
    // Fields
    fftw_plan c2r;
    fftw_plan r2c;

  double* allocate_memory(const std::array<int, 3> shp);

  void free_memory(double* grid_eval);

  fftw_complex* construct_plans(double* grid_eval, const std::array<int, 3> shp);

  void destroy_plans();

public:
  RandomScalarField() : RandomField<number, double*>() {};

  RandomScalarField(std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment);

  ~RandomScalarField();

  // Fields

  const int ndim = 1;
  // methods

  // implemented/hidden in child classes, rms amplitude and anisotropy direction
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  double* random_numbers_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed);

  double* on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed);

  double* on_grid(const int seed);

  void _on_grid(double* val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed);
  
  double* profile_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rfp, const std::array<double, 3> &inc);
};

#endif /* RANDOMSCALARFIELD_H */
