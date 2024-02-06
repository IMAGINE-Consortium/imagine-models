#ifndef RANDOMVECTORFIELD_H
#define RANDOMVECTORFIELD_H

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
#include "RegularField.h"


class RandomVectorField : public RandomField<vector, std::array<double*, 3>>  {
protected:
  // protected fields

  std::array<fftw_plan, 3> r2c;
  std::array<fftw_plan, 3> c2r;

  // protected member functions 

  // fftw plan (con/de)struction
  std::array<fftw_complex*, 3> construct_plans(std::array<double*, 3> grid_eval, std::array<int, 3> shp);
  void destroy_plans();

public:
  // constructors
  RandomVectorField() : RandomField() {};

  RandomVectorField(std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment);

  // destructor
  ~RandomVectorField();

  // FIELDS

  const int ndim = 3;
  bool clean_divergence = true;  
  bool apply_anisotropy = true;  

  double anisotropy_rho = 1.;

  // METHODS

  // memory (de-)allocation
  std::array<double*, 3> allocate_memory(std::array<int, 3> shp);
  void free_memory(std::array<double*, 3> grid_eval);
  
  // implemented/hidden in child classes, rms amplitude and anisotropy direction
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;
  
  vector anisotropy_direction(const double &x, const double &y, const double &z) const {
    vector a{{0., 0., 0.}}; 
    return a;
  }

  // interpolator (TBD)
  vector at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
  }

  // internal on_grid function, combining seeding of random numbers, applying the spatial profile and divergence cleaning
  void _on_grid(std::array<double*, 3> val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed); 

  // divergence cleaner
  void divergence_cleaner(fftw_complex* bx, fftw_complex* by, fftw_complex* bz,  const std::array<int, 3> &shp, const std::array<double, 3> &inc) const;

  std::array<double*, 3> random_numbers_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed);

  std::array<double*, 3> on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed);

  std::array<double*, 3> on_grid(const int seed);

  double* profile_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rfp, const std::array<double, 3> &inc);
};


#endif /* RANDOMVECTORFIELD_H */
