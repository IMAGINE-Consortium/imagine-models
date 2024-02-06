#ifndef RANDOMFIELD_H
#define RANDOMFIELD_H

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
#include "Field.h"



template<typename POSTYPE, typename GRIDTYPE>
class RandomField : public Field<POSTYPE, GRIDTYPE>  {
protected:

  bool created_fftw_plans = false;
  bool has_fftw_wisdom = false;
  bool no_profile = false;

public:
  // Constructors
  using Field<POSTYPE, GRIDTYPE> :: Field;

  // Fields

  bool apply_spectrum = true;


  // methods
  POSTYPE at_position(const double &x, const double &y, const double &z) const {
    throw NotImplementedException();
    // Here comes the interpolator
  }
  
  virtual double spatial_profile(const double &x, const double &y, const double &z) const = 0;

  virtual double calculate_fourier_sigma(const double &abs_k, const double &dk) const = 0;

    // This function is the place where the global routine should be implemented, i.e. how the spatial profile modifies the random field, and if divergence cleaning needs to be performed. 
  virtual void _on_grid(GRIDTYPE val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) = 0;

  void seed_complex_random_numbers(fftw_complex* vec,  const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed); 

  GRIDTYPE on_grid(const std::vector<double>  &grid_x, const std::vector<double>  &grid_y, const std::vector<double>  &grid_z, const int seed) {
    throw NotImplementedException();
  }

  virtual double* profile_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rfp, const std::array<double, 3> &inc) = 0;

  void remove_padding(double* val, const std::array<int, 3> &shp, const int pad);

  double simple_spectrum(const double &abs_k, const double &dk, const double &k0, const double &s) const;

  double hammurabi_spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const;

};

//include the random field method implementations (at this position, due to the use of templates)
#include "random.tpp"

#endif /* RANDOMFIELD_H */
