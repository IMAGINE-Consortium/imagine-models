#include <cmath>
#include <iostream>

#include "RandomScalarField.h"

RandomScalarField::RandomScalarField(std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment) : RandomField<number, double*>(shape, reference_point, increment) {
  //accumulate wisdom
  double* grid_eval = allocate_memory(shape);
  fftw_complex* grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
  fftw_plan r2c_temp = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], grid_eval, grid_eval_comp, FFTW_MEASURE);
  fftw_plan c2r_temp = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], grid_eval_comp, grid_eval, FFTW_MEASURE);
  //save wisdom
  const char *filename = "ImagineModelsRandomScalarField";
  int fftw_export_wisdom_to_filename(*filename);
  //destroy plans, but wisdom will be kept!
  has_fftw_wisdom = true;
  fftw_destroy_plan(c2r_temp);
  fftw_destroy_plan(r2c_temp); 
  fftw_free(grid_eval);
};

RandomScalarField::~RandomScalarField() {
//delete wisdom
  destroy_plans();
  fftw_forget_wisdom();
  #ifdef _OPENMP
    fftw_cleanup_threads();
  #else
    fftw_cleanup();
  #endif
};

double* RandomScalarField::allocate_memory(const std::array<int, 3> shp) {
  int newshp2;
  if (shp[2] % 2) {
    newshp2 = shp[2] + 1;
  }
  else {
    newshp2 = shp[2] + 2;
  }
  double* grid_eval = (double*) fftw_alloc_real(shp[0]*shp[1]*newshp2);
  return grid_eval;  
};  

void RandomScalarField::free_memory(double* grid_eval) {
  fftw_free(grid_eval);
}

fftw_complex* RandomScalarField::construct_plans(double* grid_eval, const std::array<int, 3> shp) {
  fftw_complex* grid_eval_comp = reinterpret_cast<fftw_complex*>(grid_eval);
  if (has_fftw_wisdom) {
    const char *filename = "ImagineModelsRandomFieldScalar"; // FIX this, currently overwrites
    int fftw_im_wisdom_to_filename(*filename);
  }
  r2c = fftw_plan_dft_r2c_3d(shp[0], shp[1], shp[2], grid_eval, grid_eval_comp, FFTW_ESTIMATE);
  c2r = fftw_plan_dft_c2r_3d(shp[0], shp[1], shp[2], grid_eval_comp, grid_eval,  FFTW_ESTIMATE);
  created_fftw_plans = true;
  return grid_eval_comp;
}

void RandomScalarField::destroy_plans() {
  if (created_fftw_plans) {
    fftw_destroy_plan(c2r);
    fftw_destroy_plan(r2c);
  }
}

double* RandomScalarField::on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {
  double* grid_eval = allocate_memory(shp);
  _on_grid(grid_eval, shp, rpt, inc, seed);
  return grid_eval;
}

double* RandomScalarField::on_grid(const int seed) {
  if (not initialized_with_grid) 
    throw GridException();
  double* grid_eval = allocate_memory(internal_shape);
  const char *filename = "ImagineModelsRandomScalarField";
  int fftw_import_wisdom_from_filename(*filename);
  _on_grid(grid_eval, internal_shape, internal_ref_point, internal_increment, seed);
  return grid_eval;
}

double* RandomScalarField::profile_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rfp, const std::array<double, 3> &inc) {
  double* grid_eval;
  size_t arr_sz = shp[0]*shp[1]*shp[2];
  grid_eval = new double[arr_sz];
  evaluate_function_on_grid<number, double*>(grid_eval, shp, rfp, inc,
                                    [this](double xx, double yy, double zz)
                                    { return spatial_profile(xx, yy, zz); });
  return grid_eval;
} 

double* RandomScalarField::random_numbers_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed) {
    double* val = allocate_memory(shp);
    fftw_complex* val_comp = construct_plans(val, shp); 
    int gs = grid_size(shp);
    double sqrt_gs = std::sqrt(gs);
    std::array<int, 3> padded_shp = {shp[0],  shp[1],  2*(shp[2]/2 + 1)}; 
    int padded_size = grid_size(padded_shp);
    int pad =  padded_shp[2] - shp[2];

    seed_complex_random_numbers(val_comp, shp, inc, seed);
    fftw_execute(c2r);
    for (int s = 0; s < padded_size; ++s)  {
        (val)[s] /= sqrt_gs;  
    }
    remove_padding(val, shp, pad);
    
    return val;
}

void RandomScalarField::_on_grid(double* val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {

  fftw_complex* val_comp = construct_plans(val, shp); 
  std::array<int, 3> padded_shp = {shp[0],  shp[1],  2*(shp[2]/2 + 1)}; 
  int gs = grid_size(shp);
  double sqrt_gs = std::sqrt(gs);
  int padded_size = grid_size(padded_shp);
  int pad =  padded_shp[2] - shp[2];
  // Step 1: draw random numbers with variance 1, possibly correlated

  seed_complex_random_numbers(val_comp, shp, inc, seed);
  fftw_execute(c2r);
  
  
  // Step 2: apply spatial amplitude, possibly introduce anisotropy depending on regular field.
  if (!no_profile) {
    auto apply_profile = [&](double &rand_val, const double xx, const double yy, const double zz) {
        
      double sp = spatial_profile(xx, yy, zz);
      // apply profile
      rand_val *= sp;
      return rand_val;      
    };
  
    apply_function_to_field<double*, double>(val, padded_shp, rpt, inc, apply_profile);
  }

  for (int s = 0; s < padded_size; ++s)  {
    (val)[s] /= sqrt_gs;  
  }
  remove_padding(val, shp, pad);
}
