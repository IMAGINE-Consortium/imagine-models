#include <cmath>
#include <iostream>

#include "RandomVectorField.h"

// Non trivial constructor
RandomVectorField::RandomVectorField(std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment) : RandomField(shape, reference_point, increment) {
    int newshp2;
    if (shape[2] % 2) {
      newshp2 = shape[2] + 1;
    }
    else {
      newshp2 = shape[2] + 2;
    }
    double* val_temp = (double*) fftw_alloc_real(shape[0]*shape[1]*newshp2);
    fftw_complex* val_temp_comp = reinterpret_cast<fftw_complex*>(val_temp);
    fftw_plan r2c_temp = fftw_plan_dft_r2c_3d(shape[0], shape[1], shape[2], val_temp, val_temp_comp, FFTW_MEASURE);
    fftw_plan c2r_temp = fftw_plan_dft_c2r_3d(shape[0], shape[1], shape[2], val_temp_comp, val_temp,  FFTW_MEASURE);
    const char *filename = "ImagineModelsRandomVectorField";
    int fftw_export_wisdom_to_filename(*filename);
    has_fftw_wisdom = true;
    fftw_free(val_temp);
    fftw_destroy_plan(c2r_temp);
    fftw_destroy_plan(r2c_temp); // plans are destroyed but wisdom is saved!
}

RandomVectorField::~RandomVectorField() {
  destroy_plans();
  fftw_forget_wisdom();
  #ifdef _OPENMP
    fftw_cleanup_threads();
  #else
    fftw_cleanup();
  #endif
};

std::array<double*, 3> RandomVectorField::allocate_memory(std::array<int, 3> shp) {
    std::array<double*, 3> grid_eval;
    int newshp2;
    if (shp[2] % 2) {
      newshp2 = shp[2] + 1;
    }
    else {
      newshp2 = shp[2] + 2;
    }
    for (int i=0; i < ndim; ++i) {
      grid_eval[i] = (double*) fftw_alloc_real(shp[0]*shp[1]*newshp2);
      }
    return grid_eval;  
}

void RandomVectorField::free_memory(std::array<double*, 3> grid_eval) {
    for (int i=0; i < ndim; ++i) { 
      fftw_free(grid_eval[i]);
    }
}


std::array<fftw_complex*, 3> RandomVectorField::construct_plans(std::array<double*, 3> grid_eval, std::array<int, 3> shp) {
    std::array<fftw_complex*, 3> grid_eval_comp;
      for (int i=0; i < ndim; ++i) {
        grid_eval_comp[i] = reinterpret_cast<fftw_complex*>(grid_eval[i]);
        r2c[i] = fftw_plan_dft_r2c_3d(shp[0], shp[1], shp[2], grid_eval[i], grid_eval_comp[i], FFTW_ESTIMATE);
        c2r[i] = fftw_plan_dft_c2r_3d(shp[0], shp[1], shp[2], grid_eval_comp[i], grid_eval[i],  FFTW_ESTIMATE);
      }
      created_fftw_plans = true;
    return grid_eval_comp;
}

void RandomVectorField::destroy_plans() {
    if (created_fftw_plans) {
      for (int i=0; i < ndim; ++i) {
        fftw_destroy_plan(c2r[i]);
        fftw_destroy_plan(r2c[i]);
      }
    }
}

std::array<double*, 3> RandomVectorField::on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {
    std::array<double*, 3> grid_eval = allocate_memory(shp);
    _on_grid(grid_eval, shp, rpt, inc, seed);
    return grid_eval;
  }

std::array<double*, 3> RandomVectorField::on_grid(const int seed) {
  if (not initialized_with_grid) 
    throw GridException();
  std::array<double*, 3> grid_eval = allocate_memory(internal_shape);
  const char *filename = "ImagineModelsRandomVectorField";
  int fftw_import_wisdom_from_filename(*filename);
  _on_grid(grid_eval, internal_shape, internal_ref_point, internal_increment, seed);
  return grid_eval;
}

double* RandomVectorField::profile_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &rfp, const std::array<double, 3> &inc) {
  double* grid_eval;
  size_t arr_sz = shp[0]*shp[1]*shp[2];
  grid_eval = new double[arr_sz];
  evaluate_function_on_grid<number, double*>(grid_eval, shp, rfp, inc,
                                    [this](double xx, double yy, double zz)
                                    { return spatial_profile(xx, yy, zz); });
  return grid_eval;
} 
std::array<double*, 3> RandomVectorField::random_numbers_on_grid(const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed) {
    std::array<double*, 3> val = allocate_memory(shp);
    std::array<fftw_complex*, 3> val_comp = construct_plans(val, shp); 
    int gs = grid_size(shp);
    double sqrt_gs = std::sqrt(gs);
    auto gen_int = std::mt19937(seed);
    std::uniform_int_distribution<int> uni(0, 1215752192);
    std::array<int, 3> padded_shp = {shp[0],  shp[1],  2*(shp[2]/2 + 1)}; 
    int padded_size = grid_size(padded_shp);
    int pad =  padded_shp[2] - shp[2];

    for (int i =0; i<3; ++i) {
      int sub_seed = uni(gen_int); 
      seed_complex_random_numbers(val_comp[i], shp, inc, sub_seed);
      fftw_execute(c2r[i]);
      for (int s = 0; s < padded_size; ++s)  {
        (val[i])[s] /= sqrt_gs;  
      }
      remove_padding(val[i], shp, pad);
    }
    return val;
}

void RandomVectorField::_on_grid(std::array<double*, 3> val, const std::array<int, 3> &shp, const std::array<double, 3> &rpt, const std::array<double, 3> &inc, const int seed) {

  std::array<fftw_complex*, 3> val_comp = construct_plans(val, shp); 
  std::array<int, 3> padded_shp = {shp[0],  shp[1],  2*(shp[2]/2 + 1)}; 
  int gs = grid_size(shp);   
  int padded_size = grid_size(padded_shp);
  int pad =  padded_shp[2] - shp[2];
  auto gen_int = std::mt19937(seed);
  std::uniform_int_distribution<int> uni(0, 1215752192);
  double sqrt_gs = std::sqrt(gs);

  // Step 1: draw random numbers with variance 1, possibly correlated

  for (int i =0; i<3; ++i) {
    int sub_seed = uni(gen_int); 
    seed_complex_random_numbers(val_comp[i], shp, inc, sub_seed);
    fftw_execute(c2r[i]);
    //for (int s = 0; s < padded_size; ++s)  {
    //    (val[i])[s] /= sqrt_gs;  
    //  }
    //remove_padding(val[i], shp, pad);
  }
  
  // Step 2: apply spatial amplitude, possibly introduce anisotropy depending on regular field.
  if (!no_profile) {
    auto apply_profile = [&](std::array<double, 3> &b_rand_val, const double xx, const double yy, const double zz) {
        
      double sp = spatial_profile(xx, yy, zz);
      // apply profile
      b_rand_val[0] *= sp;
      b_rand_val[1] *= sp;
      b_rand_val[2] *= sp;


      if (apply_anisotropy) {
        vector b_reg_val = anisotropy_direction(xx, yy, zz);
        double b_reg_x = static_cast<double>(b_reg_val[0]); 
        double b_reg_y = static_cast<double>(b_reg_val[1]);
        double b_reg_z = static_cast<double>(b_reg_val[2]);

        double b_reg_length = std::sqrt(std::pow(b_reg_x, 2) + std::pow(b_reg_y, 2) + std::pow(b_reg_z, 2));

        if (b_reg_length > 1e-10) { // non zero regular field, -> prefered anisotropy
          
          b_reg_x /= b_reg_length;
          b_reg_y /= b_reg_length;
          b_reg_z /= b_reg_length;
          const double rho2 = anisotropy_rho * anisotropy_rho;
          const double rhonorm = 1. / std::sqrt(0.33333333 * rho2 + 0.66666667 / rho2);
          double reg_dot_rand  = b_reg_x*b_rand_val[0] + b_reg_y*b_rand_val[1] + b_reg_z*b_rand_val[2];

          for (int ii=0; ii==3; ++ii) {
            double b_rand_par = b_rand_val[ii] / reg_dot_rand;
            double b_rand_perp = b_rand_val[ii]  - b_rand_par;
            b_rand_val[ii] = (b_rand_par * anisotropy_rho + b_rand_perp / anisotropy_rho) * rhonorm;
          } 
        }
      }
      return b_rand_val;
    };
    apply_function_to_field<std::array<double*, 3>, std::array<double, 3>>(val, padded_shp, rpt, inc, apply_profile);
  }
  // Step 3 (optional): divergence cleaning using Gram Schmidt process
  if (clean_divergence) {
  
    for (int i =0; i<3; ++i) {
      fftw_execute(r2c[i]);
    }
    divergence_cleaner(val_comp[0], val_comp[1], val_comp[2], shp, inc);
    
    for (int i =0; i<3; ++i) {
      fftw_execute(c2r[i]);
      for (int s = 0; s < padded_size; ++s)  {
        (val[i])[s] /= (gs*sqrt_gs);  
      }
      remove_padding(val[i], shp, pad);
    }
  }
  else {
    for (int i =0; i<3; ++i) {
      double sqrt_gs = std::sqrt(gs);
      for (int s = 0; s < padded_size; ++s)  {
        (val[i])[s] /= sqrt_gs;  
      }
      remove_padding(val[i], shp, pad);
    }
  }
}

// this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
// original author: https://github.com/gioacchinowang
void RandomVectorField::divergence_cleaner(fftw_complex* bx, fftw_complex* by, fftw_complex* bz,  const std::array<int, 3> &shp, const std::array<double, 3> &inc) const {
    double lx = shp[0]*inc[0];
    double ly = shp[1]*inc[1];
    double lz = shp[2]*inc[2];
  
    #ifdef _OPENMP
      #pragma omp parallel for schedule(static)
    #endif
      for (int i = 0; i < shp[0]; ++i) {
        double kx = i / lx;
        if (i >= (shp[0] + 1) / 2)
          kx -= 1. /  inc[0];
          // it's faster to calculate indices manually
        const int idx_lv1 = i * shp[1] * shp[2];
        for (int j = 0; j < shp[1]; ++j) {
          double ky = j /  ly;
          if (j >= (shp[1] + 1) / 2)
            ky -= 1. /  inc[1];
          const int idx_lv2 = idx_lv1 + j * shp[2];
          for (int l = 0; l < (int)shp[2]/2 + 1; ++l) {
            // 0th term is fixed to zero in allocation
            if (i == 0 and j == 0 and l == 0)
              continue;
            double kz = l /  lz;
            const int idx = idx_lv2 + l;
            double k_length = 0;
            double b_length = 0;
            double b_dot_k = 0;
            std::array<double, 3> k{kx, ky, kz};
            std::array<double, 3> b{(*bx)[idx], (*by)[idx], (*bz)[idx]};
            b_length = static_cast<double>(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
            k_length = static_cast<double>(k[0]*k[0] + k[1]*k[1] + k[2]*k[2]);

            if (k_length == 0 or b_length == 0) {
              continue;
              }
            k_length = std::sqrt(k_length);
            b_dot_k = (b[0]*k[0] + b[1]*k[1] + b[2]*k[2]);

            const double bk_over_k = b_dot_k / k_length;
            // multiply \sqrt(3) for preserving spectral power statistically
            (*bx)[idx] = 1.73205081 * ((*bx)[idx] - k[0] * bk_over_k);
            (*by)[idx] = 1.73205081 * ((*by)[idx] - k[1] * bk_over_k);
            (*bz)[idx] = 1.73205081 * ((*bz)[idx] - k[2] * bk_over_k);
          } // l
        } // j
      } // i
    }