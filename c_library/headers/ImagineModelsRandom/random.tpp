

template<typename POSTYPE, typename GRIDTYPE>
double RandomField<POSTYPE, GRIDTYPE>::hammurabi_spectrum(const double &abs_k, const double &rms, const double &k0, const double &k1, const double &a0, const double &a1) const {
  // this function is adapted from https://github.com/hammurabi-dev/hammurabiX/blob/master/source/field/b/brnd_jf12.cc
  // original author: https://github.com/gioacchinowang
  const double p0 = rms*rms;
  double pi = 3.141592653589793;
  const double unit = 1. / (4 * pi * abs_k * abs_k);   // units fixing, wave vector in 1/kpc units
  // power laws
  const double band1 = double(abs_k < k1);
  const double band2 = double(abs_k > k1) * double(abs_k < k0);
  const double band3 = double(abs_k > k0);
  const double P = band1 * std::pow(k0 / k1, a1) * std::pow(abs_k / k1, 6.0) +
                  band2 / std::pow(abs_k / k0, a1) +
                  band3 / std::pow(abs_k / k0, a0);
  return P * p0 * unit;
  }

template<typename POSTYPE, typename GRIDTYPE>
double RandomField<POSTYPE, GRIDTYPE>::simple_spectrum(const double &abs_k, const double &dk, const double &k0, const double &s) const {
  double pi = 3.141592653589793;
  const double unit = 1. / (4 * pi * abs_k * abs_k); 
  const double dP = unit / std::pow(abs_k + k0, s);
  //double norm = 1. / ((s - 1.) * std::pow(k0, (s - 1))); // normalize to unity
  return dP; // / norm;
}

template<typename POSTYPE, typename GRIDTYPE>
void RandomField<POSTYPE, GRIDTYPE>::remove_padding(double* val, const std::array<int, 3> &shp, const int pad) {
  int start = 0;
  int sz = shp[2];
  int n = shp[0]*shp[1];
  for (int i = 1; i<n; i++) {
      std::copy(val + i*(sz + pad), val + i*(sz + pad) + sz, val + i*sz);
  }
}


template<typename POSTYPE, typename GRIDTYPE>
void RandomField<POSTYPE, GRIDTYPE>::seed_complex_random_numbers(fftw_complex* vec,  const std::array<int, 3> &shp, const std::array<double, 3> &inc, const int seed)  {

  bool debug_random = false;
  auto gen = std::mt19937(seed);
  int no_of_real = 0;
  int no_of_free = 0;
  int no_of_rand = 0;

  double lx = shp[0]*inc[0];
  double ly = shp[1]*inc[1];
  double lz = shp[2]*inc[2];
  
  float nyquist_x = shp[0]/2.;
  float nyquist_y = shp[1]/2.;
  float nyquist_z = shp[2]/2.;

  int size_z = static_cast<int>(nyquist_z) + 1;

  for (int i = 0; i < shp[0]; ++i) {
    double kx = (double)i / lx;
    if (i > nyquist_x)
      kx -= 1./ inc[0];
    const int idx_lv1 = i * shp[1] * size_z;
    for (int j = 0; j < shp[1]; ++j) {
      double ky = (double)j / ly;
      if (j > nyquist_y)
        ky -= 1./ inc[1];
      const int idx_lv2 = idx_lv1 + j * size_z;
      for (int l = 0; l < size_z; ++l) {
        const int idx = idx_lv2 + l;
        double kz = (double)l / lz;
        if (debug_random) {
          std::cout << "At Index (i, j, k): (" << i << j << l << ")" << std::endl;
          std::cout << "flat array index " << idx <<  std::endl;
        }
        if (l == 0 and j == 0 and i == 0) {
          // Full Monopole is set to zero, dealt with seperately in the outer scope
          vec[0][0] = 0.;
          vec[0][1] = 0.;
          if (debug_random) {
          std::cout << "Type: Monopole" << std::endl;
          std::cout << "Array val (real, imag): " << vec[idx][0] << ", " << vec[idx][1] << std::endl;
          std::cout <<  "\n";
          }
          continue;
        }
        
        double sigma;

        if (apply_spectrum) {
          const double ks = std::sqrt(kx * kx + ky * ky + kz * kz);
          sigma = calculate_fourier_sigma(ks, 1./(lx*ly*lz));
        }
        else {
          sigma = 1.;
        }
        std::normal_distribution<double> nd{0, sigma};

        bool l_is_zero_or_nyquist = (l == 0 or l == nyquist_z);
        bool j_is_zero_or_nyquist = (j == 0 or j == nyquist_y);
        bool i_is_zero_or_nyquist = (i == 0 or i == nyquist_x);
        int cg_idx;


        if (l_is_zero_or_nyquist) { // real z planes
          if (j_is_zero_or_nyquist) { // real y_lines 
            if (i_is_zero_or_nyquist) {  // line monopole or nyquist, ->draw real numbers (global momopole is dealt with earlier)
              vec[idx][0] = nd(gen);
              vec[idx][1] = 0.;     
              if (debug_random) {
                no_of_rand += 1;
                no_of_real += 1;
              }  
            }
            else if (i < nyquist_x) {  // line values below nyqist
              vec[idx][0] = nd(gen);
              vec[idx][1] = nd(gen);
              if (debug_random) {
                no_of_rand += 2;
                no_of_free += 1;
              }
            }
            else {                     // complex conjugate online
            cg_idx = (shp[0] - i)  * shp[1] * size_z + j* size_z  + l;
              vec[idx][0] = vec[cg_idx][0];
              vec[idx][1] = - vec[cg_idx][1];
            }  
          }
          else {                 // complex lines below nyqist  
            if (i_is_zero_or_nyquist) {
              if (j  < nyquist_y) {
                vec[idx][0] = nd(gen);
                vec[idx][1] = nd(gen);
                if (debug_random) {
                  no_of_rand += 2;
                  no_of_free += 1;
                }
              }
              else {
                cg_idx = i * shp[1] * size_z + (shp[1] - j) * size_z  + l;
                vec[idx][0] = vec[cg_idx][0];
                vec[idx][1] = - vec[cg_idx][1];
              }
            } 
            else if (i < nyquist_x) {
              vec[idx][0] = nd(gen);
              vec[idx][1] = nd(gen);
              if (debug_random) {
                no_of_rand += 2;
                no_of_free += 1;
              }
            }
            else {
              cg_idx = (shp[0] - i) * shp[1] * size_z + (shp[1] - j) * size_z  + l;
              vec[idx][0] = vec[cg_idx][0];
              vec[idx][1] = - vec[cg_idx][1];
            }
          }
        }
        else { //  complex z planes
          vec[idx][0] = nd(gen);
          vec[idx][1] = nd(gen);
          if (debug_random) {
            no_of_rand += 2;
            no_of_free += 1;
            }
        }
      }
    }
  }
  if (debug_random) {
    std::cout <<  "\n" << std::endl;

    std::cout <<  "\nnumber of real: " << no_of_real << "\nnumber of free: " << no_of_free<<  std::endl;    std::cout <<  "\ndegrees of freedom: " << (shp[1] * nyquist_z*2 * shp[0]) << "\nnumber of random numbers drawn: " << no_of_rand << std::endl;
    std::cout <<  "\n" << std::endl;
  }

}
