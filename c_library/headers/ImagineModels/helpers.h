#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>
#include <array>

template<typename V>
V Cyl2Cart(const double phi, V& invec) {
	V outvec{{0., 0., 0.}};
	double cosphi = std::cos(phi);
	double sinphi = std::sin(phi);
	auto inv = cosphi * invec[0];
	outvec[0] = cosphi * invec[0] - sinphi * invec[1];
	outvec[1] = sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
	return outvec;
 }

 template<typename V>
V Cyl2Cart(V& invec, const double cosphi, const double sinphi) {
	V outvec{{0., 0., 0.}};
	auto inv = cosphi * invec[0];
	outvec[0] = cosphi * invec[0] - sinphi * invec[1];
	outvec[1] = sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
	return outvec;
 }

 template<typename V>
V Cart2Cyl(const double phi, V& invec) {
	V outvec{{0., 0., 0.}};
	double cosphi = std::cos(phi);
	double sinphi = std::sin(phi);
	auto inv = cosphi * invec[0];
	outvec[0] = cosphi * invec[0] + sinphi * invec[1];
	outvec[1] = - sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
	return outvec;
 }

 template<typename V>
V Cart2Cyl(V& invec, const double cosphi, const double sinphi) {
	V outvec{{0., 0., 0.}};
	auto inv = cosphi * invec[0];
	outvec[0] = cosphi * invec[0] + sinphi * invec[1];
	outvec[1] = - sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
	return outvec;
 }

template<typename V>
 vector addVector(std::initializer_list<vector> vs) {
    V outvec{{0., 0., 0.}}; 
	
	for (vector v : vs) {
		outvec[0] += v[0];
		outvec[1] += v[1];
		outvec[2] += v[2];
	}
	return outvec;
 }


template<typename s>
  inline
  s
  Sigmoid(s& x, const double x0, const double w)
  {
    return 1 / (1 + std::exp(-(x-x0)/w));
  }

  // angle between v0 = (cos(phi0), sin(phi0)) and v1 = (cos(phi1), sin(phi1))
template<typename out, typename in1, typename in2>
  inline
  out
  DeltaPhi(const in1 phi0, const in2 phi1)
  {
    return std::acos(std::cos(phi1)*std::cos(phi0) + std::sin(phi1)*std::sin(phi0));
  }


#endif

