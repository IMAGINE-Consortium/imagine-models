#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>
#include <array>


inline void Cyl2Cart( double phi, std::array<double, 3>  invec, std::array<double, 3>  outvec){
	double cosphi = std::cos(phi);
	double sinphi = std::sin(phi);
	outvec[0] = cosphi * invec[0] - sinphi * invec[1];
	outvec[1] = sinphi * invec[0] + cosphi * invec[1];
	outvec[2] = invec[2];
}// Cyl2Cart

#endif