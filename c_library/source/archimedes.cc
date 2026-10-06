#include "ImagineModels/Archimedes.h"

namespace imagine {

// J. L. Han et al 2018 ApJS 234 11
template <typename T>
Vec3<T> ArchimedeanMagneticField::field(const double &x, const double &y, const double &z, const ArchimedeanParameters<T> &p) const
{

    Vec3<T> B_cart{{0., 0., 0.}};
    const double r = sqrt(x * x + y * y + z * z);
    if (r == 0.)
        return B_cart;

	double theta = atan2(sqrt(x * x + y * y), z);
	double phi = std::atan2(y, x); 
	
	double cos_phi = cos(phi);
	double sin_phi = sin(phi);
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);

	// radial direction
	auto c1 = p.R_0*p.R_0/r/r;
	B_cart[0] += c1 * cos_phi * sin_theta;
	B_cart[1] += c1 * sin_phi * sin_theta;
	B_cart[2] += c1 * cos_theta;
	
	// azimuthal direction	
	auto c2 = - (p.Omega*p.R_0*p.R_0*sin_theta) / (r*p.v_w);
	B_cart[0] += c2 * (-sin_phi);
	B_cart[1] += c2 * cos_phi;

	// magnetic field switch at z = 0
	auto B_0 = p.B_0;

	if (z<0.) {
		B_0 *= -1;
	}

	// overall scaling
	B_cart[0] *= B_0;
	B_cart[1] *= B_0;
	B_cart[2] *= B_0;

	return B_cart;

}


IMAGINE_INSTANTIATE_VECTOR_MODEL(ArchimedeanMagneticField)

}
