#include "Archimedes.h"

// J. L. Han et al 2018 ApJS 234 11
vector ArchimedeanMagneticField::_at_position(const double &x, const double &y, const double &z, const ArchimedeanMagneticField &p) const
{

    vector B_cart{{0., 0., 0.}};
    const double r = sqrt(x * x + y * y + z * z);

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

#if autodiff_FOUND

Eigen::MatrixXd ArchimedeanMagneticField::_jac(const double &x, const double &y, const double &z, ArchimedeanMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, ArchimedeanMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.R_0, p.Omega, p.v_w, p.B_0), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif


