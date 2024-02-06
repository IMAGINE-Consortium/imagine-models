#include <cmath>
#include "hamunits.h"
#include "Pshirkov.h"

vector PshirkovMagneticField::_at_position(const double &x, const double &y, const double &z, const PshirkovMagneticField &p) const
{
	const double r = std::sqrt(x * x + y * y); // radius in cylindrical coordinates
	const double phi = atan2(y, x);
	vector b{{0.0, 0.0, 0.0}};

	if ((x == 0.) && (y == 0.)) {
		return b;
	}

	auto pitch = p.pitch * M_PI / 180;

	auto cos_pitch = cos(pitch);
	auto sin_pitch = sin(pitch);
	auto PHI = cos_pitch / sin_pitch * log(1. + p.d / p.R_sun) - M_PI / 2;
	auto cos_PHI = cos(PHI);

	// disk field
	if ((p.useASS) or (p.useBSS)) {
	// CRPROPA COMMENT:
		// PT11 paper has B_theta = B * cos(p) but this seems because they define azimuth clockwise, while we have anticlockwise.
		// see Tinyakov 2002 APh 18,165: "local field points to l=90+p" so p=-5 deg gives l=85 and hence clockwise from above.
		// so to get local B clockwise in our system, need minus (like Sun etal).
		// Ps base their system on Han and Qiao 1994 A&A 288,759 which has a diagram with azimuth clockwise, hence confirmed.
		// PT11 paper define Earth position at (+8.5, 0, 0) kpc; but usual convention is (-8.5, 0, 0)
		// thus we have to rotate our position by 180 degree in azimuth
		auto theta = M_PI - phi;  // // CRPROPA COMMENT: azimuth angle theta: PT11 paper uses opposite convention for azimuth
		// // CRPROPA COMMENT: the following is equivalent to sin(pi - phi) and cos(pi - phi) which is computationally slower
		double cos_theta = - x / r;
		double sin_theta = y / r;

	// CRPROPA COMMENT:
		// After some geometry calculations (on whiteboard) one finds:
		// Bx = +cos(theta) * B_r - sin(theta) * B_{theta}
		// By = -sin(theta) * B_r - cos(theta) * B_{theta}
		// Use from paper: B_theta = B * cos(pitch)	and B_r = B * sin(pitch)
		b[0] = - sin_pitch * cos_theta + cos_pitch * sin_theta;
		b[1] = sin_pitch * sin_theta + cos_pitch * cos_theta;
	// ADAPTED CRPROPA COMMENT: flipped in eq above magnetic field direction, as B_{theta} and B_{phi} refering to 180 degree rotated field

		auto bMag = cos(theta - cos_pitch / sin_pitch * log(r / p.R_sun) + PHI);  // eq. 3 / 4
		if ((p.useASS) and (bMag < 0))
			bMag *= -1.;
		bMag *= p.B0_D * p.R_sun / std::max(r, p.R_c) / cos_PHI * exp(-fabs(z) / p.z0_D);  // eq. 5, eq. 4
		b[0] *= bMag;
		b[1] *= bMag;
		b[2] *= bMag; // does not do anything as b[2] was zero
		}

	// halo field
	if (p.useHalo) {
		auto bMag = (z > 0 ? p.B0_Hn : - p.B0_Hs);
		auto z1 = (fabs(z) < p.z0_H ? p.z11_H : p.z12_H);
		bMag *= r / p.R0_H * exp(1 - r / p.R0_H) / (1 + pow((fabs(z) - p.z0_H) / z1, 2.));
	// CRPROPA COMMENT:
		// equation (8) in paper: theta uses now the conventional azimuth definition in contrast to equation (3)
		// cos(phi) = pos.x / r (phi going counter-clockwise)
		// sin(phi) = pos.y / r
		// unitvector of phi in polar coordinates: (-sin(phi), cos(phi), 0)
		b[0] += -y / r * bMag;
		b[1] += x / r * bMag;
	}

	return b;
}

#if autodiff_FOUND

Eigen::MatrixXd PshirkovMagneticField::_jac(const double &x, const double &y, const double &z, PshirkovMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, PshirkovMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.pitch, p.d, p.R_sun, p.z0_D, p.B0_D, p.z0_H, p.R0_H, p.B0_Hn, p.B0_Hs, p.z11_H, p.z12_H), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif
