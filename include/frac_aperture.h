#ifndef FRAC_APERTURE_H_
#define FRAC_APERTURE_H_

#include "export_import.h"
#include "frac_aperture_impl.h"

class API FracAperture
{
public:
	FracAperture(
		double dil_fact, double D_c, double dil_ang, double b_0,
		const double u_end[], size_t num_u_end,
		const double u_ini[], size_t num_u_ini,
		const double v[], size_t num_v
	);

	~FracAperture();

	FracAperture(const FracAperture& frac_apt);

	FracAperture& operator=(const FracAperture& frac_apt);

	void solve();

	size_t get_len_b();

	double* get_b_mod();

	double* get_b_slip();

private:
	FracApertureImpl* p_impl;
};

#endif  // FRAC_APERTURE_H_