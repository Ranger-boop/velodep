#include "frac_aperture.h"

FracAperture::FracAperture(
	double dil_fact, double D_c, double dil_ang, double b_0,
	const double u_end[], size_t num_u_end,
	const double u_ini[], size_t num_u_ini,
	const double v[], size_t num_v
) {
	this->p_impl = new FracApertureImpl(
		dil_fact, D_c, dil_ang, b_0,
		u_end, num_u_end,
		u_ini, num_u_ini,
		v, num_v
	);
}

FracAperture::FracAperture(const FracAperture& frac_apt)
{
	this->p_impl = new FracApertureImpl(*frac_apt.p_impl);
}

FracAperture& FracAperture::operator=(const FracAperture& frac_apt)
{
	if (this != &frac_apt) {
		delete this->p_impl;

		this->p_impl = new FracApertureImpl(*frac_apt.p_impl);
	}

	return *this;
}

void FracAperture::solve()
{
	this->p_impl->solve();
}

size_t FracAperture::get_len_b()
{
	return this->p_impl->get_len_b();
}

double* FracAperture::get_b_mod()
{
	return this->p_impl->get_b_mod();
}

double* FracAperture::get_b_slip()
{
	return this->p_impl->get_b_slip();
}
