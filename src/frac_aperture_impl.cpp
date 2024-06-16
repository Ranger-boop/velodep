#include <cmath>
#include <numbers>
#include <stdexcept>
#include "frac_aperture_impl.h"

FracApertureImpl::FracApertureImpl(
	double dil_fact, double D_c, double dil_ang, double b_0,
	const double u_end[], size_t num_u_end,
	const double u_ini[], size_t num_u_ini,
	const double v[], size_t num_v
) {
	this->dil_fact = dil_fact;
	this->D_c = D_c;
	this->dil_ang = dil_ang;
	this->b_0 = b_0;

	// Set independent values
	if ((num_u_end != num_u_ini) || (num_u_ini != num_v) || (num_v != num_u_end)) {
		throw std::logic_error("The lengths of u_end, u_ini, and v must be the same.");
	}
	else {
		this->u_end = std::valarray<double>(u_end, num_u_end);
		this->u_ini = std::valarray<double>(u_ini, num_u_ini);
		this->v = std::valarray<double>(v, num_v);
	}

	//this->b_slip = std::valarray<double>(0., num_u_end);
	//this->b_mod = std::valarray<double>(0., num_u_end);
	// Allocate memory for dilation parameters
	//this->d_phi_2dim.resize(num_v);
}

//void FracApertureImpl::aperture_slip_disp()
//{
//	std::valarray<double> delta_u = this->u_end - this->u_ini[0];
//	this->b_slip = this->b_0 + delta_u * std::tan(this->dil_ang * std::numbers::pi / 180);
//}

std::valarray<double> FracApertureImpl::dil_para(
	double un_end,
	const std::valarray<double> &ui_ini,
	const std::valarray<double> &v
) {
	std::valarray<double> d_phi = std::valarray<double>(0., v.size());
	//d_phi[0] = 0;
	if (v.size() > 1) {
		std::valarray<double> v_b = v[std::slice(0, v.size() - 1, 1)];
		std::valarray<double> v_f = v[std::slice(1, v.size() - 1, 1)];

		//std::valarray<double> dt = std::valarray<double>(0., v.size());
		std::valarray<double> dt = (un_end - ui_ini) / v;
		//dt[0] = 0;

		std::slice l_slice(1, d_phi.size() - 1, 1);
		std::valarray<double> dt_slice(dt[l_slice]);
		d_phi[l_slice] = -this->dil_fact * std::log(
			v_b / v_f * (1 + (v_f / v_b - 1) * std::exp(-v_f * dt_slice / this->D_c))
		);

		/*for (size_t i = 1; i < d_phi.size(); i++) {
			d_phi[i] = -dil_fact * std::log(
				v_b[i - 1] / v_f[i - 1] * 
				(1 + (v_f[i - 1] / v_b[i - 1] - 1) * std::exp(-v_f[i - 1] * dt[i] / D_c))
				);
		}*/
	}
	return d_phi;
}

double FracApertureImpl::b_mod_n(
	double b_slip_n, 
	const std::valarray<double> &d_phi_1dim
) {
	double ratio_bslip_cumprod = 1;
	for (const double& d_phi : d_phi_1dim) {
		ratio_bslip_cumprod *= 1 + d_phi;
	}
	return b_slip_n * ratio_bslip_cumprod;
}

void FracApertureImpl::solve()
{
	//this->aperture_slip_disp();  // Compute b_slip
	std::valarray<double> delta_u = this->u_end - this->u_ini[0];
	this->b_slip = this->b_0 + delta_u * std::tan(this->dil_ang * std::numbers::pi / 180);

	this->d_phi_2dim.resize(this->b_slip.size());
	for (size_t i = 0; i < this->b_slip.size(); i++) {
		std::slice l_slice(0, i + 1, 1);
		this->d_phi_2dim[i] = dil_para(this->u_end[i], this->u_ini[l_slice], this->v[l_slice]);
	}

	this->b_mod = std::valarray<double>(0., this->b_slip.size());
	for (size_t i = 0; i < this->b_slip.size(); i++) {
		this->b_mod[i] = b_mod_n(this->b_slip[i], this->d_phi_2dim[i]);
	}
}

size_t FracApertureImpl::get_len_b()
{
	return this->b_mod.size();
}

double* FracApertureImpl::get_b_mod()
{
	return &(this->b_mod[0]);
}