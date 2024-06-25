#ifndef FRAC_APERTURE_IMPL_H_
#define FRAC_APERTURE_IMPL_H_

//#include <memory>
#include <vector>
#include <valarray>

//#ifndef M_PI
//constexpr double m_pi = 3.14159265358979323846;
//#endif

class FracApertureImpl
{
public:
	FracApertureImpl(
		double dil_fact, double D_c, double dil_ang, double b_0,
		const double u_end[], size_t num_u_end,
		const double u_ini[], size_t num_u_ini,
		const double v[], size_t num_v
	);

	~FracApertureImpl() = default;

	FracApertureImpl(const FracApertureImpl& frac_apt_impl) = default;

	FracApertureImpl & operator=(const FracApertureImpl& frac_apt_impl) = default;

	void solve();

	size_t get_len_b();

	double* get_b_mod();

	double* get_b_slip();

private:
	//void aperture_slip_disp();

	std::valarray<double> dil_para(
		double un_end, 
		const std::valarray<double> &ui_ini, 
		const std::valarray<double> &v
	);

	double b_mod_n(double b_slip_n, const std::valarray<double> &d_phi_1dim);

	//std::valarray<double> aperture_shear_dil();

	// Consitutive parameters:
	// dilation factor (dil_fact), characteristic slip distance (D_c), 
	// dilation angle (dil_ang)
	double dil_fact, D_c, dil_ang;

	// Initial aperture at the onset of reactivation (the beginning of the first 
	// velocity step)
	double b_0;

	// Slip displacement after the reactivation of fracture at the end of each 
	// velocity step
	std::valarray<double> u_end;
	// Slip displacement after the reactivation of fracture at the initiation 
	// of each velocity step
	std::valarray<double> u_ini;
	// Slip velocity at each velocity step
	std::valarray<double> v;

	// Fracture aperture dependent on slip displacement with the change of 
	// slip velocity not considered.
	std::valarray<double> b_slip;
	// Sequences of dilation parameter (i.e., incremental porosity) for each 
	// velocity step.
	// The lengths of these sequences should be different and increasing.
	std::vector<std::valarray<double>> d_phi_2dim;
	// Aperture influenced by both slip displacement and velocity of the fracture 
	// at each velocity step
	std::valarray<double> b_mod;
};

#endif  // FRAC_APERTURE_IMPL_H_
