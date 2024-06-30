#ifndef C_API_H_
#define C_API_H_

#include "export_import.h"
#include "rsf_tau_fix.h"
#include "rsf_tau_chg.h"
#include "frac_aperture.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creat a RsfTauFix object
API RsfTauFix* new_RsfTauFix(
    double a, double b, double Dc, double alpha, 
    double mu_0, double v_0, double disp_0,
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
);

// Creat a RsfTauChg object
API RsfTauChg* new_RsfTauChg(
    double a, double b, double Dc, double alpha, 
    double k_s, double mu_0, double v_0, double disp_0,
    const double vlp[], size_t num_vlp,
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
);

// Destruct the pointed object
API void delete_RSF(RsfTauFix* p_rsf);

// Call the solver to solve the RSF equation group
API void solve_RSF(RsfTauFix* p_rsf);

// Get the length of the solution array of the RSF equation group
API size_t get_len_RSF_sol(RsfTauFix* p_rsf);

API double* get_disp_RSF(RsfTauFix* p_rsf);

API double* get_theta_RSF(RsfTauFix* p_rsf);

API double* get_tau_RSF(RsfTauFix* p_rsf);

API double* get_vel_RSF(RsfTauFix* p_rsf);

API double* get_mu_RSF(RsfTauFix* p_rsf);


// C API for the class FracAperture
// Creat a FracAperture object
API FracAperture* new_FracAperture(
    double dil_fact, double D_c, double dil_ang, double b_0,
    const double u_end[], size_t num_u_end,
    const double u_ini[], size_t num_u_ini,
    const double v[], size_t num_v
);

// Destruct the pointed object
API void delete_FracAperture(FracAperture* p_frac_apt);

// Solve the evolving fracture aperture
API void solve_FracAperture(FracAperture* p_frac_apt);

API size_t get_len_b_FracAperture(FracAperture* p_frac_apt);

API double* get_b_mod_FracAperture(FracAperture* p_frac_apt);

API double* get_b_slip_FracAperture(FracAperture* p_frac_apt);

#ifdef __cplusplus
}
#endif

#endif  // C_API_H_
