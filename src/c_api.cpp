#include <stdexcept>
#include "c_api.h"

#ifdef __cplusplus
extern "C" {
#endif

RsfTauFix* new_RsfTauFix(
    double a, double b, double Dc, double alpha, 
    double mu_0, double v_0, double disp_0,
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    try {
        return new RsfTauFix(
            a, b, Dc, alpha,
            mu_0, v_0, disp_0,
            sig_n, num_sig_n,
            pp, num_pp,
            t, num_t
        );
    }
    catch (const std::logic_error& e) {
        std::cerr << "Logic error: " << e.what() << std::endl;
        std::abort();
    }
}

RsfTauChg* new_RsfTauChg(
    double a, double b, double Dc, double alpha, 
    double k_s, double mu_0, double v_0, double disp_0,
    const double vlp[], size_t num_vlp,
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    try {
        return new RsfTauChg(
            a, b, Dc, alpha,
            k_s, mu_0, v_0, disp_0,
            vlp, num_vlp,
            sig_n, num_sig_n,
            pp, num_pp,
            t, num_t
        );
    }
    catch (const std::logic_error& e) {
        std::cerr << "Logic error: " << e.what() << std::endl;
        std::abort();
    }
    
}

void delete_RSF(RsfTauFix* p_rsf)
{
    delete p_rsf;
    p_rsf = nullptr;
}

void solve_RSF(RsfTauFix* p_rsf)
{
    try {
        p_rsf->solve_rk4();
    }
    catch (const std::logic_error& e) {
        std::cerr << "Logic error: " << e.what() << std::endl;
        std::abort();
    }
    
}

size_t get_len_RSF_sol(RsfTauFix* p_rsf)
{
    return p_rsf->get_len_out_para();
}

double* get_disp_RSF(RsfTauFix* p_rsf)
{
    return p_rsf->get_disp();
}

double* get_theta_RSF(RsfTauFix* p_rsf)
{
    return p_rsf->get_theta();
}

double* get_tau_RSF(RsfTauFix* p_rsf)
{
    return p_rsf->get_tau();
}

double* get_vel_RSF(RsfTauFix* p_rsf)
{
    return p_rsf->get_vel();
}

double* get_mu_RSF(RsfTauFix* p_rsf)
{
    return p_rsf->get_mu();
}


// C API for the class FracAperture
FracAperture* new_FracAperture(
    double dil_fact, double D_c, double dil_ang, double b_0,
    const double u_end[], size_t num_u_end,
    const double u_ini[], size_t num_u_ini,
    const double v[], size_t num_v
) {
    try {
        return new FracAperture(
            dil_fact, D_c, dil_ang, b_0,
            u_end, num_u_end,
            u_ini, num_u_ini,
            v, num_v
        );
    }
    catch (const std::logic_error& e) {
        std::cerr << "Logic error: " << e.what() << std::endl;
        std::abort();
    }
}

void delete_FracAperture(FracAperture* p_frac_apt)
{
    delete p_frac_apt;
    p_frac_apt = nullptr;
}

void solve_FracAperture(FracAperture* p_frac_apt)
{
    p_frac_apt->solve();
}

size_t get_len_b_FracAperture(FracAperture* p_frac_apt)
{
    return p_frac_apt->get_len_b();
}

double* get_b_mod_FracAperture(FracAperture* p_frac_apt)
{
    return p_frac_apt->get_b_mod();
}

double* get_b_slip_FracAperture(FracAperture* p_frac_apt)
{
    return p_frac_apt->get_b_slip();
}

#ifdef __cplusplus
}
#endif
