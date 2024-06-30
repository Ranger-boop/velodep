#include <cmath>
#include <iostream>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include "rsf_tau_fix.h"

void RsfTauFix::set_const_paras(double a, double b, double Dc, double alpha)
{
    this->p_impl->set_const_paras(a, b, Dc, alpha);
}

void RsfTauFix::set_ini_vals(double mu_0, double v_0, double disp_0)
{
    this->p_impl->set_ini_vals(mu_0, v_0, disp_0);
}

void RsfTauFix::set_indpt_val(
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    this->p_impl->set_indpt_val(
        sig_n, num_sig_n,
        pp, num_pp,
        t, num_t
    );
}

void RsfTauFix::compute_sig_neff()
{
    this->p_impl->compute_sig_neff();
}

void RsfTauFix::compute_dSig_neff_dt()
{
    this->p_impl->compute_dSig_neff_dt();
}

RsfTauFix::RsfTauFix()
{
    this->p_impl = new RsfTauFixImpl();
}

RsfTauFix::RsfTauFix(double a, double b, double Dc, double alpha)
{
    this->p_impl = new RsfTauFixImpl(a, b, Dc, alpha);
}

RsfTauFix::RsfTauFix(
    double a, double b, double Dc, double alpha, 
    double mu_0, double v_0, double disp_0
)
{
    this->p_impl = new RsfTauFixImpl(
        a, b, Dc, alpha,
        mu_0, v_0, disp_0
    );
}

RsfTauFix::RsfTauFix(
    double a, double b, double Dc, double alpha, 
    double mu_0, double v_0, double disp_0,
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    this->p_impl = new RsfTauFixImpl(
        a, b, Dc, alpha,
        mu_0, v_0, disp_0,
        sig_n, num_sig_n,
        pp, num_pp,
        t, num_t
    );
}

RsfTauFix::RsfTauFix(const RsfTauFix& rsf)
{
    this->p_impl = new RsfTauFixImpl(*rsf.p_impl);
}

RsfTauFix& RsfTauFix::operator=(const RsfTauFix& rsf)
{
    if (this != &rsf) {
        delete this->p_impl;
        this->p_impl = new RsfTauFixImpl(*rsf.p_impl);
    }

    return *this;
}

void RsfTauFix::solve_rk4()
{
    this->p_impl->solve_rk4();
}

size_t RsfTauFix::get_len_out_para() const
{
    return this->p_impl->get_len_out_para();
}

double* RsfTauFix::get_disp()
{
    return this->p_impl->get_disp();
}

double* RsfTauFix::get_theta()
{
    return this->p_impl->get_theta();
}

double* RsfTauFix::get_tau()
{
    return this->p_impl->get_tau();
}

double* RsfTauFix::get_vel()
{
    return this->p_impl->get_vel();
}

double* RsfTauFix::get_mu()
{
    return this->p_impl->get_mu();
}

RsfTauFix::~RsfTauFix()
{
    if (this->p_impl) {
        delete this->p_impl;
        this->p_impl = nullptr;
    }
}
