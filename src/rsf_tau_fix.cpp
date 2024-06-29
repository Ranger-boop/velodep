#include <cmath>
#include <iostream>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include "rsf_tau_fix.h"

void RSF_tau_fix::set_const_paras(double a, double b, double Dc, double alpha)
{
    this->p_impl->set_const_paras(a, b, Dc, alpha);
}

void RSF_tau_fix::set_ini_vals(double mu_0, double v_0, double disp_0)
{
    this->p_impl->set_ini_vals(mu_0, v_0, disp_0);
}

void RSF_tau_fix::set_indpt_val(
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

void RSF_tau_fix::compute_sig_neff()
{
    this->p_impl->compute_sig_neff();
}

void RSF_tau_fix::compute_dSig_neff_dt()
{
    this->p_impl->compute_dSig_neff_dt();
}

RSF_tau_fix::RSF_tau_fix()
{
    this->p_impl = new RsfTauFixImpl();
}

RSF_tau_fix::RSF_tau_fix(double a, double b, double Dc, double alpha)
{
    this->p_impl = new RsfTauFixImpl(a, b, Dc, alpha);
}

RSF_tau_fix::RSF_tau_fix(
    double a, double b, double Dc, double alpha, 
    double mu_0, double v_0, double disp_0
)
{
    this->p_impl = new RsfTauFixImpl(
        a, b, Dc, alpha,
        mu_0, v_0, disp_0
    );
}

RSF_tau_fix::RSF_tau_fix(
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

RSF_tau_fix::RSF_tau_fix(const RSF_tau_fix& rsf)
{
    this->p_impl = new RsfTauFixImpl(*rsf.p_impl);
}

RSF_tau_fix& RSF_tau_fix::operator=(const RSF_tau_fix& rsf)
{
    if (this != &rsf) {
        delete this->p_impl;
        this->p_impl = new RsfTauFixImpl(*rsf.p_impl);
    }

    return *this;
}

void RSF_tau_fix::solve_rk4()
{
    this->p_impl->solve_rk4();
}

size_t RSF_tau_fix::get_len_out_para() const
{
    return this->p_impl->get_len_out_para();
}

double* RSF_tau_fix::get_disp()
{
    return this->p_impl->get_disp();
}

double* RSF_tau_fix::get_theta()
{
    return this->p_impl->get_theta();
}

double* RSF_tau_fix::get_tau()
{
    return this->p_impl->get_tau();
}

double* RSF_tau_fix::get_vel()
{
    return this->p_impl->get_vel();
}

double* RSF_tau_fix::get_mu()
{
    return this->p_impl->get_mu();
}

RSF_tau_fix::~RSF_tau_fix()
{
    if (this->p_impl) {
        delete this->p_impl;
        this->p_impl = nullptr;
    }
}
