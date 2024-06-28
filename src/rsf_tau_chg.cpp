#include <cmath>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include <typeinfo>
#include "rsf_tau_chg.h"
#include "rsf_tau_chg_impl.h"

void RSF_tau_chg::set_ini_vals(double k_s, double mu_0, double v_0, double disp_0)
{
    if (RsfTauChgImpl* l_p_impl = dynamic_cast<RsfTauChgImpl*>(this->p_impl)) {
        l_p_impl->set_ini_vals(k_s, mu_0, v_0, disp_0);
    }
    else {
        std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
    }
}

void RSF_tau_chg::set_indpt_val(
    const double vlp[], size_t num_vlp,
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    if (RsfTauChgImpl* l_p_impl = dynamic_cast<RsfTauChgImpl*>(this->p_impl)) {
        l_p_impl->set_indpt_val(
            vlp, num_vlp,
            sig_n, num_sig_n,
            pp, num_pp,
            t, num_t
        );
    }
    else {
        std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
    }
}

RSF_tau_chg::RSF_tau_chg()
{
    this->p_impl = new RsfTauChgImpl();
}

RSF_tau_chg::RSF_tau_chg(
    double a, double b, double Dc, double alpha
)
{
    this->p_impl = new RsfTauChgImpl(a, b, Dc, alpha);
}

RSF_tau_chg::RSF_tau_chg(
    double a, double b, double Dc, double alpha, 
    double k_s, double mu_0, double v_0, double disp_0
)
{
    this->p_impl = new RsfTauChgImpl(
        a, b, Dc, alpha,
        k_s, mu_0, v_0, disp_0
    );
}

RSF_tau_chg::RSF_tau_chg(
    double a, double b, double Dc, double alpha, 
    double k_s, double mu_0, double v_0, double disp_0,
    const double vlp[], size_t num_vlp,
    const double sig_n[], size_t num_sig_n, 
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    this->p_impl = new RsfTauChgImpl(
        a, b, Dc, alpha,
        k_s, mu_0, v_0, disp_0,
        vlp, num_vlp,
        sig_n, num_sig_n,
        pp, num_pp,
        t, num_t
    );
}

RSF_tau_chg::RSF_tau_chg(const RSF_tau_chg& rsf)
{
    if (RsfTauChgImpl* rsf_p_impl = dynamic_cast<RsfTauChgImpl*>(rsf.p_impl)) {
        this->p_impl = new RsfTauChgImpl(*rsf_p_impl);
    }
    else {
        std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
    }
}

RSF_tau_chg& RSF_tau_chg::operator=(const RSF_tau_chg& rsf)
{
    if (this != &rsf) {
        delete this->p_impl;
        if (RsfTauChgImpl* rsf_p_impl = dynamic_cast<RsfTauChgImpl*>(rsf.p_impl)) {
            this->p_impl = new RsfTauChgImpl(*rsf_p_impl);
        }
        else {
            std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
        }
    }

    return *this;
}

void RSF_tau_chg::solve_rk4()
{
    if (RsfTauChgImpl* l_p_impl = dynamic_cast<RsfTauChgImpl*>(this->p_impl)) {
        l_p_impl->solve_rk4();
    }
    else {
        std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
    }
}

RSF_tau_chg::~RSF_tau_chg() 
{
    if (this->p_impl) {
        delete this->p_impl;
    }
}
