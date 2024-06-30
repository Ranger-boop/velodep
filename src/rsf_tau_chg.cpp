#include <cmath>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include <typeinfo>
#include "rsf_tau_chg.h"
#include "rsf_tau_chg_impl.h"

void RsfTauChg::set_ini_vals(double k_s, double mu_0, double v_0, double disp_0)
{
    if (RsfTauChgImpl* l_p_impl = dynamic_cast<RsfTauChgImpl*>(this->p_impl)) {
        l_p_impl->set_ini_vals(k_s, mu_0, v_0, disp_0);
    }
    else {
        std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
    }
}

void RsfTauChg::set_indpt_val(
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

RsfTauChg::RsfTauChg()
{
    this->p_impl = new RsfTauChgImpl();
}

RsfTauChg::RsfTauChg(
    double a, double b, double Dc, double alpha
)
{
    this->p_impl = new RsfTauChgImpl(a, b, Dc, alpha);
}

RsfTauChg::RsfTauChg(
    double a, double b, double Dc, double alpha, 
    double k_s, double mu_0, double v_0, double disp_0
)
{
    this->p_impl = new RsfTauChgImpl(
        a, b, Dc, alpha,
        k_s, mu_0, v_0, disp_0
    );
}

RsfTauChg::RsfTauChg(
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

RsfTauChg::RsfTauChg(const RsfTauChg& rsf)
{
    if (RsfTauChgImpl* rsf_p_impl = dynamic_cast<RsfTauChgImpl*>(rsf.p_impl)) {
        this->p_impl = new RsfTauChgImpl(*rsf_p_impl);
    }
    else {
        std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
    }
}

RsfTauChg& RsfTauChg::operator=(const RsfTauChg& rsf)
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

//void RsfTauChg::solve_rk4()
//{
//    //if (RsfTauChgImpl* l_p_impl = dynamic_cast<RsfTauChgImpl*>(this->p_impl)) {
//    //    l_p_impl->solve_rk4();
//    //}
//    //else {
//    //    std::cerr << "Error: dynamic cast failed in class RsfTauChg." << std::endl;
//    //}
//    this->p_impl->solve_rk4();
//}

RsfTauChg::~RsfTauChg() 
{
    //if (this->p_impl) {
    //    delete this->p_impl;
    //    this->p_impl = nullptr;
    //}
}
