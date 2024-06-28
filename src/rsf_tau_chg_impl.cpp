#include <cmath>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include "rsf_tau_chg_impl.h"

void RsfTauChgImpl::set_ini_vals(double k_s, double mu_0, double v_0, double disp_0)
{
    RsfTauFixImpl::set_ini_vals(mu_0, v_0, disp_0);
    this->k_s = k_s;
}

void RsfTauChgImpl::set_indpt_val(
    const double vlp[], size_t num_vlp,
    const double sig_n[], size_t num_sig_n,
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    if (num_vlp != num_t) {
        throw std::logic_error("The lengths of vlp, sig_n, pp, and time must be the same!");
    }
    else {
        RsfTauFixImpl::set_indpt_val(
            sig_n, num_sig_n,
            pp, num_pp,
            t, num_t
        );
        this->vlp = std::valarray<double>(vlp, num_vlp);
        this->vlp_evl = vlp[0];
    }
}

RsfTauChgImpl::RsfTauChgImpl() : RsfTauFixImpl()
{
    this->k_s = 0;
    //this->vlp = nullptr;
    this->vlp_evl = 0;
}

RsfTauChgImpl::RsfTauChgImpl(
    double a, double b, double Dc, double alpha
) : RsfTauFixImpl(a, b, Dc, alpha)
{
    this->k_s = 0;
    //this->vlp = nullptr;
    this->vlp_evl = 0;
}

RsfTauChgImpl::RsfTauChgImpl(
    double a, double b, double Dc, double alpha,
    double k_s, double mu_0, double v_0, double disp_0
) : RsfTauFixImpl(a, b, Dc, alpha, mu_0, v_0, disp_0)
{
    this->k_s = k_s;
    //this->vlp = nullptr;
    this->vlp_evl = 0;
}

RsfTauChgImpl::RsfTauChgImpl(
    double a, double b, double Dc, double alpha,
    double k_s, double mu_0, double v_0, double disp_0,
    const double vlp[], size_t num_vlp,
    const double sig_n[], size_t num_sig_n,
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
) : RsfTauFixImpl(
    a, b, Dc, alpha,
    mu_0, v_0, disp_0,
    sig_n, num_sig_n,
    pp, num_pp,
    t, num_t
)
{
    this->k_s = k_s;
    if (num_vlp != num_t) {
        throw std::logic_error("The lengths of vlp, sig_n, pp, and time must be the same!");
    }
    else {
        this->vlp = std::valarray<double>(vlp, num_vlp);
        this->vlp_evl = vlp[0];
    }
}

RsfTauChgImpl::RsfTauChgImpl(const RsfTauChgImpl& rsf) : RsfTauFixImpl(rsf)
{
    this->k_s = rsf.k_s;
    this->vlp = std::valarray<double>(rsf.vlp);
    this->vlp_evl = rsf.vlp_evl;
}

RsfTauChgImpl& RsfTauChgImpl::operator=(const RsfTauChgImpl& rsf)
{
    if (this == &rsf) {
        return *this;
    }

    RsfTauFixImpl::operator=(rsf);
    this->k_s = rsf.k_s;
    this->vlp = std::valarray<double>(rsf.vlp);
    this->vlp_evl = rsf.vlp_evl;

    return *this;
}

double RsfTauChgImpl::dTau_dt(double disp_, double theta_, double tau_) const
{
    double vel = this->velocity(theta_, tau_);
    return this->k_s * (this->vlp_evl - vel);
}

void RsfTauChgImpl::solve_rk4()
{
    if (!this->mu_0 || !this->v_0) {
        throw std::logic_error("mu_0 or v_0 is not set (at least one of them is 0).");
    }
    else if (!this->vlp.size() || !this->sig_n.size() || !this->pp.size() || !this->t.size()) {
        throw std::logic_error("Independent values are not set (at least one of them is NULL).");
    }
    else if (!this->sig_neff.size() || !this->dSig_neff_dt.size()) {
        throw std::logic_error("Effective normal stress (sig_neff) or the time derivative of effective normal stress (dSig_neff_dt) is not computed.");
    }
    else {
        this->theta_0 = this->Dc / this->v_0;
        this->tau_0 = this->sig_neff[0] * this->mu_0;

        this->len_out_para = this->t.size();
        this->vel = new double[this->len_out_para];
        this->disp = new double[this->len_out_para];
        this->theta = new double[this->len_out_para];
        this->tau = new double[this->len_out_para];
        this->mu = new double[this->len_out_para];
    }
    this->vel[0] = this->v_0;
    this->theta[0] = this->theta_0;
    this->tau[0] = this->tau_0;
    this->mu[0] = this->mu_0;
    this->disp[0] = this->disp_0;

    double dt;
    double k1_disp, k1_theta, k1_tau;
    double k2_disp, k2_theta, k2_tau;
    double k3_disp, k3_theta, k3_tau;
    double k4_disp, k4_theta, k4_tau;
    // Runge-Kutta 4th Order Iteration Loop
    for (size_t i = 1; i < this->len_out_para; i++) {
        // update the evolving independent variables
        dt = this->t[i] - this->t[i - 1];
        this->sig_neff_evl = this->sig_neff[i];
        this->vlp_evl = this->vlp[i];
        this->dSig_neff_dt_evl = this->dSig_neff_dt[i];
        // k1------------------------------------------------------------------------
        k1_disp = this->dD_dt(
            this->disp[i - 1], this->theta[i - 1], this->tau[i - 1]
        );
        k1_theta = this->dTheta_dt(
            this->disp[i - 1], this->theta[i - 1], this->tau[i - 1]
        );
        k1_tau = this->dTau_dt(
            this->disp[i - 1], this->theta[i - 1], this->tau[i - 1]
        );
        // k2------------------------------------------------------------------------
        k2_disp = this->dD_dt(
            this->disp[i - 1] + dt * k1_disp / 2,
            this->theta[i - 1] + dt * k1_theta / 2,
            this->tau[i - 1] + dt * k1_tau / 2
        );
        k2_theta = this->dTheta_dt(
            this->disp[i - 1] + dt * k1_disp / 2,
            this->theta[i - 1] + dt * k1_theta / 2,
            this->tau[i - 1] + dt * k1_tau / 2
        );
        k2_tau = this->dTau_dt(
            this->disp[i - 1] + dt * k1_disp / 2,
            this->theta[i - 1] + dt * k1_theta / 2,
            this->tau[i - 1] + dt * k1_tau / 2
        );
        // k3------------------------------------------------------------------------
        k3_disp = this->dD_dt(
            this->disp[i - 1] + dt * k2_disp / 2,
            this->theta[i - 1] + dt * k2_theta / 2,
            this->tau[i - 1] + dt * k2_tau / 2
        );
        k3_theta = this->dTheta_dt(
            this->disp[i - 1] + dt * k2_disp / 2,
            this->theta[i - 1] + dt * k2_theta / 2,
            this->tau[i - 1] + dt * k2_tau / 2
        );
        k3_tau = this->dTau_dt(
            this->disp[i - 1] + dt * k2_disp / 2,
            this->theta[i - 1] + dt * k2_theta / 2,
            this->tau[i - 1] + dt * k2_tau / 2
        );
        // k4------------------------------------------------------------------------
        k4_disp = this->dD_dt(
            this->disp[i - 1] + dt * k3_disp,
            this->theta[i - 1] + dt * k3_theta,
            this->tau[i - 1] + dt * k3_tau
        );
        k4_theta = this->dTheta_dt(
            this->disp[i - 1] + dt * k3_disp,
            this->theta[i - 1] + dt * k3_theta,
            this->tau[i - 1] + dt * k3_tau
        );
        k4_tau = this->dTau_dt(
            this->disp[i - 1] + dt * k3_disp,
            this->theta[i - 1] + dt * k3_theta,
            this->tau[i - 1] + dt * k3_tau
        );
        // Assemble dependent variables----------------------------------------------
        this->disp[i] = this->disp[i - 1] + dt / 6 * (
            k1_disp + 2 * k2_disp + 2 * k3_disp + k4_disp
            );
        this->theta[i] = this->theta[i - 1] + dt / 6 * (
            k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta
            );
        this->tau[i] = this->tau[i - 1] + dt / 6 * (
            k1_tau + 2 * k2_tau + 2 * k3_tau + k4_tau
            );
        // Compute other variables---------------------------------------------------
        this->vel[i] = (this->disp[i] - this->disp[i - 1]) / dt;
        this->mu[i] = this->tau[i] / this->sig_neff_evl;
    }
}

RsfTauChgImpl::~RsfTauChgImpl() { }