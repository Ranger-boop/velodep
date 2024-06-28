#include <cmath>
#include <valarray>
#include <algorithm>
#include <stdexcept>
#include "rsf_tau_fix_impl.h"

void RsfTauFixImpl::set_const_paras(double a, double b, double Dc, double alpha)
{
    this->a = a;
    this->b = b;
    this->Dc = Dc;
    this->alpha = alpha;
}

void RsfTauFixImpl::set_ini_vals(double mu_0, double v_0, double disp_0)
{
    this->mu_0 = mu_0;
    this->v_0 = v_0;
    this->disp_0 = disp_0;
}

void RsfTauFixImpl::set_indpt_val(
    const double sig_n[], size_t num_sig_n,
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    if ((num_sig_n != num_pp) || (num_sig_n != num_t) || (num_pp != num_t)) {
        throw std::logic_error("The lengths of sig_n, pp, and time must be the same.");
    }
    else {
        this->sig_n = std::valarray<double>(sig_n, num_sig_n);
        this->pp = std::valarray<double>(pp, num_pp);
        this->t = std::valarray<double>(t, num_t);
    }
}

void RsfTauFixImpl::compute_sig_neff()
{
    if (!this->sig_n.size()) {
        throw std::logic_error("Please set independent values (sig_n, pp, and t).");
    }
    else {
        this->sig_neff = std::valarray<double>(this->sig_n - this->pp);
        this->sig_neff_evl = this->sig_neff[0];
    }
}

void RsfTauFixImpl::compute_dSig_neff_dt()
{
    if (!this->sig_neff.size()) {
        throw std::logic_error("Please compute effective normal stress (sig_neff).");
    }
    else if (!this->t.size()) {
        throw std::logic_error("Please set independent values (sig_n, pp, and t).");
    }
    else {
        this->dSig_neff_dt = std::valarray<double>(0., this->t.size());
        for (size_t i = 1; i < this->t.size(); i++) {
            this->dSig_neff_dt[i] = (this->sig_neff[i] - this->sig_neff[i - 1])
                / (this->t[i] - this->t[i - 1]);
        }

        this->dSig_neff_dt_evl = this->dSig_neff_dt[0];
    }
}

RsfTauFixImpl::RsfTauFixImpl()
{
    this->a = this->b = this->Dc = this->alpha = 0;

    this->mu_0 = this->v_0 = this->disp_0 = this->theta_0 = this->tau_0 = 0;

    //this->sig_n = this->pp = this->t = NULL;

    //this->sig_neff = this->dSig_neff_dt = NULL;
    this->sig_neff_evl = this->dSig_neff_dt_evl = 0;

    this->len_out_para = 0;
    this->vel = nullptr;
    this->disp = nullptr;
    this->theta = nullptr;
    this->tau = nullptr;
    this->mu = nullptr;
}

RsfTauFixImpl::RsfTauFixImpl(double a, double b, double Dc, double alpha)
{
    this->a = a;
    this->b = b;
    this->Dc = Dc;
    this->alpha = alpha;

    this->mu_0 = this->v_0 = this->disp_0 = this->theta_0 = this->tau_0 = 0;

    //this->sig_n = this->pp = this->t = nullptr;

    //this->sig_neff = this->dSig_neff_dt = nullptr;
    this->sig_neff_evl = this->dSig_neff_dt_evl = 0;

    this->len_out_para = 0;
    this->vel = nullptr;
    this->disp = nullptr;
    this->theta = nullptr;
    this->tau = nullptr;
    this->mu = nullptr;
}

RsfTauFixImpl::RsfTauFixImpl(
    double a, double b, double Dc, double alpha,
    double mu_0, double v_0, double disp_0
)
{
    this->a = a;
    this->b = b;
    this->Dc = Dc;
    this->alpha = alpha;

    //Set initial values
    this->mu_0 = mu_0;
    this->v_0 = v_0;
    this->disp_0 = disp_0;

    this->theta_0 = 0;
    this->tau_0 = 0;

    //this->sig_n = this->pp = this->t = nullptr;

    //this->sig_neff = this->dSig_neff_dt = nullptr;
    this->sig_neff_evl = this->dSig_neff_dt_evl = 0;

    this->len_out_para = 0;
    this->vel = nullptr;
    this->disp = nullptr;
    this->theta = nullptr;
    this->tau = nullptr;
    this->mu = nullptr;
}

RsfTauFixImpl::RsfTauFixImpl(
    double a, double b, double Dc, double alpha,
    double mu_0, double v_0, double disp_0,
    const double sig_n[], size_t num_sig_n,
    const double pp[], size_t num_pp,
    const double t[], size_t num_t
)
{
    this->a = a;
    this->b = b;
    this->Dc = Dc;
    this->alpha = alpha;

    // Set initial values
    this->mu_0 = mu_0;
    this->v_0 = v_0;
    this->disp_0 = disp_0;

    this->theta_0 = 0;
    this->tau_0 = 0;

    // Set independent values
    if ((num_sig_n != num_pp) || (num_sig_n != num_t) || (num_pp != num_t)) {
        throw std::logic_error("The lengths of sig_n, pp, and time must be the same.");
    }
    else {
        this->sig_n = std::valarray<double>(sig_n, num_sig_n);
        this->pp = std::valarray<double>(pp, num_pp);
        this->t = std::valarray<double>(t, num_t);
    }

    this->compute_sig_neff();
    this->compute_dSig_neff_dt();

    this->len_out_para = 0;
    this->vel = nullptr;
    this->disp = nullptr;
    this->theta = nullptr;
    this->tau = nullptr;
    this->mu = nullptr;
}

RsfTauFixImpl::RsfTauFixImpl(const RsfTauFixImpl& rsf)
{
    this->a = rsf.a;
    this->b = rsf.b;
    this->Dc = rsf.Dc;
    this->alpha = rsf.alpha;

    this->mu_0 = rsf.mu_0;
    this->v_0 = rsf.v_0;
    this->disp_0 = rsf.disp_0;
    this->theta_0 = rsf.theta_0;
    this->tau_0 = rsf.tau_0;

    this->sig_n = rsf.sig_n;
    this->pp = rsf.pp;
    this->t = rsf.t;

    this->sig_neff = rsf.sig_neff;
    this->dSig_neff_dt = rsf.dSig_neff_dt;

    this->sig_neff_evl = rsf.sig_neff_evl;
    this->dSig_neff_dt_evl = rsf.dSig_neff_dt_evl;

    this->len_out_para = rsf.len_out_para;
    this->vel = new double[rsf.len_out_para];
    this->disp = new double[rsf.len_out_para];
    this->theta = new double[rsf.len_out_para];
    this->tau = new double[rsf.len_out_para];
    this->mu = new double[rsf.len_out_para];
    std::copy(rsf.vel, rsf.vel + rsf.len_out_para, this->vel);
    std::copy(rsf.disp, rsf.disp + rsf.len_out_para, this->disp);
    std::copy(rsf.theta, rsf.theta + rsf.len_out_para, this->theta);
    std::copy(rsf.tau, rsf.tau + rsf.len_out_para, this->tau);
    std::copy(rsf.mu, rsf.mu + rsf.len_out_para, this->mu);
}

RsfTauFixImpl& RsfTauFixImpl::operator=(const RsfTauFixImpl& rsf)
{
    if (this == &rsf) {
        return *this;
    }

    delete[] this->vel;
    delete[] this->disp;
    delete[] this->theta;
    delete[] this->tau;
    delete[] this->mu;

    this->a = rsf.a;
    this->b = rsf.b;
    this->Dc = rsf.Dc;
    this->alpha = rsf.alpha;

    this->mu_0 = rsf.mu_0;
    this->v_0 = rsf.v_0;
    this->disp_0 = rsf.disp_0;
    this->theta_0 = rsf.theta_0;
    this->tau_0 = rsf.tau_0;

    this->sig_n = rsf.sig_n;
    this->pp = rsf.pp;
    this->t = rsf.t;

    this->sig_neff = rsf.sig_neff;
    this->dSig_neff_dt = rsf.dSig_neff_dt;

    this->sig_neff_evl = rsf.sig_neff_evl;
    this->dSig_neff_dt_evl = rsf.dSig_neff_dt_evl;

    this->len_out_para = rsf.len_out_para;
    this->vel = new double[rsf.len_out_para];
    this->disp = new double[rsf.len_out_para];
    this->theta = new double[rsf.len_out_para];
    this->tau = new double[rsf.len_out_para];
    this->mu = new double[rsf.len_out_para];
    std::copy(rsf.vel, rsf.vel + rsf.len_out_para, this->vel);
    std::copy(rsf.disp, rsf.disp + rsf.len_out_para, this->disp);
    std::copy(rsf.theta, rsf.theta + rsf.len_out_para, this->theta);
    std::copy(rsf.tau, rsf.tau + rsf.len_out_para, this->tau);
    std::copy(rsf.mu, rsf.mu + rsf.len_out_para, this->mu);

    return *this;
}

double RsfTauFixImpl::velocity(double theta_, double tau_) const
{
    return this->v_0 * std::exp(
        tau_ / this->sig_neff_evl / this->a
        - this->mu_0 / this->a
        - this->b / this->a * std::log(theta_ * this->v_0 / this->Dc)
    );
}

double RsfTauFixImpl::dD_dt(double disp_, double theta_, double tau_) const
{
    return this->velocity(theta_, tau_);
}

double RsfTauFixImpl::dTheta_dt(double disp_, double theta_, double tau_) const
{
    double vel = this->velocity(theta_, tau_);
    return 1
           - theta_ * vel / this->Dc
           - this->alpha
           * (theta_ / (this->b * this->sig_neff_evl))
           * this->dSig_neff_dt_evl;  // Aging law
}

double RsfTauFixImpl::dTau_dt(double disp_, double theta_, double tau_) const
{
    return 0.;
}

void RsfTauFixImpl::solve_rk4()
{
    if (!this->mu_0 || !this->v_0) {
        throw std::logic_error("mu_0 or v_0 is not set (at least one of them is 0).");
    }
    else if (!this->sig_n.size() || !this->pp.size() || !this->t.size()) {
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

size_t RsfTauFixImpl::get_len_out_para() const
{
    return this->len_out_para;
}

double* RsfTauFixImpl::get_disp()
{
    return this->disp;
}

double* RsfTauFixImpl::get_theta()
{
    return this->theta;
}

double* RsfTauFixImpl::get_tau()
{
    return this->tau;
}

double* RsfTauFixImpl::get_vel()
{
    return this->vel;
}

double* RsfTauFixImpl::get_mu()
{
    return this->mu;
}

RsfTauFixImpl::~RsfTauFixImpl()
{
    if (this->len_out_para) {
        delete[] this->vel;
        // this->vel = nullptr;
        delete[] this->disp;
        // this->disp = nullptr;
        delete[] this->theta;
        // this->theta = nullptr;
        delete[] this->tau;
        // this->tau = nullptr;
        delete[] this->mu;
        // this->mu = nullptr;
        this->vel = this->disp = this->theta = this->tau = this->mu = nullptr;
    }
}