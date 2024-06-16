#ifndef RSF_TAU_CHG_H_
#define RSF_TAU_CHG_H_

#include <iostream>
#include <valarray>
#include <stdexcept>
#include "export_import.h"
#include "rsf_tau_fix.h"

class API RSF_tau_chg : public RSF_tau_fix
{
public:
    RSF_tau_chg();

    RSF_tau_chg(double a, double b, double Dc, double alpha);

    RSF_tau_chg(
        double a, double b, double Dc, double alpha, 
        double k_s, double mu_0, double v_0, double disp_0
    );

    RSF_tau_chg(
        double a, double b, double Dc, double alpha, 
        double k_s, double mu_0, double v_0, double disp_0,
        const double vlp[], size_t num_vlp,
        const double sig_n[], size_t num_sig_n, 
        const double pp[], size_t num_pp,
        const double t[], size_t num_t
    );

    RSF_tau_chg(const RSF_tau_chg& rsf);

    RSF_tau_chg& operator=(const RSF_tau_chg& rsf);

    // Set initial values of medium stiffness, friction coefficient, slip velocity, and slip displacement
    virtual void set_ini_vals(double mu_0, double v_0, double disp_0) override
    {
        throw std::logic_error("Must provide four arguments: k_s, mu_0, v_0 and disp_0.");
    }
    void set_ini_vals(double k_s, double mu_0, double v_0, double disp_0);

    // Set sequential independent variables
    virtual void set_indpt_val(
        const double sig_n[], size_t num_sig_n, 
        const double pp[], size_t num_pp,
        const double t[], size_t num_t
    ) override
    {
        throw std::logic_error("Must provide eight arguments: vlp[], num_vlp, sig_n[], num_sig_n, pp[], num_pp, t[], num_t.");
    }
    void set_indpt_val(
        const double vlp[], size_t num_vlp,
        const double sig_n[], size_t num_sig_n, 
        const double pp[], size_t num_pp,
        const double t[], size_t num_t
    );

    // Solve the equation group of the rate-and-state friction law
    // based on 4th Order Runge-Kutta method
    virtual void solve_rk4() override;

    virtual ~RSF_tau_chg();

protected:
    // Calculate the time derivative of shear stress (tau)
    virtual double dTau_dt(double disp_, double theta_, double tau_) const override;

private:
    double k_s;  // Stiffness of the medium around the fault

    // Sequential independent variable
    std::valarray<double>* vlp;  // Load point velocity
    
    // Evolving independent variable
    double vlp_evl;  // Evolving load point velocity
};

#endif  // RSF_TAU_CHG_H_
