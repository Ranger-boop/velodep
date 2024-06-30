#ifndef RSF_TAU_CHG_H_
#define RSF_TAU_CHG_H_

#include <iostream>
#include <valarray>
#include <stdexcept>
#include "export_import.h"
#include "rsf_tau_fix.h"

class API RsfTauChg : public RsfTauFix
{
public:
    RsfTauChg();

    RsfTauChg(double a, double b, double Dc, double alpha);

    RsfTauChg(
        double a, double b, double Dc, double alpha, 
        double k_s, double mu_0, double v_0, double disp_0
    );

    RsfTauChg(
        double a, double b, double Dc, double alpha, 
        double k_s, double mu_0, double v_0, double disp_0,
        const double vlp[], size_t num_vlp,
        const double sig_n[], size_t num_sig_n, 
        const double pp[], size_t num_pp,
        const double t[], size_t num_t
    );

    RsfTauChg(const RsfTauChg& rsf);

    RsfTauChg& operator=(const RsfTauChg& rsf);

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
    // virtual void solve_rk4() override;

    virtual ~RsfTauChg();
};

#endif  // RSF_TAU_CHG_H_
