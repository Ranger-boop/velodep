#ifndef RSF_TAU_FIX_H_
#define RSF_TAU_FIX_H_

#include <iostream>
#include <valarray>
#include "export_import.h"

class API RSF_tau_fix
{
public:
    RSF_tau_fix();

    RSF_tau_fix(double a, double b, double Dc, double alpha);

    RSF_tau_fix(
        double a, double b, double Dc, double alpha, 
        double mu_0, double v_0, double disp_0
    );

    RSF_tau_fix(
        double a, double b, double Dc, double alpha, 
        double mu_0, double v_0, double disp_0,
        const double sig_n[], size_t num_sig_n, 
        const double pp[], size_t num_pp,
        const double t[], size_t num_t
    );

    RSF_tau_fix(const RSF_tau_fix& rsf);

    RSF_tau_fix& operator=(const RSF_tau_fix& rsf);

    // Set constitutive patameters of the rate-and-state friciton law,
    // including a, b, Dc and alpha
    void set_const_paras(double a, double b, double Dc, double alpha);

    // Set initial values of friction coefficient, slip velocity, and slip displacement
    virtual void set_ini_vals(double mu_0, double v_0, double disp_0);

    // Set sequential independent variables
    virtual void set_indpt_val(
        const double sig_n[], size_t num_sig_n, 
        const double pp[], size_t num_pp,
        const double t[], size_t num_t
    );

    // Compute effective normal stress
    void compute_sig_neff();

    // Compute the time derivative of the effective normal stress
    void compute_dSig_neff_dt();

    // // Allocate memory for output parameters
    // void creat_output_paras();

    // Solve the equation group of the rate-and-state friction law
    // based on 4th Order Runge-Kutta method
    virtual void solve_rk4();

    size_t get_len_out_para();

    double* get_disp();

    double* get_theta();

    double* get_tau();

    double* get_vel();

    double* get_mu();

    virtual ~RSF_tau_fix();

protected:
    // Calculate slip velocity based on theta and tau
    double velocity(double theta_, double tau_) const;

    // Calculate the time derivative of slip displacement (i.e., slip velocity)
    double dD_dt(double disp_, double theta_, double tau_) const;

    // Calculate the time derivative of state variable (theta)
    double dTheta_dt(double disp_, double theta_, double tau_) const;

    // Calculate the time derivative of shear stress (tau)
    virtual double dTau_dt(double disp_, double theta_, double tau_) const;

    // Consitutive parameters
    double a;
    double b;
    double Dc;
    double alpha;
    // Predetermined parameters, initial values
    double mu_0;
    double v_0;
    double disp_0;
    double theta_0;
    double tau_0;
    // Sequential independent variables
    std::valarray<double>* sig_n;
    std::valarray<double>* pp;
    std::valarray<double>* t;
    // Variables computed based on independent variables
    std::valarray<double>* sig_neff;
    std::valarray<double>* dSig_neff_dt;
    // Evolving independent variables
    double sig_neff_evl;
    double dSig_neff_dt_evl;
    // output parameters
    double* vel;
    double* disp;
    double* theta;
    double* tau;
    double* mu;
    size_t len_out_para;  // Length of each output parameter
};

#endif // RSF_TAU_FIX_H_
