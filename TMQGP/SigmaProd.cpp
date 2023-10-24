#include "SigmaProd.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>

double x_cm_onshell_integrand(double x, double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l){
    double s = pow(eps1(p) + eps2(k), 2.) - (p*p + k*k + 2*p*k*x);
    double Ecm = (pow(omp + om, 2.) - (p*p + k*k + 2*p*k*x));
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test disabling the spacelike integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Ecm < 0)
        return 0;
    Ecm = sqrt(Ecm);
    double m1sq = pow(eps1(0), 2.);
    double m2sq = pow(eps2(0), 2.);
    double k2 = 1/s * (pow(s - m1sq - m2sq, 2) / 4 - m1sq*m2sq);
    // if (om + omp < 0){
    //     return 0;
    // }

    if (debug){
        printf("m1sq = %.2f, m2sq = %.2f, k2 = %.2f, s=%.2f", m1sq, m2sq, k2, s);
    }
    if (k2 < 0) return 0;

    // double res = k*k * iImT(sqrt(k2), Ecm) / 4 / M_PI/ M_PI;
    double res = k*k * iImT(sqrt(k2), 
                ((om + omp > 0) - (om + omp < 0)) * Ecm)  / 4 / M_PI/ M_PI;

    // double res = k*k * iImT(sqrt(k2), sqrt(s))  / 4 / M_PI/ M_PI;

    ////// Angluar part //////

    // double Pl = gsl_sf_legendre_Pl(l, x);
    // All the legendre polynomials at x_cm = 1 are = 1

    return res;
}

double x_integral_cm_onshell(double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l){
    double res, err;

    funct i_func_x = [&](double x) -> double {
        // return x_cm_onshell_integrand(x, omp, om, p, k, iImT, iImG, eps1, eps2, debug, l);
            double s = pow(eps1(p) + eps2(k), 2.) - (p*p + k*k + 2*p*k*x);
    double Ecm = (pow(omp + om, 2.) - (p*p + k*k + 2*p*k*x));
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test disabling the spacelike integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Ecm < 0)
        return 0;
    Ecm = sqrt(Ecm);
    double m1sq = pow(eps1(0), 2.);
    double m2sq = pow(eps2(0), 2.);
    double k2 = 1/s * (pow(s - m1sq - m2sq, 2) / 4 - m1sq*m2sq);
    // if (om + omp < 0){
    //     return 0;
    // }

    if (debug){
        printf("m1sq = %.2f, m2sq = %.2f, k2 = %.2f, s=%.2f", m1sq, m2sq, k2, s);
    }
    if (k2 < 0) return 0;

    // double res = k*k * iImT(sqrt(k2), Ecm) / 4 / M_PI/ M_PI;
    double res = k*k * iImT(sqrt(k2), 
                ((om + omp > 0) - (om + omp < 0)) * Ecm)  / 4 / M_PI/ M_PI;

    // double res = k*k * iImT(sqrt(k2), sqrt(s))  / 4 / M_PI/ M_PI;

    ////// Angluar part //////

    // double Pl = 1;
    // if (l == 1) Pl = x;
    // if (l == 2) Pl = 1 - 3*pow(x, 2);
    
    return res;
    };

    integ_x.integrate(&i_func_x, -1, 1, res, err);
    return res;
}

double k_integral_onshell(double omp, double om, double p, 
Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
int l){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    // printf("Hello!!!\n");
    funct i_func_k = [&](double k) -> double {
        return x_integral_cm_onshell(omp, om, p, k, iImT, iImG, eps1, eps2, 0, l) *  iImG(k, omp);
    };
    integ_k.integrate(&i_func_k, 0, 3, res, err);
    return res;
}

double sigma_ff_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
    int l){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        double res = k_integral_onshell(omp, om, p, iImT, iImG, eps1, eps2, l) *
         (n_f(omp, T) + n_b(omp + om, T));
        return res;
    };

    integ_E.integrate(&i_func_e, -1, 4, res, err);
    return res;
}

double sigma_bb_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        double res = k_integral_onshell(omp, om, p, iImT, iImG, eps1, eps2, l) *
         (n_b(omp, T) - n_b(omp + om, T));
        return res;
    };

    integ_E.integrate(&i_func_e, 0, 4, res, err);
    return res;
}

double sigma_bf_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        double res = k_integral_onshell(omp, om, p, iImT, iImG, eps1, eps2, l) *
         (n_f(omp, T) - n_f(omp + om, T));
        return res;
    };

    integ_E.integrate(&i_func_e, -1, 4, res, err);
    return res;
}

double sigma_fb_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        double res = k_integral_onshell(omp, om, p, iImT, iImG, eps1, eps2, l) *
         (n_b(omp, T) + n_f(omp + om, T));
        return res;
    };

    integ_E.integrate(&i_func_e, 0, 4, res, err);
    return res;
}

double sigma_tot(double om, double p, double T, TMArray TMs, Interpolator2D & iImG_Q, Interpolator2D & iImG_G){
    // double res, err;


    // funct i_func_e = [&](double omp) -> double {
    //     double r = 0;
    //     double n_j, n_2b;
    //     for (auto TM : TMs){
    //         double mult = 1;
    //         if (TM.stat_j == "f"){
    //             n_j = n_f(omp, T);
    //         }
    //         else{
    //             n_j = n_b(omp, T);
    //         }

    //         if (TM.stat_i == "b") mult = -1;

    //         if ((TM.stat_i == "f" && TM.stat_j == "f") || (TM.stat_i == "b" && TM.stat_j == "b")){
    //             n_2b = n_b(omp + om, T);
    //         }
    //         else{
    //             n_2b = n_f(omp + om, T);
    //         }

    //         if (TM.stat_j == "b"){
    //             // if (omp < 0){
    //             //     r += 0;
    //             //     cout << "Im here\n" << endl;
    //             // }
    //             // else{
    //                 double rr = TM.Nf * TM.da * TM.ds / TM.d * k_integral_onshell(omp, om, p, *TM.iImT, iImG_G, TM.eps_i, TM.eps_j, 0)
    //                 * (n_j + mult*n_2b);
    //                 // cout << "rr = " << rr << endl;
    //                 r += rr;
    //             // }
    //         }
    //         else{
    //             r += TM.Nf * TM.da * TM.ds / TM.d * k_integral_onshell(omp, om, p, *TM.iImT, iImG_Q, TM.eps_i, TM.eps_j, 0)
    //             * (n_j + mult*n_2b);
    //         }
    //     }

    //     return r;
    // };
    
    // integ_E.integrate(&i_func_e, 0, 4, res, err);
    // return res;
    return -1;
}