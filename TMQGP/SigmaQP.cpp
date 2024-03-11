#include "SigmaQP.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>

double QM_integrand(double om, double p, 
            double k, double x, double T, Interpolator2D & iImT, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l){
    double omp = eps2(k); // delta-function integral

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
    // double res = iImT(sqrt(k2), 
    //             ((om + omp > 0) - (om + omp < 0)) * Ecm) * iImG(k, omp);

    double res = iImT(sqrt(k2), Ecm 
    * ((om + omp > 0) - (om + omp < 0))
    ) 
    * (n_f(omp, T) + n_b(om + omp, T))
                * k * k / 4 / M_PI/ M_PI;
    // double res = iImT(sqrt(k2), Ecm) * (n_f(omp, T) + n_b(om + omp, T))
    //             * k * k / 4 / M_PI/ M_PI;
    // double res = k*k * iImT(sqrt(k2), sqrt(s))  / 4 / M_PI/ M_PI;

    ////// Angluar part //////

    // double Pl = gsl_sf_legendre_Pl(l, x);
    // All the legendre polynomials at x_cm = 1 are = 1

    return res;
}

double QP_x_int(double om, double p,
            double k, double T, Interpolator2D & iImT, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l){
    double res, err, res_m, err_m, res_p, err_p;

    funct i_func_x = [&](double x) -> double {
        return QM_integrand(om, p, k, x, T, iImT, eps1, eps2, debug, l) ;
    };

    
    integ_x.integrate(&i_func_x, -1, 1, res_p, err_p);
    return res_p;
}

double SigmaQP(double om, double p,
            double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l, double UL){
    double res, err;

    funct i_func_x = [&](double k) -> double {
        return QP_x_int(om, p, k, T, iImT, eps1, eps2, debug, l);
    };

    integ_x.integrate(&i_func_x, 0, UL, res, err);
    return res;
}

double ReSigmaQP(double om, double p,
            double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l, double UL){
    double res, err;

    funct i_func_x = [&](double k) -> double {
        return QP_x_int(om, p, k, T, iImT, eps1, eps2, debug, l);
    };

    integ_x.integrate(&i_func_x, 0, UL, res, err);
    return res;
}

double Re_QM_integrand(double om, double p, 
            double k, double x, double T, Interpolator2D & iImT, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l){
    double omp = eps2(k); // delta-function integral
    double s = pow(eps1(p) + eps2(k), 2.) - (p*p + k*k + 2*p*k*x);
    double Ecm = (pow(omp + om, 2.) - (p*p + k*k + 2*p*k*x));
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test disabling the spacelike integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Ecm < 0)
        return 0;
    Ecm = sqrt(Ecm);
    double m1sq = pow(eps1(0), 2.);
    double m2sq = pow(eps2(0), 2.);
    double k2 = 1/s * (pow(s - m1sq - m2sq, 2) / 4 - m1sq*m2sq);

    if (debug){
        printf("m1sq = %.2f, m2sq = %.2f, k2 = %.2f, s=%.2f", m1sq, m2sq, k2, s);
    }
    
    if (k2 < 0) return 0;

    double res = iImT(sqrt(k2), Ecm * ((om + omp > 0) - (om + omp < 0))) * (n_f(omp, T) + n_b(om + omp, T))
                * k * k / 4 / M_PI/ M_PI;

    return res;
}