#include "G2.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>


double G2_conv_ff(double om, double q, double T, Interpolator2D & R1, Interpolator2D & R2, double Lambda){
    double res, err;

    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        return (1 - n_f(omp, T) - n_f(om - omp, T)) * R1(q, omp) * R2(q, om - omp);
        // return res;
    };

    integ_Om.integrate(&i_func_e, -Lambda, Lambda, res, err);
    return res;
}

double ReG2_conv_ff_integrand(double om1, double om, double q, double T, Interpolator2D & R1, Interpolator2D & R2, double Lambda){
    double res, err;

    gsl_set_error_handler_off();
    funct i_func_e = [&](double om2) -> double {
        return (1 - n_f(om1, T) - n_f(om2, T)) * R1(q, om1) * R2(q, om2);
        // return res;
    };

    inter_cauchy.integrate(&i_func_e, -Lambda, Lambda, om - om1, res, err);
    return res;
}

double ReG2_conv_ff(double om, double q, double T, Interpolator2D & R1, Interpolator2D & R2, double Lambda){
    double res, err, res_l, err_l, res_r, err_r;
    funct i_func_e = [&](double om1) -> double {
        return ReG2_conv_ff_integrand(om1, om, q, T, R1, R2);
        // return res;
    };
    double pole = sqrt(0.6*0.6 + q*q);
    double width = 2 * 0.05;
    integ_E.integrate(&i_func_e, -Lambda, pole - width, res_l, err_l);
    integ_E.integrate(&i_func_e, pole - width, pole + width, res, err);
    integ_E.integrate(&i_func_e, pole + width, Lambda, res_r, err_r);
    return res_l + res + res_r;
}

double ReG2_pole(double om, double q, double T, PoleInterpolator & R1, PoleInterpolator & R2, double Lambda){
    double res, err, res_l, err_l, res_r, err_r;
    funct i_func_e = [&](double om1) -> double {
        return ReG2_conv_ff_integrand(om1, om, q, T, R1, R2);
        // return res;
    };
    double pole = R1.iPole->operator()(q);
    double width = R1.iWidth->operator()(q);

    // cout << pole << "  " << width << endl;
    integ_E.integrate(&i_func_e, -Lambda, pole - width, res_l, err_l);
    integ_E.integrate(&i_func_e, pole - width, pole + width, res, err);
    integ_E.integrate(&i_func_e, pole + width, Lambda, res_r, err_r);
    return res_l + res + res_r;
}

double ReG2_conv_ff_integrand_subtr(double om1, double om, double q, double T, PoleInterpolator & R1, PoleInterpolator & R2, double Lambda){
    double res_l, err_l, res, err, res_r, err_r;

    funct i_func = [&](double om2) -> double {
        return (1 - n_f(om1, T) - n_f(om2, T)) * R1(q, om1) * R2(q, om2);
        // return res;
    };

    gsl_set_error_handler_off();
    funct i_func_e = [&](double om2) -> double {
        double r = (i_func(om2) - i_func(om - om1)) / (om - om1 - om2);
        return r;
        // return res;
    };

    double pole = R1.iPole->operator()(q);
    double width = R1.iWidth->operator()(q);
    // integ_x.integrate(&i_func_e, -Lambda, pole - width, res_l, err_l);
    // integ_x.integrate(&i_func_e, pole - width, pole + width, res, err);
    // integ_x.integrate(&i_func_e, pole + width, Lambda, res_r, err_r);
    integ_Om.integrate(&i_func_e, -Lambda, Lambda, res, err);
    double residue = - i_func(om - om1) * log(abs((Lambda - (om - om1)) / (-Lambda - (om - om1))));

    // return res_l + res + res_r + residue;
    return res + residue;
}


double ReG2_subtr(double om, double q, double T, PoleInterpolator & R1, PoleInterpolator & R2, double Lambda){
    double res, err, res_l, err_l, res_r, err_r;
    funct i_func_e = [&](double om1) -> double {
        return ReG2_conv_ff_integrand_subtr(om1, om, q, T, R1, R2, Lambda);
        // return res;
    };
    double pole = R1.iPole->operator()(q);
    double width = R1.iWidth->operator()(q);

    // cout << pole << "  " << width << endl;
    integ_E.integrate(&i_func_e, -Lambda, pole - width, res_l, err_l);
    integ_E.integrate(&i_func_e, pole - width, pole + width, res, err);
    integ_E.integrate(&i_func_e, pole + width, Lambda, res_r, err_r);
    return res_l + res + res_r;
}