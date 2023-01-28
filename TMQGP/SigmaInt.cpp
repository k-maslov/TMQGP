#include "SigmaInt.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>


double x_integrand(double x, double k, double E, double p, Interpolator2D iImT){
    double k2 = k*k + p*p + 2*p*k*x;
    if (k2 < 0) return 0;
    return iImT(sqrt(k2), E);
}

double x_integral(double k, double E, double p, Interpolator2D iImT){
    double res, err;

    funct i_func_x = [&](double x) -> double {return x_integrand(x, k, E, p, iImT);};
    integ_x.integrate(&i_func_x, -1, 1, res, err);
    return res;
}

double k_integral(double E, double om, double p, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    // printf("Hello!!!\n");
    funct i_func_k = [&](double k) -> double {
        //std::cerr << "Hello  " << k << std::endl;
        double res = k*k * x_integral(k, E, p, iImT) * iImG(k, E - om);
        //std::cerr << "res  " << res << std::endl;
        // printf("%.3e \n", res);
        // printf('')
        return res;
    };

    integ_k.integrate(&i_func_k, 1e-3, 2, res, err);
    return res;
}

double n_f(double om, double T){
    return 1/(exp(om/T) + 1);
}

double n_b(double om, double T){
    return real(1./(exp(om/T) - 1 + std::complex<double>(0, 1e-8)));
}

double E_integral(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();
    // printf("Hello!!!\n");
    funct i_func_e = [&](double e) -> double {
        //std::cerr << "Hello  " << k << std::endl;
        double res = k_integral(e, om, p, iImT, iImG) * (n_f(e - om, T) + 0*n_b(e, T)) / M_PI / 2; // 1/2 from omega' integration
        //std::cerr << "res  " << res << std::endl;
        // printf("%.3e \n", res);
        // printf('')
        return res;
    };

    integ_E.integrate(&i_func_e, 1e-3, 2, res, err);
    return res;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double x_integral_cm(double E, double om, double p, double k, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;

    funct i_func_x = [&](double x) -> double {
        double s = E*E - (p*p + k*k + 2*p*k*x);
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test disabling the spacelike integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (s < 0)
            return 0;
        double m1sq = pow(E - om, 2) - k*k;
        double m2sq = pow(om, 2) - p*p;
        double k2 = 1/s * (pow(s - m1sq - m2sq, 2) / 4 - m1sq*m2sq);
        if (k2 < 0) return 0;
        double res = k*k * iImT(sqrt(k2), sqrt(s))  / 4 / M_PI/ M_PI;
        return res;
    };
    integ_x.integrate(&i_func_x, -1, 1, res, err);
    return res;
}



double k_integral_cm(double E, double om, double p, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    // printf("Hello!!!\n");
    funct i_func_k = [&](double k) -> double {
        return x_integral_cm(E, om, p, k, iImT, iImG) *  iImG(k, E - om);
    };
    integ_k.integrate(&i_func_k, 0, 5, res, err);
    return res;
}

double E_integral_cm(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();
    // printf("Hello!!!\n");
    funct i_func_e = [&](double e) -> double {
        //std::cerr << "Hello  " << k << std::endl;
        double res = k_integral_cm(e, om, p, iImT, iImG) * (n_f(e - om, T) + n_b(e, T));
        //std::cerr << "res  " << res << std::endl;
        // printf("%.3e \n", res);
        // printf('')
        return res;
    };

    integ_E.integrate(&i_func_e, 0, 5, res, err);
    return res;
}


double sigma_bb(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG){
    double res, err, res2, err2;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double e) -> double {
        double res = k_integral_cm(e, om, p, iImT, iImG) * (n_b(e - om, T) - n_b(e, T));
        return res;
    };

    double eps = 1e-2;
    if (om > 0){
        // integ_E.integrate(&i_func_e, 0, om-eps, res, err);
        integ_E.integrate(&i_func_e, om, 5, res, err2);
        res2 = 0;
    }
    else{
        return 0;
        // integ_E.integrate(&i_func_e, 0, 10, res, err);
        // res2 = 0;
    }
    // integ_E.integrate(&i_func_e, om, 10, res, err);
    

    return res + res2;
}

double sigma_bf(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;
    if (om < 0) return 0;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double e) -> double {
        double res = k_integral_cm(e, om, p, iImT, iImG) * (n_f(e - om, T) - n_f(e, T));
        return res;
    };

    integ_E.integrate(&i_func_e, 0, 5, res, err);
    return res;
}

double sigma_fb(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG){
    double res, err, res2, err2;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double e) -> double {
        double res = k_integral_cm(e, om, p, iImT, iImG) * (n_b(e - om, T) + n_f(e, T));
        return res;
    };

    if (om > 0){
        // integ_E.integrate(&i_func_e, 0, om-eps, res, err);
        integ_E.integrate(&i_func_e, om, 10, res, err2);
        res2 = 0;
    }
    else{
        integ_E.integrate(&i_func_e, 0, 10, res, err);
        res2 = 0;
    }
    
    return res;
}

double sigma_ff(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double e) -> double {
        double res = k_integral_cm(e, om, p, iImT, iImG) * (n_f(e - om, T) + n_b(e, T));
        return res;
    };

    integ_E.integrate(&i_func_e, 0, 5, res, err);
    return res;
}



double func_E_int_first(double e, double om, double p, double T, double k, Interpolator2D iImT, Interpolator2D iImG){
    return x_integral_cm(e, om, p, k, iImT, iImG) * (n_f(e - om, T) + n_b(e, T)) * iImG(k, e - om);
}

double E_int_first(double om, double p, double T, double k, Interpolator2D iImT, Interpolator2D iImG){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double e) -> double {
        double res = func_E_int_first(e, om, p, T, k, iImT, iImG);
        return res;
    };

    integ_E.integrate(&i_func_e, 1e-3, 5, res, err);
    return res;
}

double Res_Kint(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG){
        double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    // printf("Hello!!!\n");
    funct i_func_k = [&](double k) -> double {
        return E_int_first(om, p, T, k, iImT, iImG);
    };
    integ_k.integrate(&i_func_k, 0, 5, res, err);
    return res;
}

void get_E_int(double om, double T, Interpolator2D iImT, Interpolator2D iImG, double * p, int dimP, double * out, int dimOut)
{
    omp_set_dynamic(0);
    omp_set_num_threads(16);
    int i;
    // #pragma omp parallel shared(p, out, iImT, iImG) private(i)
    {
        // #pragma omp parallel for
        #pragma omp parallel for schedule(dynamic) shared(p, out, iImT, iImG) private(i)
        for (i = 0; i < dimP; i++){
            out[i] = E_integral_cm(om, p[i], T, iImT, iImG);
        }
    }
}

void get_T(double E, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, double * p, int dimP, 
            // std::complex<double> *out, int dimOut){//, double * out2, int dimOut2){
                double *out, int dimOut){//, double * out2, int dimOut2){
    omp_set_dynamic(1);
    omp_set_num_threads(16);
    int i, thr;
    double res;
    // #pragma omp parallel shared(p, out, iImT, iImG) private(i)
    {
        // #pragma omp parallel for
        #pragma omp parallel private(thr)
        #pragma omp for private(i)
        for (i = 0; i < dimP; i++){
            out[i] = T_solveRe(E, p[i], p[i], T, iVK, iOmK, iReGqq, iImGqq);
            // cout << i << "   " << p[i] << "    " <<  out[i] << endl;
            // out2[i] = imag(res);
            thr = omp_get_thread_num();

            // cout << "thread " << thr << " doing" << i << endl;
        }
    }
}


std::complex<double> T_solve(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
            double Lambda){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -iVK(k)*iVK(k) * iReGqq(k, E); 
    };

    funct i_func_im = [&](double k) -> double {
        return 2/M_PI * k*k * -iVK(k)*iVK(k) * iImGqq(k, E); 
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    std::complex<double> res(res1, res2);
    return (-iVK(q) * iVK(q1) / (1.0 - res));
}

double T_solveRe(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
            double Lambda){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -iVK(k)*iVK(k) * iReGqq(k, E); 
    };

    funct i_func_im = [&](double k) -> double {
        return 2/M_PI * k*k * -iVK(k)*iVK(k) * iImGqq(k, E); 
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    std::complex<double> res(res1, res2);
    return real(-iVK(q) * iVK(q1) / (1.0 - res));
}
double ReSigmaKK_2D(double E, double q, Interpolator2D iImS){
    gsl_set_error_handler_off();
    double res, err;

    funct i_func = [&] (double z) -> double{
        return iImS(q, z) / M_PI;
    };
    
    inter_cauchy.integrate(&i_func, -5, 5, E, res, err);
    return res;
}

double ReSigmaKK(double E,Interpolator iImS){
    gsl_set_error_handler_off();
    double res, err;

    funct i_func = [&] (double z) -> double{
        return iImS(z) / M_PI;
    };
    
    inter_cauchy.integrate(&i_func, -5, 5, E, res, err);
    return res;
}
// double k_integral()