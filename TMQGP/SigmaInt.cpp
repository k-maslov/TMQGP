#include "SigmaInt.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>



// Define the integrators once
IntGSL<std::function<double(double)>> integ_x;
// Int_gsl_adaptive integ_x;
IntGSL<std::function<double(double)>> integ_k;
IntGSL<std::function<double(double)>> integ_E;
// Int_gsl_adaptive integ_E;

Int_gsl_adaptive integ_Om;
Int_gsl_adaptive integ_Om2;
IntGSL<std::function<double(double)>> integ_T;
// Int_gsl_fixed * integ_T = new Int_gsl_fixed(1e-3, 5);
Int_gsl_adaptive integ_T_adaptive;
Int_gsl_cauchy inter_cauchy;

double x_integrand(double x, double k, double E, double p, Interpolator2D & iImT){
    double k2 = k*k + p*p + 2*p*k*x;
    if (k2 < 0) return 0;
    return iImT(sqrt(k2), E);
}

double x_integral(double k, double E, double p, Interpolator2D & iImT){
    double res, err;

    funct i_func_x = [&](double x) -> double {return x_integrand(x, k, E, p, iImT);};
    integ_x.integrate(&i_func_x, -1, 1, res, err);
    return res;
}

double k_integral(double E, double om, double p, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    // printf("Hello!!!\n");
    funct i_func_k = [&](double k) -> double {
        //std::cerr << "Hello  " << k << std::endl;
        double res = k*k * x_integral(k, E, p, iImT) * iImG(k, E - om) /4/M_PI/M_PI;
        //std::cerr << "res  " << res << std::endl;
        // printf("%.3e \n", res);
        // printf('')
        return res;
    };

    integ_k.integrate(&i_func_k, 1e-3, 2, res, err);
    return res;
}

double n_f(double om, double T){
    // if (om < 0){
    //     return - n_f(-om, T);
    // }
    return 1/(exp(om/T) + 1);
}

double n_b(double om, double T){
    return real(1./(exp(om/T) - 1 + std::complex<double>(0, 1e-8)));
}

double E_integral(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();
    // printf("Hello!!!\n");
    funct i_func_e = [&](double e) -> double {
        //std::cerr << "Hello  " << k << std::endl;
        double res = k_integral(e, om, p, iImT, iImG) * (n_f(e - om, T) + n_b(e, T)) / M_PI / 2; // 1/2 from omega' integration
        //std::cerr << "res  " << res << std::endl;
        // printf("%.3e \n", res);
        // printf('')
        return res;
    };

    integ_E.integrate(&i_func_e, 1e-3, 2, res, err);
    return res;
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double x_integral_cm(double E, double om, double p, double k, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err;

    funct i_func_x = [&](double x) -> double {
        double s = E*E - (p*p + k*k + 2*p*k*x);
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test disabling the spacelike integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // if (s < 0)
            // return 0;
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

double x_integral_cm2(double omp, double om, double p, double k, 
        Interpolator2D & iImT, Interpolator2D & iImG, int debug){
    double res, err;

    funct i_func_x = [&](double x) -> double {
        double s = (om + omp)*(om + omp) - (p*p + k*k + 2*p*k*x);
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test disabling the spacelike integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (s < 0)
            return 0;
        double m1sq = pow(omp, 2) - k*k;
        double m2sq = pow(om, 2) - p*p;
        double k2 = 1/s * (pow(s - m1sq - m2sq, 2) / 4 - m1sq*m2sq);

        if (debug){
            printf("m1sq = %.2f, m2sq = %.2f, k2 = %.2f, s=%.2f", m1sq, m2sq, k2, s);
        }
        
        // if (om + omp < 0){
        //     return 0;
        // }
        if (k2 < 0) return 0;
        double res = k*k * iImT(sqrt(k2), ((om + omp > 0) - (om + omp < 0)) * sqrt(s))  / 4 / M_PI/ M_PI;
        // double res = k*k * iImT(sqrt(k2), sqrt(s))  / 4 / M_PI/ M_PI;
        return res;
    };
    integ_x.integrate(&i_func_x, -1, 1, res, err);
    return res;
}



double k_integral_cm2(double omp, double om, double p, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    // printf("Hello!!!\n");
    funct i_func_k = [&](double k) -> double {
        return x_integral_cm2(omp, om, p, k, iImT, iImG) *  iImG(k, omp);
    };
    integ_k.integrate(&i_func_k, 0, 5, res, err);
    return res;
}

double k_integral_cm(double E, double om, double p, Interpolator2D & iImT, Interpolator2D & iImG){
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

double E_integral_cm(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
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

double sigma_integrand_bb(double omp, double om, double p, double T, Interpolator2D & iImT, 
                        Interpolator2D & iImG){
    double res = k_integral_cm2(omp, om, p, iImT, iImG) * (n_b(omp, T) - n_b(omp + om, T));
    return res;
}

double sigma_bb(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err, res2, err2;
    gsl_set_error_handler_off();

    funct i_func_e = [&](double e) -> double {
        return sigma_integrand_bb(e, om, p, T, iImT, iImG);
    };

    double eps = 1e-2;
    // if (om > 0){
        // integ_E.integrate(&i_func_e, 0, om-eps, res, err);
        integ_E.integrate(&i_func_e, -5, 5, res, err2);
        res2 = 0;
    // }
    // else{
    //     return 0;
    //     // integ_E.integrate(&i_func_e, 0, 10, res, err);
    //     // res2 = 0;
    // }
    // integ_E.integrate(&i_func_e, om, 10, res, err);
    

    return res + res2;
}

void get_test(double * p, int dimP, double * out, int dimOut){
    int i, j;
    double h = 10./10000.;
    double res;

    omp_set_dynamic(1);
    omp_set_num_threads(16);
    #pragma omp parallel for private(i, j)// reduction(+:res)
    for (i = 0; i < dimP; i++){
        res = 0;
        for (j = 0; j < 10000; j++){
            res += h * (1 / (p[i]*p[i] + pow(-5 + h*j, 2)));
            // printf( "%.2f\n", res);
            // cout << p[i] << "   " << res << endl;
        }
        out[i] = res;
    }
}

void get_test_gsl(int N, double * p, int dimP, double * out, int dimOut){
    int i, j;
    double h = 10./10000.;
    double res, err;
    gsl_set_error_handler_off();

    omp_set_dynamic(0);
    omp_set_num_threads(N);
    IntGSL<std::function<double(double)>> integ;

    // std::vector<IntGSL<funct>> ints;
    // for (int i = 0; i < N; i++){
    //     ints.push_back(IntGSL<funct>());
    // }

    std::vector<Int_gsl_adaptive> ints;
    for (int i = 0; i < N; i++){
        ints.push_back(Int_gsl_adaptive());
    }

    int thr=0;
    #pragma omp parallel for schedule(static) shared(out) private(i, thr)
    for (i = 0; i < dimP; i++){
        // res = 0;
        double _p = p[i];
        thr = omp_get_thread_num();
        funct i_func = [=](double x) -> double {
            return 1/(_p*_p + x*x);
        };
        ints[thr].integrate(&i_func, -5, 5, res, err);
        out[i] = res;
    }
}


void get_test_gsl_interp(int N, double * p, int dimP, double * out, int dimOut){
    int i, j;
    double h = 10./10000.;
    double res, err;
    gsl_set_error_handler_off();

    omp_set_dynamic(0);
    omp_set_num_threads(N);
    IntGSL<std::function<double(double)>> integ;

    std::vector<IntGSL<funct>> ints;
    for (int i = 0; i < N; i++){
        ints.push_back(IntGSL<funct>());
    }

    double prange[100];
    double xrange[100];
    for (int i = 0; i< 100; i++){
        prange[i] = 4.8/100 * i;
        xrange[i] = -6 + 6/100 * i;
    }

    double func_range[100][100];

    auto func = [](double p, double x) -> double {
            return 1/(p*p + x*x);
    };
    
    for (int i = 0; i < 100; i++){
        for (int j = 0; j < 100; j++){
            func_range[i][j] = func(prange[i], xrange[j]);
        }
    }

    Interpolator2D interp(prange, 100, xrange, 100, *func_range, 100, 100);

    int thr=0;
    #pragma omp parallel for schedule(static) shared(out) private(i, thr)
    for (i = 0; i < dimP; i++){
        // res = 0;
        double _p = p[i];
        funct i_func = [&](double x) -> double {
            return interp(_p, x);
        };

        thr = omp_get_thread_num();
        ints[thr].integrate(&i_func, -5, 5, res, err);
        out[i] = res;
    }
}

double sigma_bb2(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err, res2, err2, res3, err3;
    res2 = 0;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        double res = k_integral_cm2(omp, om, p, iImT, iImG) * (n_b(omp, T) - n_b(omp + om, T));
        return res;
    };

    double eps = 1e-2;

    double l1, l2;
    if (om > 0){
        l1 = 0;
        l2 = om;
    }
    else{
        l1 = om;
        l2 = 0;
    }

    // if (om > 0){
        // integ_E.integrate(&i_func_e, 0, om-eps, res, err);
    

    // if (om > 0){
    integ_E.integrate(&i_func_e, -5, l1 - eps, res, err2);
    integ_E.integrate(&i_func_e, l1 + eps, l2 - eps, res2, err2);
    integ_E.integrate(&i_func_e, l2 + eps, 5, res3, err2);
    // }
    // else{
    //     integ_E.integrate(&i_func_e, -5, 0, res, err2);
    //     return -res;
    // }

    // integ_E.integrate(&i_func_e, om + eps, 5, res2, err2);
    // res2 = 0;
    // }
    // else{
        // return 0;
        // integ_E.integrate(&i_func_e, 0, 10, res, err);
        // res2 = 0;
    // }
    // integ_E.integrate(&i_func_e, om, 10, res, err);
    return res + res2 + res3;
}

double sigma_bb3(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err, res2, err2, res3, err3;
    res2 = 0;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        double res = k_integral_cm2(omp, om, p, iImT, iImG) * (n_b(omp, T) - n_b(omp + om, T));
        return res;
    };

    double eps = 1e-2;

    double l1, l2;
    if (om > 0){
        l1 = 0;
        l2 = om;
    }
    else{
        l1 = om;
        l2 = 0;
    }

    integ_E.integrate(&i_func_e, -5, 5, res, err2);
    // integ_E.integrate(&i_func_e, om + eps, 5, res2, err2);
    // integ_E.integrate(&i_func_e, l2 + eps, 5, res3, err2);
    return res;
}

double sigma_integrand_bf(double omp, double om, double p, double T, Interpolator2D & iImT, 
                        Interpolator2D & iImG){
    double res = k_integral_cm2(omp, om, p, iImT, iImG) * (n_f(omp, T) - n_f(omp + om, T));
    return res;
}


double sigma_bf(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err;
    // if (om < 0) return 0;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        return sigma_integrand_bf(omp, om, p, T, iImT, iImG);
    };

    integ_E.integrate(&i_func_e, -5, 5, res, err);
    return res;
}

double sigma_integrand_fb(double omp, double om, double p, double T, Interpolator2D & iImT, 
                        Interpolator2D & iImG){
    double res = k_integral_cm2(omp, om, p, iImT, iImG) * (n_b(omp, T) + n_f(omp + om, T));
    return res;
}

double sigma_fb(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err, res2, err2;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        return sigma_integrand_fb(omp, om, p, T, iImT, iImG);
    };


    // double eps = 1e-3;
    // if (om > 0){
    //     integ_E.integrate(&i_func_e, 0, om-eps, res2, err);
    //     integ_E.integrate(&i_func_e, om, 10, res, err2);
    //     // res2 = 0;
    // }
    // else{
    //     integ_E.integrate(&i_func_e, 0, 10, res, err);
    //     res2 = 0;
    // }
    
    // return res + res2;

    integ_E.integrate(&i_func_e, -5, 5, res, err);
    return res;
}

double sigma_ff(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double omp) -> double {
        double res = k_integral_cm2(omp, om, p, iImT, iImG) * (n_f(omp, T) + n_b(omp + om, T));
        return res;
    };

    integ_E.integrate(&i_func_e, -0, 5, res, err);
    return res;
}

double func_E_int_first(double e, double om, double p, double T, double k, Interpolator2D & iImT, Interpolator2D & iImG){
    return x_integral_cm(e, om, p, k, iImT, iImG) * (n_f(e - om, T) + n_b(e, T)) * iImG(k, e - om);
}

double E_int_first(double om, double p, double T, double k, Interpolator2D & iImT, Interpolator2D & iImG){
    double res, err;
    gsl_set_error_handler_off();
    funct i_func_e = [&](double e) -> double {
        double res = func_E_int_first(e, om, p, T, k, iImT, iImG);
        return res;
    };

    integ_E.integrate(&i_func_e, 1e-3, 5, res, err);
    return res;
}

double Res_Kint(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG){
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

void get_E_int(double om, double T, Interpolator2D & iImT, Interpolator2D & iImG, double * p, int dimP, double * out, int dimOut)
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

std::complex<double> T_solve(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iReGqq(k, E); 
    };

    funct i_func_im = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iImGqq(k, E); 
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);

    // cout << res1 << "  " << res2 << endl;
    std::complex<double> res(res1, res2);
    return (-sign*iVK(q) * iVK(q1) / (1.0 - res));
}


std::complex<double> T_solve_explicit(
    double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iReGqq(k, E); 
    };

    funct i_func_im = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iImGqq(k, E);
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);

    // cout << res1 << "  " << res2 << endl;
    std::complex<double> res(res1, res2);
    return (-sign*iVK(q) * iVK(q1) / ((E*E - pow(0.7, 2)) - res));
}

std::complex<double> J_solve_explicit(
    double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iReGqq(k, E);
    };

    funct i_func_im = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iImGqq(k, E);
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);

    // cout << res1 << "  " << res2 << endl;
    std::complex<double> res(res1, res2);
    return res;
}



std::complex<double> T_solve_test(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            IntGSL<funct> integ_T, double Lambda){
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

std::complex<double> T_solve_BB(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        double omk = iOmK(k); 
        // return 2/M_PI * k*k * omk*omk * -sign*iVK(k)*iVK(k) * iReGqq(k, E); 
        return 1/2/M_PI/M_PI * k*k * -sign*iVK(k)*iVK(k) * iReGqq(k, E); 
    };

    funct i_func_im = [&](double k) -> double {
        double omk = iOmK(k); 
        // return 2/M_PI * k*k * omk*omk* -sign*iVK(k)*iVK(k) * iImGqq(k, E); 
        return 1/2/M_PI/M_PI * k*k * -sign*iVK(k)*iVK(k) * iImGqq(k, E); 
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    std::complex<double> res(res1, res2);
    return (-sign*iVK(q) * iVK(q1) / (1.0 - res));
}

std::complex<double> T_solve_BF(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -sign*iVK(k)*iVK(k) * iReGqq(k, E); 
    };

    funct i_func_im = [&](double k) -> double {
        return 2/M_PI * k*k * -sign*iVK(k)*iVK(k) * iImGqq(k, E); 
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    std::complex<double> res(res1, res2);
    return (-sign*iVK(q) * iVK(q1) / (1.0 - res));
}

void get_T(double E, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, double * p, int dimP, 
            // std::complex<double> *out, int dimOut){//, double * out2, int dimOut2){
                // double *out, int dimOut){
                std::complex<double> *out, int dimOut){
    omp_set_dynamic(0);
    int N = 8;
    omp_set_num_threads(N);
    int i, thr;
    double res;

    // std::vector<IntGSL<funct>> vec_integ;

    // for (int i =0; i < N; i++){
    //     vec_integ.push_back(IntGSL<funct> ());
    // }
    
    std::vector<Interpolator> ints_O;
    std::vector<Interpolator> ints_V;
    std::vector<Interpolator2D> ints_Re;
    std::vector<Interpolator2D> ints_Im;
    std::vector<IntGSL<funct>> integs;

    for (int i = 0; i < N; i++){
        ints_O.push_back(Interpolator(iOmK.interp->x, iOmK.interp->size, iOmK.interp->y, iOmK.interp->size, "linear"));
        ints_V.push_back(Interpolator(iVK.interp->x, iVK.interp->size, iVK.interp->y, iVK.interp->size, "linear"));
        ints_Re.push_back(Interpolator2D(iReGqq));
        ints_Im.push_back(Interpolator2D(iImGqq));
        integs.push_back(IntGSL<funct>());
    }

    
    #pragma omp default(none)
    {
        // #pragma omp parallel private(i, thr) firstprivate(iVK, iOmK, iReGqq, iImGqq) shared(out)
        // Interpolator ivk = Interpolator(iVK.interp->x, iVK.interp->size, iVK.interp->y, iVK.interp->size, "linear");
        // Interpolator iomk = Interpolator(iVK.interp->x, iVK.interp->size, iVK.interp->y, iVK.interp->size, "linear");
        #pragma omp parallel for schedule(static) private(i, thr) shared(out, ints_V, ints_O, ints_Re, ints_Im, integs) //firstprivate(iVK, iOmK, iReGqq, iImGqq)
        for (i = 0; i < dimP; i++){    
            thr = omp_get_thread_num();
            out[i] = T_solve_test(E, p[i], p[i], T, ints_V[thr], ints_O[thr], ints_Re[thr], ints_Im[thr], integs[thr], 5.0);
        }
    }
}


double T_solveRe(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda){
    double res1, res1_l, res1_r, res2, err;

    // gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -iVK(k)*iVK(k) * iReGqq(k, E); 
    };

    // funct i_func_im = [&](double k) -> double {
    //     return 2/M_PI * k*k * -iVK(k)*iVK(k) * iImGqq(k, E); 
    // };

    IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    // integ.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    // std::complex<double> res(res1, res2);
    
    // auto out = real(-iVK(q) * iVK(q1) / (1.0 - res));
    return res1;
}

double ReSigmaKK_2D(double E, double q, Interpolator2D & iImS){
    gsl_set_error_handler_off();
    double res, err;

    funct i_func = [&] (double z) -> double{
        return iImS(q, z) / M_PI;
    };
    
    inter_cauchy.integrate(&i_func, -5, 5, E, res, err);
    return res;
}

double ReSigmaKK(double E,Interpolator & iImS){
    gsl_set_error_handler_off();
    double res, err;

    funct i_func = [&] (double z) -> double{
        return iImS(z) / M_PI;
    };
    
    inter_cauchy.integrate(&i_func, -5, 5, E, res, err);
    return res;
}
// double k_integral()

double delta(double om, double q, Interpolator2D & iImG, Interpolator2D & iReG){
    double im = iImG(q, om);
    double re = iReG(q, om);

    double res = atan(im/re);
    if (re > 0){
        res += M_PI;
    }
    return res;
}



double OmS2_F_om_int2(double omp, double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double om) -> double{
        return real(1./(omp - om + std::complex<double>(0, 1)*1e-3))*(n_f(om, T) - n_f(omp, T)) * (iImG(q, om) * iImS(q, omp)) / M_PI;
    };

    double res, err;
    integ_Om2.integrate(&func, 0, 5, res, err);
    return res;
}

double OmS2_F_om_int1(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double omp) -> double{
        return OmS2_F_om_int2(omp, q, T, iImG, iReG, iImS, iReS) / M_PI;
    };

    double res, err;
    integ_Om.integrate(&func, 0, 5, res, err);
    return res;
}

double OmS2_F(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double q) -> double{
        return q*q /2/M_PI/M_PI * OmS2_F_om_int1(q, T, iImG, iReG, iImS, iReS);
    };

    double res, err;
    integ_k.integrate(&func, 0, 5, res, err);
    return res;
}



double OmS_B_qfirst_q_int(double om, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    
    funct func = [&] (double q) -> double{
        return q*q /2/M_PI/M_PI * (iImG(q, om) * iReS(q, om) + iReG(q, om) * iImS(q, om));
    };

    double res, err;
    integ_Om.integrate(&func, -5, 5, res, err);
    return res;
}

double OmS_B_qfirst(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double om) -> double{
        return n_b(om, T) * OmS_B_qfirst_q_int(om, T, iImG, iReG, iImS, iReS) / M_PI;
    };

    double res, err;
    integ_k.integrate(&func, 0, 5, res, err);
    return res;
}


///////// Trying interpolating the denominator of quantites of interest


std::complex<double> T_solve_den(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, 
                    Interpolator2D & iReGqq, 
                    Interpolator2D & iImGqq, 
            double Lambda){

    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        double re = iReGqq(k, E);
        double im = iImGqq(k, E);
        return 2/M_PI * k*k * -iVK(k)*iVK(k) * (re / (re * re + im * im)); 
    };

    funct i_func_im = [&](double k) -> double {
        double re = iReGqq(k, E);
        double im = iImGqq(k, E);
        return 2/M_PI * k*k * -iVK(k)*iVK(k) * (im / (re * re + im * im)); 
    };

    // IntGSL<std::function<double(double)>> integ;

    integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
    integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    std::complex<double> res(res1, res2);
    return (-iVK(q) * iVK(q1) / (1.0 - res));
}


double test_integration(double a){
// 1. Define a lambda of the function to be integrated:
    funct i_func = [&] (double x){
        return 1 / (1 + x*x);
    };

// 2. Create the integrator
    Int_gsl_adaptive integ;

    double res, err;
    integ.integrate(&i_func, -a, a, res, err);

    return res;
}


double k_integral_QQ_func(double k, double E, double om, double p, Interpolator2D & iImG){
    return k*k*iImG(k, om) *  iImG(k, E - om);
}

double k_integral_QQ(double E, double om, double p, Interpolator2D & iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    // printf("Hello!!!\n");
    funct i_func_k = [&](double k) -> double {
        return k_integral_QQ_func(k, E, om, p, iImG);
    };
    integ_k.integrate(&i_func_k, 0, 5, res, err);
    return res;
}

double E_integral_QQ(double om, double p, double T, Interpolator2D & iImG){
    double res, err;
    // gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();
    // printf("Hello!!!\n");
    funct i_func_e = [&](double e) -> double {
        //std::cerr << "Hello  " << k << std:;
        double res = k_integral_QQ(e, om, p, iImG) * (tanh((e-om)/(2*T)) + tanh(e/2/T));
        //std::cerr << "res  " << res << std::endl;
        // printf("%.3e \n", res);
        // printf('')
        return res;
    };

    integ_Om.integrate(&i_func_e, -5, 5, res, err);
    return res;
}