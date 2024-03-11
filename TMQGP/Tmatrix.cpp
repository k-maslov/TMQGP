#include "Tmatrix.h"
#include "SigmaInt.h"
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>

std::complex<double> x_solve(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign, int adaptive){
    double res1, res1_l, res1_r, res2, err;

    gsl_set_error_handler_off();
    funct i_func_re = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iReGqq(k, E);
    };

    funct i_func_im = [&](double k) -> double {
        return 2/M_PI * k*k * -sign * iVK(k)*iVK(k) * iImGqq(k, E); 
    };

    // IntGSL<std::function<double(double)>> integ;
    if (!adaptive){
        integ_T.integrate(&i_func_re, 1e-3, Lambda, res1, err);
        integ_T.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    }
    else{
        integ_T_adaptive.integrate(&i_func_re, 1e-3, Lambda, res1, err);
        integ_T_adaptive.integrate(&i_func_im, 1e-3, Lambda, res2, err);
    }
    // cout << res1 << "  " << res2 << endl;
    std::complex<double> res(res1, res2);
    return res;
}

std::complex<double> G2_an(double E, double q, double m){

}

std::complex<double> x_solve_analytic(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
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
