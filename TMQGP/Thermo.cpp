#include "Thermo.h"
#include "SigmaInt.h"
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>


double OmQ_F_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG){
    funct func = [&] (double om) -> double{
        return n_f(om, T) * delta(om, q, iImG, iReG) / M_PI;
    };
    gsl_set_error_handler_off();

    double res, err;
    integ_Om.integrate(&func, -1, 5, res, err);
    return res;
}

double OmQ_F(double T, Interpolator2D & iImG, Interpolator2D & iReG){
    funct func = [&] (double q) -> double{
        return q*q /2/M_PI/M_PI * OmQ_F_om_int(q, T, iImG, iReG);
    };

    double res, err;
    integ_k.integrate(&func, 0, 3, res, err);
    return res;
}

double OmQ_B_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG){
    funct func = [&] (double om) -> double{
        return n_b(om, T) * delta(om, q, iImG, iReG) / M_PI;
    };

    double res, err;
    integ_Om.integrate(&func, -5, 5, res, err);
    return res;
}

double OmQ_B(double T, Interpolator2D & iImG, Interpolator2D & iReG){
    funct func = [&] (double q) -> double{
        return q*q /2/M_PI/M_PI * OmQ_B_om_int(q, T, iImG, iReG);
    };

    double res, err;
    integ_k.integrate(&func, 0, 5, res, err);
    return res;
}

double OmS_F_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double om) -> double{
        return n_f(om, T) * (iImG(q, om) * iReS(q, om) + iReG(q, om) * iImS(q, om)) / M_PI;
    };

    double res, err;
    integ_Om.integrate(&func, -5, 5, res, err);
    return res;
}

double OmS_F(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double q) -> double{
        return q*q /2/M_PI/M_PI * OmS_F_om_int(q, T, iImG, iReG, iImS, iReS);
    };

    double res, err;
    integ_k.integrate(&func, 0, 5, res, err);
    return res;
}

double OmS_B_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double om) -> double{
        return n_b(om, T) * (iImG(q, om) * iReS(q, om) + iReG(q, om) * iImS(q, om)) / M_PI;
    };

    double res, err;
    integ_Om.integrate(&func, -5, 5, res, err);
    return res;
}

double OmS_B(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS){
    funct func = [&] (double q) -> double{
        return q*q /2/M_PI/M_PI * OmS_B_om_int(q, T, iImG, iReG, iImS, iReS);
    };

    double res, err;
    integ_k.integrate(&func, 0, 5, res, err);
    return res;
}