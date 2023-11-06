#ifndef MP_H_
#define MP_H_
#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>

void get_sigma_ff(double * omrange, int dimOmrange, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
    int l, double * out, int dimOut);


#endif