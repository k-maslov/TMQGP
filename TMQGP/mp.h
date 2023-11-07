#ifndef MP_H_
#define MP_H_
#include "Interpolator.h"
#include "Integrator.h"
#include "IntGSL.h"
#include <complex>

Int_gsl_adaptive integ_test;

void get_sigma_ff(double * omrange, int dimOmrange, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
    int l, double * out, int dimOut);

void try_integration(double * omrange, int dimOmrange, int N, double * out, int dimOut);

double get_integ(double a, Int_gsl_adaptive integ);


void mpi_hw();

#endif