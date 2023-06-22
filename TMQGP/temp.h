
#ifndef TEMP_H_
#define TEMP_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>

// VACUUM CALCULATIONS FUNCTIONS

double p0_integral(double omega, double p, double m, double eps);

double FInt_no_inter(double omega, double p, double m, double eps);

double p_integral(double omega, double p0, double m, double eps);



double Integrand(double omega, double p0, double p, double m, double eps);

double full_int_using_interpo(double omega, double m, Interpolator2D & ImG);


double full_int(double omega, double m, double eps);



// 
//
//
// THERMAL CALCULATIONS FUNCTIONS
double Fint_integrand(double omega,double omega_1, double p, double p_1, double T, double m, double eps);

double FInt_no_inter_t(double omega, double p, double T, double m, double eps);

double full_int_t(double omega, double T, double m, double eps);

double FInt_t(double omega, double p, double T, Interpolator2D & ImG);

#endif