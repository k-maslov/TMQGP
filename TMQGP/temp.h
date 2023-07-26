
#ifndef TEMP_H_
#define TEMP_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>



IntGSL<std::function<double(double)>> integ_kt;
IntGSL<std::function<double(double)>> integ_Et;

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

double FInt_t(double omega, double p, double T, Interpolator2D & ImG);


double F_p(double omega,double q,double T,double m,double eps, double lambda = 0.651);

double F_s(double omega,double q,double T,double m,double eps, double lambda = 0.651);

//double Fp_real(double omega, double q, double T, double m , double eps);

double Re_meson(double omega,double lambda4, Interpolator & iImS);
//double full_integrand_t(double p0, double omega, double p, double T, double m, double eps);



double imag_pi_inter1( double omega,double T,double m, double lambda, Interpolator2D & ImG);
//
//
// KITAZAWA 2014 CALCULATION
double Sigma_re_plus1(double p0, double T, double mf, double mb, Interpolator2D & ImG);

double Sigma_re_minus1(double p0, double T, double mf, double mb, Interpolator2D & ImG);

// QUARK LOOPS
complex<double> G0(double p0, double p, double m, double eps);

complex<double> Gv(double p0, double p, double m, double eps);

complex<double> Gm(double p0, double p, double m, double eps);

double Impi(double p, double omega, double T, double m, double eps);

double Impi_int(double omega, double T, double m, double eps);

// Harada Nemoto

double Impi2(double p, double q, double omega, double T, double m, double eps);

double Impi_int2(double q , double omega, double T, double m, double eps);
#endif