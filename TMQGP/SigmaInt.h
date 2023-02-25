
#ifndef COMMANDS_H_
#define COMMANDS_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>

// Python code:
// def x_integrand(E, p, k, x):
//     k2 = k**2 + p**2 + 2*k*p*x
//     if k2 < 0:
//         return 0
//     return iImTs(E, sqrt(k2))

// def x_integral(E, p, k):
//     return quad(lambda x: x_integrand(E, p, k, x), -1, 1, epsrel=1e-2, epsabs=1e-2)[0]

double x_integrand(double x, double k, double E, double p, Interpolator2D iImT);

// auto * integ_x = new IntGSL<std::function<double(double)>>;
// auto * integ_k = new IntGSL<std::function<double(double)>>;
// auto * integ_E = new IntGSL<std::function<double(double)>>;

// Int_gsl_fixed * integ_x = new Int_gsl_fixed(-1, 1);
// Int_gsl_fixed * integ_k = new Int_gsl_fixed(1e-3, 5);
// Int_gsl_fixed * integ_E = new Int_gsl_fixed(1e-3, 5); 

IntGSL<std::function<double(double)>> integ_x;
IntGSL<std::function<double(double)>> integ_k;
IntGSL<std::function<double(double)>> integ_E;

Int_gsl_adaptive integ_Om;
IntGSL<std::function<double(double)>> integ_T;
// Int_gsl_fixed * integ_T = new Int_gsl_fixed(1e-3, 5);
Int_gsl_cauchy inter_cauchy;

double x_integral(double k, double E, double p, Interpolator2D iImT);
double k_integral(double E, double om, double p, Interpolator2D iImT, Interpolator2D iImG);
double k_integral_cm(double E, double om, double p, Interpolator2D iImT, Interpolator2D iImG);
double E_integral(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);
double E_integral_cm(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);
double E_int_first(double om, double p, double T, double k, Interpolator2D iImT, Interpolator2D iImG);
double Res_Kint(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);
double x_integral_cm(double E, double om, double p, double k, Interpolator2D iImT, Interpolator2D iImG);
double func_E_int_first(double e, double om, double p, double T, double k, Interpolator2D iImT, Interpolator2D iImG);

void get_E_int(double om, double T, Interpolator2D iImT, Interpolator2D iImG, double * p, int dimP, double * out, int dimOut);
// void get_T(double E, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, double * p, int dimP, 
//             double *out, int dimOut, double * out2, int dimOut2);
void get_T(double E, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, double * p, int dimP, 
            // double *out, int dimOut);
            std::complex<double> *out, int dimOut);
std::complex<double> T_solve(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
            double Lambda = 5);

std::complex<double> T_solve_test(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
            IntGSL<funct> integ_T, double Lambda = 5);

// double T_solveRe(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
//             double Lambda, IntGSL<funct> integ);

double T_solveRe(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
            double Lambda);

double ReSigmaKK_2D(double E, double q, Interpolator2D iImS);
double ReSigmaKK(double E, Interpolator iImS);

double sigma_bb(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);
double sigma_bf(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);
double sigma_ff(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);
double sigma_fb(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);

void get_test_gsl_interp(int N, double * p, int dimP, double * out, int dimOut);
//////////////////////////////////////////////////////////////////////

// Another integration variable choice
double sigma_bb2(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);
double sigma_bb3(double om, double p, double T, Interpolator2D iImT, Interpolator2D iImG);

double k_integral_cm2(double omp, double om, double p, Interpolator2D iImT, Interpolator2D iImG);
double x_integral_cm2(double omp, double om, double p, double k, Interpolator2D iImT, Interpolator2D iImG);


void get_test(double * p, int dimP, double * out, int dimOut);
void get_test_gsl(int N, double * p, int dimP, double * out, int dimOut);

double OmQ_F_om_int(double q, double T, Interpolator2D iImG, Interpolator2D iReG);
double OmQ_F(double T, Interpolator2D iImG, Interpolator2D iReG);
double OmQ_B_om_int(double q, double T, Interpolator2D iImG, Interpolator2D iReG);
double OmQ_B(double T, Interpolator2D iImG, Interpolator2D iReG);
double delta(double om, double q, Interpolator2D iImG, Interpolator2D iReG);

double OmS_F_om_int(double q, double T, Interpolator2D iImG, Interpolator2D iReG,
                                        Interpolator2D iImS, Interpolator2D iReS);

double OmS_F(double T, Interpolator2D iImG, Interpolator2D iReG,
                                        Interpolator2D iImS, Interpolator2D iReS);

std::complex<double> T_solve_BB(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
            double Lambda = 5);
std::complex<double> T_solve_BF(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D iReGqq, Interpolator2D iImGqq, 
            double Lambda = 5);
// class Runner {
//     public:
//         Runner();
        
//         Interpolator2D iImS;
//         Interpolator2D iGq;
//         Interpolator2D iImGqq;
//         Interpolator2D iReGqq;
// };

#endif