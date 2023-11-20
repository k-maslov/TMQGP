
#ifndef SIGMAINT_H_
#define SIGMAINT_H_

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

double x_integrand(double x, double k, double E, double p, Interpolator2D & iImT);

// auto * integ_x = new IntGSL<std::function<double(double)>>;
// auto * integ_k = new IntGSL<std::function<double(double)>>;
// auto * integ_E = new IntGSL<std::function<double(double)>>;

// Int_gsl_fixed * integ_x = new Int_gsl_fixed(-1, 1);
// Int_gsl_fixed * integ_k = new Int_gsl_fixed(1e-3, 5);
// Int_gsl_fixed * integ_E = new Int_gsl_fixed(1e-3, 5); 

extern IntGSL<std::function<double(double)>> integ_x;
extern IntGSL<std::function<double(double)>> integ_k;
extern IntGSL<std::function<double(double)>> integ_E;
 
extern Int_gsl_adaptive integ_Om;
extern Int_gsl_adaptive integ_Om2;
// extern IntGSL<std::function<double(double)>> integ_T;
extern Int_gsl_adaptive integ_T;
 // Int_gsl_fixed * integ_T = new Int_gsl_fixed(1e-3, 5);
extern Int_gsl_cauchy inter_cauchy;

double x_integral(double k, double E, double p, Interpolator2D & iImT);
double k_integral(double E, double om, double p, Interpolator2D & iImT, Interpolator2D & iImG);
double k_integral_cm(double E, double om, double p, Interpolator2D & iImT, Interpolator2D & iImG);
double E_integral(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);
double E_integral_cm(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);
double E_int_first(double om, double p, double T, double k, Interpolator2D & iImT, Interpolator2D & iImG);
double Res_Kint(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);
double x_integral_cm(double E, double om, double p, double k, Interpolator2D & iImT, Interpolator2D & iImG);
double func_E_int_first(double e, double om, double p, double T, double k, Interpolator2D & iImT, Interpolator2D & iImG);

void get_E_int(double om, double T, Interpolator2D & iImT, Interpolator2D & iImG, double * p, int dimP, double * out, int dimOut);
// void get_T(double E, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, double * p, int dimP, 
//             double *out, int dimOut, double * out2, int dimOut2);
void get_T(double E, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, double * p, int dimP, 
            // double *out, int dimOut);
            std::complex<double> *out, int dimOut);
std::complex<double> T_solve(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda = 5, int sign=1);
std::complex<double> T_solve_den(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda = 5);
std::complex<double> T_solve_test(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            IntGSL<funct> integ_T, double Lambda = 5);

// double T_solveRe(double E, double q, double q1, double T, Interpolator iVK, Interpolator iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
//             double Lambda, IntGSL<funct> integ);

double T_solveRe(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda);

double ReSigmaKK_2D(double E, double q, Interpolator2D & iImS);
double ReSigmaKK(double E, Interpolator & iImS);

double sigma_bb(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);
double sigma_bf(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);
double sigma_ff(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);
double sigma_fb(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);

void get_test_gsl_interp(int N, double * p, int dimP, double * out, int dimOut);
//////////////////////////////////////////////////////////////////////

// Another integration variable choice
double sigma_bb2(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);
double sigma_bb3(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG);

double k_integral_cm2(double omp, double om, double p, Interpolator2D & iImT, 
    Interpolator2D & iImG);
double x_integral_cm2(double omp, double om, double p, double k, Interpolator2D & iImT, 
Interpolator2D & iImG, int debug=0);


void get_test(double * p, int dimP, double * out, int dimOut);
void get_test_gsl(int N, double * p, int dimP, double * out, int dimOut);



std::complex<double> T_solve_BB(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda = 5, int sign=1);
std::complex<double> T_solve_BF(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda = 5, int sign=1);

double OmS2_F_om_int2(double omp, double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);
double OmS2_F_om_int1(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);
double OmS2_F(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);

std::complex<double> T_solve_explicit(
    double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign);

std::complex<double> J_solve_explicit(
    double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign);

double OmS_B_qfirst_q_int(double om, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                            Interpolator2D & iImS, Interpolator2D & iReS);
                            
double OmS_B_qfirst(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);

double n_f(double om, double T);
double n_b(double om, double T);
// class Runner {
//     public:
//         Runner();
        
//         Interpolator2D & iImS;
//         Interpolator2D & iGq;
//         Interpolator2D & iImGqq;
//         Interpolator2D & iReGqq;
// };

double test_integration(double a);

double k_integral_QQ(double E, double om, double p, Interpolator2D & iImG);
double E_integral_QQ(double om, double p, double T, Interpolator2D & iImG);
double k_integral_QQ_func(double k, double E, double om, double p, Interpolator2D & iImG);

double sigma_integrand_bb(double omp, double om, double p, double T, Interpolator2D & iImT, 
                        Interpolator2D & iImG);

double sigma_integrand_bf(double omp, double om, double p, double T, Interpolator2D & iImT, 
                        Interpolator2D & iImG);

double sigma_integrand_fb(double omp, double om, double p, double T, Interpolator2D & iImT, 
                        Interpolator2D & iImG);

////////////////////////// Phi-functional stuff ////////////////////////////////




//////////////////////////// On-shell integrals ////////////////////////////////






#endif