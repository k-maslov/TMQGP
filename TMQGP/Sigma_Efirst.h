#ifndef SIGMA_Efirst_H_
#define SIGMA_Efirst_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>
#include "SigmaInt.h"
#include <string>
#include <vector>
using namespace std;



double Efirst_cm_onshell_integrand(double omp, double om, double p, 
            double k, double x,  double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double Efirst_integral_cm_onshell(double om, double p, 
            double k, double x, double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double Efirst_x_integral(double om, double p, 
            double k, double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double Efirst_k_integral(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug=0, int l=0);

double Efirst_kfirst_int(double om, double p, 
            double x, double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double Efirst_kfirst_x(double om, double p, double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

#endif