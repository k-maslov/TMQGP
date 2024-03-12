#ifndef SIGMAQP_H_
#define SIGMAQP_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>
#include "SigmaInt.h"
#include <string>
#include <vector>
using namespace std;

double QM_integrand(double om, double p, 
            double k, double x, double T, Interpolator2D & iImT, 
            Interpolator & eps1, Interpolator & eps2, int debug=0, int l=0);

double QP_x_int(double om, double p,
            double k, double T, Interpolator2D & iImT, 
            Interpolator & eps1, Interpolator & eps2, int debug=0, int l=0);

double SigmaQP(double om, double p,
            double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l, double UL=5);

#endif