
#ifndef SIGMAPROD_H_
#define SIGMAPROD_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>
#include "SigmaInt.h"


double x_cm_onshell_integrand(double x, double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator eps1, Interpolator eps2, int debug, int l);

double x_integral_cm_onshell(double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator eps1, Interpolator eps2, int l=0, int debug=0);

double k_integral_onshell(double omp, double om, double p, 
Interpolator2D & iImT, Interpolator2D & iImG, Interpolator eps1, Interpolator eps2, 
int l=0);


double sigma_bf_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator eps1, Interpolator eps2,
    int l=0);

double sigma_fb_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator eps1, Interpolator eps2,
    int l=0);

double sigma_ff_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator eps1, Interpolator eps2, 
    int l=0);

double sigma_bb_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator eps1, Interpolator eps2,
    int l=0);



#endif