

#ifndef SIGMAPROD_H_
#define SIGMAPROD_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>
#include "SigmaInt.h"
#include <string>
#include <vector>
using namespace std;

class TMChannel{
    public:
        TMChannel(){};
        TMChannel(double da, double ds, double d, double Nf, string stat_i, string stat_j, Interpolator2D * iImT, Interpolator * eps_i, Interpolator * eps_j){
            this->da = da;
            this->ds = ds;
            this->d = d;
            this->Nf = Nf;
            this->stat_i = stat_i;
            this->stat_j = stat_j;
            this->iImT = iImT;
            this->eps_i = eps_i;
            this->eps_j = eps_j;
        }
        double da;
        double ds;
        double d;
        double Nf;
        string stat_i;
        string stat_j;
        Interpolator2D * iImT;
        Interpolator * eps_i;
        Interpolator * eps_j;
};



typedef vector<TMChannel> TMArray;

double x_cm_onshell_integrand(double x, double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double x_integral_cm_onshell(double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int l=0, int debug=0);

double k_integral_onshell(double omp, double om, double p, 
Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
int l=0);


double sigma_bf_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l=0);

double sigma_fb_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l=0);

double sigma_ff_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
    int l=0);

double sigma_bb_onshell(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l=0);

double sigma_tot(double om, double p, double T, TMArray TMs, Interpolator2D & iImG_Q, Interpolator2D & iImG_G);


double x_cm_onshell_integrand2(double x, double omp, double om, double p, 
            double k, double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double x_integral_cm_onshell2(double omp, double om, double p, 
            double k, double T, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double k_integral_onshell2(double omp, double om, double p, double T,
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
    int l);

double sigma_ff_onshell2(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
    int l);

double x_cm_onshell_integrand_allom(double x, double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double x_integral_cm_onshell_allom(double omp, double om, double p, 
            double k, Interpolator2D & iImT, Interpolator2D & iImG, 
            Interpolator & eps1, Interpolator & eps2, int debug, int l);

double k_integral_onshell_allom(double omp, double om, double p, 
Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
int l);

double sigma_ff_onshell_allom(double om, double p, double T, 
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2, 
    int l);




#endif