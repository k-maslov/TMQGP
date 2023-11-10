#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <gsl/gsl_spline2d.h>
#include "Interpolator.h"
class Particle{
    double * erange;
    double * qrange;
    int dimE, dimQ;
    double m;

    Interpolator2D * iR;
    Interpolator2D * iReG;
    Interpolator2D * iImG;

    Particle(double m, double * erange, int dimE, double * qrange, int dimQ, bool stat, double eps=2e-2, double d=6);

    double R(double q, double e);

    double om0(double q);
    double G0(double q, double e);
        
};


#endif