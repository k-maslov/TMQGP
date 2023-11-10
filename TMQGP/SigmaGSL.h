#ifndef SIGMAGSL_H_
#define SIGMAGSL_H_


#include <complex>
#include <string>
#include <vector>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>

using namespace std;

class SigmaWorkspace{
    SigmaWorkspace(
        double *qrange, int dimQ,
        double *erange, int dimE, 
		double * T, int dimT1, int dimT2,
        double * R, int dimR1, int dimR2,
        double m1, double m2){
            
    }
};

#endif