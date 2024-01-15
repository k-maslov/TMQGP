#ifndef TMATRIX_H_
#define TMATRIX_H_
#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>

std::complex<double> x_solve(double E, double q, double q1, double T, Interpolator & iVK, Interpolator & iOmK, Interpolator2D & iReGqq, Interpolator2D & iImGqq, 
            double Lambda, int sign, int adaptive=0);



#endif