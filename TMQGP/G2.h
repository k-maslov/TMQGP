

#ifndef G2_H_
#define G2_H_

#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>
#include "SigmaInt.h"
#include <string>
#include <vector>
using namespace std;

double G2_conv_ff(double om, double q, double T, Interpolator2D & R1, Interpolator2D & R2, double Lambda=5);
double ReG2_conv_ff_integrand_subtr(double om1, double om, double q, double T, PoleInterpolator & R1, PoleInterpolator & R2, double Lambda);
double ReG2_conv_ff_integrand(double om1, double om, double q, double T, Interpolator2D & R1, Interpolator2D & R2, double Lambda=5);
double ReG2_conv_ff(double om, double q, double T, Interpolator2D & R1, Interpolator2D & R2, double Lambda=5);
double ReG2_pole(double om, double q, double T, PoleInterpolator & R1, PoleInterpolator & R2, double Lambda=5);
double ReG2_subtr(double om, double q, double T, PoleInterpolator & R1, PoleInterpolator & R2, double Lambda);

#endif