#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H

#include "gsl/gsl_spline.h"
#include "gsl/gsl_spline2d.h"
#include <map>
#include <string>

using namespace std;

class Interpolator {
public:
	Interpolator(double * x, int dimX, double * y, int dimY, string kind);

	Interpolator();
// 	Interpolator(Interpolator &) = delete;
	~Interpolator();

	gsl_spline * interp;
	gsl_interp_accel * accel;
	map<string, const gsl_interp_type *> methods;

	double operator()(double x);
	double D(double x) {
		return gsl_spline_eval_deriv(interp, x, accel);
	}
};

class Interpolator2D{
	public:
		Interpolator2D(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2);

		gsl_spline2d * interp;
		gsl_interp_accel * accX;
		gsl_interp_accel * accY;
		bool debug;

		double operator()(double x, double y);
};


// ########## Write the integrations here temporarily




#endif