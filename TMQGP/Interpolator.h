#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H

#include "gsl/gsl_spline.h"
#include "gsl/gsl_spline2d.h"
#include <map>
#include <string>
#include <vector>

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

	vector<double> x;
	vector<double> data;
	string kind;

	double operator()(double x);
	double D(double x) {
		return gsl_spline_eval_deriv(interp, x, accel);
	}

	// template <class Archive> void  serialize(Archive & archive);
};

class Interpolator2D{
	public:
		Interpolator2D(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2);

		gsl_spline2d * interp;
		gsl_interp_accel * accX;
		gsl_interp_accel * accY;
		bool debug;

		vector<double> x;
		vector<double> y;
		vector<double> z;

		double operator()(double x, double y);
};

// class InterDenom2D{
// 	public:
// 		InterDenom2D(){};
// 		InterDenom2D(double *x, int dimX, double *y, int dimY, 
// 			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4);

// 		double real(double x, double y);
// 		double imag(double x, double y);
// 		complex<double> operator()(double x, double y);
// };
// ########## Write the integrations here temporarily




#endif