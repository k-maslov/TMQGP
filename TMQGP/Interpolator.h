#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H

#include "gsl/gsl_spline.h"
#include "gsl/gsl_spline2d.h"
#include <map>
#include <string>
#include <vector>
#include <complex>

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
		Interpolator2D(){};
		Interpolator2D(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2);

		~Interpolator2D();

		gsl_spline2d * interp;
		gsl_interp_accel * accX;
		gsl_interp_accel * accY;
		bool debug;

		vector<double> x;
		vector<double> y;
		vector<double> z;

		virtual double operator()(double x, double y);
};

class InterDenom2D : public Interpolator2D{
	public:
		InterDenom2D(){};
		InterDenom2D(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, 
			double * ImZ2, int dimZ3, int dimZ4, string what);

		~InterDenom2D();

		// Assume base class members interpolate the real part

		// gsl_interp_accel * accReX;
		// gsl_interp_accel * accReY;
		gsl_interp_accel * accImX;
		gsl_interp_accel * accImY;

		// gsl_spline2d * iRe;
		gsl_spline2d * iIm;

		vector<double> z2;
		string what;
		virtual double real(double x, double y);
		virtual double imag(double x, double y);
		double operator()(double x, double y) override;
};
// ########## Write the integrations here temporarily


class PoleInterpolator: public InterDenom2D{
	public:
		PoleInterpolator(){};
		PoleInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what);

		Interpolator * iPole;
		Interpolator * iWidth;

		~PoleInterpolator();
};

// Interpolating Sigma and return G on evaluation
class GFInterpolator: public PoleInterpolator{
	public:
		GFInterpolator(){};
		GFInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m);

		~GFInterpolator();

		double m;
		double real(double x, double y) override;
		double imag(double x, double y) override;
		double operator()(double x, double y) override;
};


class RhoInterpolator: public GFInterpolator{
	public:
		RhoInterpolator(){};
		RhoInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m);

		double operator()(double x, double y) override;
};

#endif