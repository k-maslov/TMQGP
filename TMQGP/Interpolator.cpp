#include "Interpolator.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
using namespace std;


Interpolator::Interpolator() {

}

Interpolator::Interpolator(double * x, int dimX, double * y, int dimY, string kind) {
	methods = { { "linear", gsl_interp_linear },{ "cubic", gsl_interp_cspline },
	{ "steffen", gsl_interp_steffen } };
	if (dimX != dimY) {
		throw;
	}
	this->kind = kind;
	if (methods.find(kind) != methods.end()) {
		interp = gsl_spline_alloc(methods[kind], dimX);
		accel = gsl_interp_accel_alloc();
		gsl_spline_init(interp, x, y, dimX);
	}
	else {
		throw;
	}

	this->x.clear();
	this->data.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
		this->data.push_back(y[i]);
	}

}

Interpolator::~Interpolator() {
// 	gsl_spline_free(interp);
// 	gsl_interp_accel_free(accel);
}

double Interpolator::operator()(double x) {
// 	cout << interp->x[0] << " ; " << interp->x[interp->size - 1] << endl;
// 	cout << interp->y[0] << " ; " << interp->y[interp->size - 1] << endl;
	return gsl_spline_eval(interp, x, accel);
}


Interpolator2D::Interpolator2D(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2){
	interp = gsl_spline2d_alloc(gsl_interp2d_bilinear, dimX, dimY);
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, z2, dimX, dimY);

	this->x.clear();
	this->y.clear();
	this->z.clear();

	// for (int i = 0; i < dimX; i++){
	// 	this->x.push_back(x[i]);
	// }

	// for (int i = 0; i < dimY; i++){
	// 	this->y.push_back(y[i]);
	// }

	// for (int i = 0; i < dimY*dimX; i++){
	// 	this->z.push_back(z2[i]);
	// }

	debug = 0;
}

double Interpolator2D::operator()(double x, double y){
	// gsl_set_error_handler_off();
	double res;
	// try
	// {
	// 	res = gsl_spline2d_eval(interp, x, y, accX, accY);
	// }
	// catch(const std::exception& e)
	// {
	// 	std::cerr << e.what() << '\n';
	// 	res = 0;
	// }
	res = gsl_spline2d_eval(interp, x, y, accX, accY);
	if (gsl_isnan(res)){
		return 0;
	}

	if (debug){
		printf("Interp2D: %.3e %.3e %.3e", x, y, res);
	}
	return res;
}