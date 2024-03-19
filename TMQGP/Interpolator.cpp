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
	gsl_spline_free(interp);
	gsl_interp_accel_free(accel);
}

double Interpolator::operator()(double x) {
// 	cout << interp->x[0] << " ; " << interp->x[interp->size - 1] << endl;
// 	cout << interp->y[0] << " ; " << interp->y[interp->size - 1] << endl;
	return gsl_spline_eval(interp, x, accel);
}


Interpolator2D::Interpolator2D(double *x, int dimX, double *y, int dimY, 
			double * z2, int dimZ1, int dimZ2){
	interp = gsl_spline2d_alloc(gsl_interp2d_bicubic, dimX, dimY);
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, z2, dimX, dimY);

	this->x.clear();
	this->y.clear();
	this->z.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(z2[i]);
	}

	debug = 0;
}

Interpolator2D::~Interpolator2D()
{
	gsl_spline2d_free(this->interp);
	gsl_interp_accel_free(this->accX);
	gsl_interp_accel_free(this->accY);


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

InterDenom2D::InterDenom2D(double *x, int dimX, double *y, int dimY, double *ReZ2, 
			int dimZ1, int dimZ2, double *ImZ2, int dimZ3, int dimZ4, string what){
	this->interp = gsl_spline2d_alloc(gsl_interp2d_bicubic, dimX, dimY);
	this->iIm = gsl_spline2d_alloc(gsl_interp2d_bicubic, dimX, dimY);
	accImX = gsl_interp_accel_alloc();
	accImY = gsl_interp_accel_alloc();
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, ReZ2, dimX, dimY);
	gsl_spline2d_init(iIm, x, y, ImZ2, dimX, dimY);

	this->what = what;

	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->z2.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(ReZ2[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z2.push_back(ImZ2[i]);
	}

	debug = 0;
}

InterDenom2D::~InterDenom2D()
{
	// gsl_spline2d_free(this->interp);
	// gsl_interp_accel_free(this->accX);
	// gsl_interp_accel_free(this->accY);
	gsl_spline2d_free(this->iIm);
	gsl_interp_accel_free(this->accImX);
	gsl_interp_accel_free(this->accImY);
}

double InterDenom2D::real(double x, double y){
	double re = gsl_spline2d_eval(interp, x, y, accX, accY);
	double im = gsl_spline2d_eval(iIm, x, y, accImX, accImY);
	return re / (re*re + im*im);
}

double InterDenom2D::imag(double x, double y){
	double re = gsl_spline2d_eval(interp, x, y, accX, accY);
	// cout << "re = " << re << endl;
	double im = gsl_spline2d_eval(iIm, x, y, accImX, accImY);
	// cout << "im = " << im << endl;
	double res = - im / (re*re + im*im);
	if (gsl_isnan(res)){
		return 0;
	}
	return res;
}

double InterDenom2D::operator()(double x, double y){
	if (this->what == "real"){
		return this->real(x, y);
	}
	if (this->what == "imag"){
		return this->imag(x, y);
	}
	else{
		return -1;
	}
}

PoleInterpolator::PoleInterpolator(double *x, int dimX, double *y, int dimY,
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4,
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what) : 
			InterDenom2D(x, dimX, y, dimY, ReZ2, dimZ1, dimZ2, ImZ2, dimZ3, dimZ4, what)
			{
	iPole = new Interpolator(q, dimQ, pole, dimPole, "cubic");
	iWidth = new Interpolator(q, dimQ, width, dimWidth, "cubic");
}

PoleInterpolator::~PoleInterpolator(){
	delete iPole;
	delete iWidth;
}

GFInterpolator::GFInterpolator(double *x, int dimX, double *y, int dimY, 
			double * ReZ2, int dimZ1, int dimZ2, double * ImZ2, int dimZ3, int dimZ4, 
			double * q, int dimQ, double * pole, int dimPole, double * width, int dimWidth, string what, double m)
{
	this->m = m;
	this->interp = gsl_spline2d_alloc(gsl_interp2d_bicubic, dimX, dimY);
	this->iIm = gsl_spline2d_alloc(gsl_interp2d_bicubic, dimX, dimY);
	accImX = gsl_interp_accel_alloc();
	accImY = gsl_interp_accel_alloc();
	accX = gsl_interp_accel_alloc();
	accY = gsl_interp_accel_alloc();
	gsl_spline2d_init(interp, x, y, ReZ2, dimX, dimY);
	gsl_spline2d_init(iIm, x, y, ImZ2, dimX, dimY);

	this->what = what;

	this->x.clear();
	this->y.clear();
	this->z.clear();
	this->z2.clear();

	for (int i = 0; i < dimX; i++){
		this->x.push_back(x[i]);
	}

	for (int i = 0; i < dimY; i++){
		this->y.push_back(y[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z.push_back(ReZ2[i]);
	}

	for (int i = 0; i < dimY*dimX; i++){
		this->z2.push_back(ImZ2[i]);
	}

	debug = 0;	
	iPole = new Interpolator(q, dimQ, pole, dimPole, "cubic");
	iWidth = new Interpolator(q, dimQ, width, dimWidth, "cubic");
}

GFInterpolator::~GFInterpolator()
{
	
}

double GFInterpolator::real(double x, double y)
{
	return gsl_spline2d_eval(interp, x, y, accX, accY);
}

double GFInterpolator::imag(double x, double y)
{
	return gsl_spline2d_eval(iIm, x, y, accImX, accImY);
}

double GFInterpolator::operator()(double x, double y)
{
	double eps = sqrt(x*x + this->m*this->m);
	double res;
	complex<double> G = 1./(y - eps - complex<double>(this->real(x, y), this->imag(x, y)));
	if (this->what == "real"){
		res = G.real();
	}
	else if (this->what == "imag"){
		res = G.imag();
	}
	else{
		res = -1;
	}
	if (gsl_isnan(res)){
		return 0;
	}
	return res;
}

