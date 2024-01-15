#include "IntGSL.h"
#include "gsl/gsl_integration.h"
#include "Integrator.h"

void Int_gsl_adaptive::integrate(funct * func, double a, double b, double & result, double & error){
    	gsl_function gsl_f;
		gsl_f.function = &INT_GSL::f<funct>;
		gsl_f.params = func;
		size_t neval;

        gsl_integration_qag(&gsl_f, a, b, 1e-7, 1e-7, 20, 2,
				w, &result, &error);
}

void Int_gsl_cauchy::integrate(funct * func, double a, double b, double wvar, double & result, double & error){
		gsl_function gsl_f;
		gsl_f.function = &INT_GSL::f<funct>;
		gsl_f.params = func;
		size_t neval;

        gsl_integration_qawc(&gsl_f, a, b, wvar, 1e-7, 1e-7, 1000, w, &result, &error);
}

Int_gsl_fixed::Int_gsl_fixed(double a, double b){
	const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
	wfix = gsl_integration_fixed_alloc(T, 150, a, b, 0, 0);
	this->a = a;
	this->b = b;
}

void Int_gsl_fixed::integrate(funct * func, double a, double b, double & result, double & error){
    	gsl_function gsl_f;
		gsl_f.function = &INT_GSL::f<funct>;
		gsl_f.params = func;
		size_t neval;
	// 				const gsl_integration_fixed_type * T = gsl_integration_fixed_legendre;
	// wfix = gsl_integration_fixed_alloc(T, 50, a, b, 0, 0);

        gsl_integration_fixed(&gsl_f, &result, wfix);
}