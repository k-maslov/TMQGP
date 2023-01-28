#pragma once

#include "gsl/gsl_integration.h"
#include "Integrator.h"

namespace INT_GSL {
	template <class C> double f(double x, void * params) {
		C * fun = (C *)params;
		return (*fun)(x);
	}
}

template <class C> class IntGSL :
public Integrator<C>
{
public:
	static gsl_integration_workspace * w;

	IntGSL() {
// 		w = nullptr;
// 		init_workspace();
	}
	
	IntGSL(C* c) : Integrator<C>(c){
// 		w = nullptr;
// 		init_workspace();
	}

	void init_workspace() {
		if (w){
			delete w;
		}
		w = gsl_integration_workspace_alloc(1e6);
	}

	virtual void integrate(C * func, double a, double b, double & result, double & error) {
		gsl_function gsl_f;
		gsl_f.function = &INT_GSL::f<C>;
		gsl_f.params = func;
		size_t neval;
		if (gsl_finite(a) && gsl_finite(b))
			// gsl_integration_qags(&gsl_f, a, b, 1e-7, 1e-7, 1000,
			// 	w, &result, &error);
			gsl_integration_qng(&gsl_f, a, b, 1e-5, 1e-5, &result, &error, &neval);
		else {
			if (!gsl_finite(a) && !gsl_finite(b)) {
				gsl_integration_qagi(&gsl_f, 1e-7, 1e-7, 1000,
					w, &result, &error);
			}
			else if (!gsl_finite(a)) {
				gsl_integration_qagil(&gsl_f, b, 1e-7, 1e-7, 1000,
					w, &result, &error);
			}
			else {
				gsl_integration_qagiu(&gsl_f, a, 1e-7, 1e-7, 1000,
					w, &result, &error);
			}
		}
	}
};

template <class C> gsl_integration_workspace * IntGSL<C>::w = gsl_integration_workspace_alloc(1e6);

class Int_gsl_adaptive : public IntGSL<funct> {
	public:
		virtual void integrate(funct * func, double a, double b, double & result, double & error) override;
};

class Int_gsl_cauchy: public IntGSL<funct>{
	public:
		void integrate(funct * func, double a, double b, double wvar, double & result, double & error);
};

class Int_gsl_fixed : public IntGSL<funct> {
	public:
		double a, b;
		gsl_integration_fixed_workspace * wfix;
		Int_gsl_fixed(double a, double b);
		virtual void integrate(funct * func, double a, double b, double & result, double & error) override;
};