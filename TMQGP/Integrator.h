#pragma once

#include <functional>
#include <map>
#include <string>
// #include "Terms/TermInput.h"
typedef std::map<std::string, double> dval;

typedef std::function<double(double)> funct;
typedef std::map<std::string, double> dval;
typedef std::function<double (double, const double[])> odefunct;
typedef std::map<std::string, odefunct> funct_arr;
typedef std::pair<std::string, double> elem;

template <class C> class Integrator
{
public:
	Integrator() {

	}
	Integrator(C * c) {
		this->func = c;
	}
	virtual void integrate(double a, double b, double & result, double & error) {
		return integrate(this->func, a, b, result, error);
	}

	virtual void integrate(C * func, double a, double b, double & result, double & error) = 0;
	C * func;

	dval params;
	dval full_output;

	virtual ~Integrator() = default;
};

