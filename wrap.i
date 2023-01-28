%module TMQGP

%include <std_map.i>
%include <std_string.i>
%include <std_unordered_map.i>
%include <std_vector.i>


%{
#define SWIG_FILE_WITH_INIT
#define PY_ARRAY_UNIQUE_SYMBOL PyArray_API
%}
%include "numpy.i"
%init %{
import_array();
%}

%feature("autodoc", 1);

typedef std::string String;
%template(vec_type) std::vector<double>;
%template(dval) std::map<std::string, double>;
%template(content_type) std::map<std::string, std::map<std::string, double>>;
%template(str_vector) std::vector<std::string>;
%template(suppress_type) std::map<std::string, int>;
%template(eos_type) std::map<std::string, std::vector<double>>;
%template(mr_type) std::vector<std::map<std::string, double>>;
%template(integfunc_type) std::map<std::string, std::function<double (double, const double[])>>;

//%template(outfunc_type) std::map<std::string, std::vector<double>>;


%naturalvar Solver::variables;
%naturalvar Model::active_fields;
%naturalvar Model::charged_fields;
%naturalvar TOVsolver::iEoS;
%naturalvar TOVsolver::variables;

%{

//#include "stacktrace/stack_exception.hpp"
//#include "stacktrace/call_stack.hpp"
#include <functional>
#include <map>

#include "TMQGP/Interpolator.h"
#include "TMQGP/SigmaInt.h"
// #include <complex>

%}


%include "std_function.i";
%include "std_complex.i"

%std_function(Functor, double, double);

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* x, int dimX)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* y, int dimY)};
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double* z2, int dimZ1, int dimZ2)};
%include "TMQGP/Interpolator.h"

%apply (double* INPLACE_ARRAY1, int DIM1) {(double * p, int dimP)};
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double * out, int dimOut)};
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double * out2, int dimOut2)};
%apply (complex<double>* ARGOUT_ARRAY1, int DIM1) {(complex<double> * out, int dimOut)};
%include "TMQGP/SigmaInt.h"

