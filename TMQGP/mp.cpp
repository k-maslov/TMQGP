#include "mp.h"
#include "SigmaProd.h"
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>

void get_sigma_ff(double * omrange, int dimOmrange, double p, double T,
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l, double * out, int dimOut){
    int i = 0;
    
    omp_set_dynamic(1);
    omp_set_num_threads(14);

    #pragma omp parallel for shared(out) private(i)
    for (i = 0; i < dimOmrange; i++){
        out[i] = sigma_ff_onshell(omrange[i], p, T, iImT, iImG, eps1, eps2);
    }

    // return 0;
}