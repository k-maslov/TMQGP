#include "mp.h"
#include "SigmaProd.h"
#include <gsl/gsl_errno.h>
#include <complex>
// #include <omp>
#include <omp.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_legendre.h>
#include <mpi.h>

void get_sigma_ff(double * omrange, int dimOmrange, double p, double T,
    Interpolator2D & iImT, Interpolator2D & iImG, Interpolator & eps1, Interpolator & eps2,
    int l, double * out, int dimOut){
    int i = 0;
    
    omp_set_dynamic(1);
    omp_set_num_threads(14);

    #pragma omp parallel for shared(out, iImT, iImG, eps1, eps2) private(i)
    for (i = 0; i < dimOmrange; i++){
        out[i] = sigma_ff_onshell(omrange[i], p, T, iImT, iImG, eps1, eps2);
    }

    // return 0;
}

double get_integ(double a, Int_gsl_adaptive integ){

    funct f = [&](double x) -> double{
            return 1. / (x*x + a*a);
        };

    double res, err;
    integ.integrate(&f, 0, 100, res, err);
    // return a*a;
}

void try_integration(double * omrange, int dimOmrange, int N, double * out, int dimOut){
     
    
    // double * res = new double[dimOmrange];
    int i = 0;

    Int_gsl_adaptive integ;

    // omp_set_dynamic(1);
    omp_set_num_threads(N);

    #pragma omp parallel for shared(out) private(i, integ)
    for (i = 0; i < dimOmrange; i++){
        out[i] = get_integ(omrange[i], integ);
    }
}


void mpi_hw(){
    MPI_Init(NULL, NULL);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Size %i, rank %i", size, rank);

    MPI_Finalize();
}