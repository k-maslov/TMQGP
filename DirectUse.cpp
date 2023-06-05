#include <iostream>
#include <fstream>
#include <cstdio>
#include "SigmaInt.h"
#include <gsl/gsl_matrix.h>
#include "Interpolator.h"
#include "gsl/gsl_vector.h"
#include <complex>

#include <chrono>
#include<omp.h>

using namespace std::chrono;

using namespace std;


void test_vec_interp2d(){

}

int main(){
    // ifstream fImG2;
    FILE * fImG2 = fopen("ImG2.dat", "r");
    FILE * fQrange = fopen("qrange.dat", "r");
    FILE * fErange = fopen("erange.dat", "r");
    FILE * fReG2 = fopen("ReG2.dat", "r");
    FILE * fOmK = fopen("OmK.dat", "r");
    FILE * fV = fopen("V.dat", "r");


    // gsl_matrix * m = gsl_matrix_alloc(51, 201);
    gsl_matrix * mImG2 = gsl_matrix_alloc(201, 51);
    gsl_matrix * mReG2 = gsl_matrix_alloc(201, 51);
    gsl_vector * qrange = gsl_vector_alloc(51);
    gsl_vector * v = gsl_vector_alloc(51);
    gsl_vector * omk = gsl_vector_alloc(51);
    gsl_vector * erange = gsl_vector_alloc(201);
    // gsl_matrix * m = gsl_matrix_alloc(3, 3);
    int res = gsl_matrix_fscanf(fImG2, mImG2);
    gsl_matrix_fscanf(fReG2, mReG2);
    gsl_vector_fscanf(fQrange, qrange);
    gsl_vector_fscanf(fErange, erange);
    gsl_vector_fscanf(fOmK, omk);
    gsl_vector_fscanf(fV, v);


    Interpolator2D iImG2(qrange->data, 51, erange->data, 201, mImG2->data, 51, 201);
    Interpolator2D iReG2(qrange->data, 51, erange->data, 201, mReG2->data, 51, 201);
    // Interpolator2D i(qrange->data, 51, erange->data, 201, mImG2->data, 51, 201);

    Interpolator iV(qrange->data, 51, v->data, 51, "linear");
    Interpolator iOmK(qrange->data, 51, omk->data, 51, "linear");
    // cout << inter(1, 2) << endl;

    // cout << res << endl;

    // cout << T_solve(0.5, 0.5, 0.5, 0.2, iV, iOmK, iReG2, iImG2) << endl;

    int N = 150000;
    

    gsl_vector * qrange_test = gsl_vector_alloc(N);

    for (int i = 0; i < N; i++){
        qrange_test->data[i] = 4.5/N * i;
    }

    double q;
    auto start = high_resolution_clock::now();

    for (int i = 0; i < N; i++){
        q = qrange_test->data[i];
        T_solve(0.5, q, q, 0.2, iV, iOmK, iReG2, iImG2);
    }

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    // // cout << "time serial : " << duration.count() << endl;

    // // std::complex<double> * out = new complex<double> [N];
    // // std::complex<double> * out = new complex<double> [N];
    // double * out = new double[N];


    // double * qrange_arr = new double[N];

    // for (int i = 0; i < N; i++){
    //     qrange_arr[i] = gsl_vector_get(qrange_test, i);
    // }

    // int Nthr = 0;

    // for (Nthr = 1; Nthr < 18; Nthr+=4)
    // {
    //     start = high_resolution_clock::now();
    //     auto start_omp = omp_get_wtime();

    //     // get_T(0.5, 0.2, iV, iOmK, iReG2, iImG2, qrange_arr, N, out, N);
    //     // get_test_gsl_interp(Nthr, qrange_arr, N, out, N);

    //     stop = high_resolution_clock::now();
    //     auto stop_omp = omp_get_wtime();

    duration = duration_cast<microseconds>(stop - start);
    cout << "time = " << " " << duration.count() << endl;//<< "   " << stop_omp - start_omp << endl;
    



    
//// Check whether the matrix was read correctly

    // for (int i = 0; i < 5; i++){
    //     for (int j = 0; j < 10; j++){
    //         // cout << m->data[i, j] << "  ";
    //         cout << gsl_matrix_get(m, i, j) << "  ";
    //     }
    //     cout << endl;
    // }

    // Interpolator2D

}