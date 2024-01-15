#include <mpi.h>
#include <cstdio>
#include "TMQGP/SigmaProd.h"
#include <gsl/gsl_matrix.h>
#include "Interpolator.h"
#include "gsl/gsl_vector.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

int main(int argc, char *argv[]){
    // define and load interpolators
    FILE * fImG = fopen("testmpi/ImG.dat", "r");
    FILE * fImT = fopen("testmpi/ImT.dat", "r");
    FILE * fQrange = fopen("testmpi/qrange.dat", "r");
    FILE * fErange = fopen("testmpi/erange.dat", "r");
    FILE * fReG = fopen("testmpi/ReG.dat", "r");
    FILE * fOmk = fopen("testmpi/omk.dat", "r");

    gsl_matrix * mImG = gsl_matrix_alloc(201, 51);
    gsl_matrix * mReG = gsl_matrix_alloc(201, 51);
    gsl_matrix * mImT = gsl_matrix_alloc(201, 51);
    gsl_vector * qrange = gsl_vector_alloc(51);
    gsl_vector * erange = gsl_vector_alloc(201);
    gsl_vector * omk = gsl_vector_alloc(51);
    // gsl_matrix * m = gsl_matrix_alloc(3, 3);

    gsl_matrix_fscanf(fImG, mImG);
    gsl_matrix_fscanf(fReG, mReG);
    gsl_matrix_fscanf(fImT, mImT);
    gsl_vector_fscanf(fQrange, qrange);
    gsl_vector_fscanf(fErange, erange);
    gsl_vector_fscanf(fOmk, omk);
    
    Interpolator2D iImG(qrange->data, 51, erange->data, 201, mImG->data, 51, 201);
    Interpolator2D iReG(qrange->data, 51, erange->data, 201, mReG->data, 51, 201);
    Interpolator2D iImT(qrange->data, 51, erange->data, 201, mImT->data, 51, 201);
    Interpolator iEps(qrange->data, 51, omk->data, 51, "cubic");

    std::complex<double> * pairs = new std::complex<double>[201*51];
    for (int i = 0; i < 51; i++){
        for (int j = 0; j < 201; j++){
            pairs[201*i + j] = {qrange->data[i], erange->data[j]};
        }
    }

    // cout << iImG(0.5, 0.5) << endl;

    double res = sigma_ff_onshell(0.5, 0.5, 0.2, iImT, iImG, iEps, iEps);

    cout << res << endl;
    auto start = high_resolution_clock::now();
    // return 0;


    int rank, size;     // for storing this process' rank, and the number of processes
    int *sendcounts;    // array describing how many elements to send to each process
    int *displs;        // array describing the displacements where each segment begins

    int sum = 0;                // Sum of counts. Used to calculate displacements
    std::complex<double> rec_buf[201*51];          // buffer where the received data should be stored

    int SIZE = 201*51;

    // the data to be distributed
    // double * data = new double[SIZE];

    double * output = new double[SIZE];

    // for (int i = 0; i < SIZE; i++){
    //     data[i] = i;
    // }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rem = (SIZE)%size; // elements remaining after division among processes

    sendcounts = new int[size];//malloc(sizeof(int)*size);
    displs = new int[size];//malloc(sizeof(int)*size);

    // calculate send counts and displacements
    for (int i = 0; i < size; i++) {
        sendcounts[i] = (SIZE)/size;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }

    // print calculated send counts and displacements for each process
    if (0 == rank) {
        printf("rem = %i \n", rem);
        for (int i = 0; i < size; i++) {
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }
    }

    // divide the data among processes as described by sendcounts and displs
    MPI_Scatterv(pairs, sendcounts, displs, MPI_CXX_DOUBLE_COMPLEX, &rec_buf, 201*51, MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // print what each process received
    printf("%d: ", rank);
    for (int i = 0; i < sendcounts[rank]; i++) {
        printf("%f\t", rec_buf[i]);
    }
    printf("\n");


    double * out = new double[sendcounts[rank]];

    for (int i = 0; i < sendcounts[rank]; i++){
        // out[i] = rec_buf[i] * rec_buf[i];
        double q = rec_buf[i].real();
        double e = rec_buf[i].imag();
        out[i] = sigma_ff_onshell(e, q, 0.2, iImT, iImG, iEps, iEps);
    }

    printf("Output of %d: ", rank);
    for (int i = 0; i < sendcounts[rank]; i++) {
        printf("%f\t", out[i]);
    }
    printf("\n");

    MPI_Gatherv(out, sendcounts[rank], MPI_DOUBLE, output, sendcounts, displs, 
        MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0){
        printf("Received data: \n");

        for (int i = 0; i < SIZE; i++){
            printf("%f\t", output[i]);
        }
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(stop - start);

        cout << "time = " << " " << duration.count() << endl;//<< "   " << stop_omp - start_omp << endl;
    }

    MPI_Finalize();


   

    // free(sendcounts);
    // free(displs);

    return 0;


    return 0;
}