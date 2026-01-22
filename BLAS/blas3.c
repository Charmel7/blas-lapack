#include <stdio.h>
#include <mkl.h>
/*

In this file we're exploring the Basic Linear Algebra Subprogram Level 3 for operations between matrices.
We'll do some current examples of use and explain how it works in the readme.
Following are just tests based on cblas, the C interface of the native BLAS in Fortran77.

Author: Prosper Affoukou
//This file is writing with Copilot Help (Thanks my guy)
Start: january 22, 11:30

*/
int main(void){

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // To perform matrix-matrix multiplication (C = alpha*A*B + beta*C) we use cblas_dgemm

    //Variables
    double A[] = {
        1.0, 2.0,
        3.0, 4.0
    }; // 2x2
    double B[] = {
        5.0, 6.0,
        7.0, 8.0
    }; // 2x2
    double C[] = {
        0.0, 0.0,
        0.0, 0.0
    }; // 2x2 result

    int M=2, N=2, K=2;
    double alpha = 1.0;
    double beta = 0.0;

    // Expected result: C = A*B = [[19,22],[43,50]]

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, alpha, A, K, B, N, beta, C, N);

    printf("C (dgemm) = [ ");
    for(int i=0;i<M*N;i++){ printf(" %f ",C[i]); }
    printf("]\n");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // To perform symmetric rank-k update (C = alpha*A*A^T + beta*C) we use cblas_dsyrk

    //Variables
    double A2[] = {
        1.0, 2.0,
        3.0, 4.0
    }; // 2x2
    double C2[] = {
        0.0, 0.0,
        0.0, 0.0
    }; // 2x2 result

    int N2=2, K2=2;
    double alpha2=1.0, beta2=0.0;

    // Expected result: C2 = A2*A2^T = [[5,11],[11,25]]

    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans,
                N2, K2, alpha2, A2, K2, beta2, C2, N2);

    printf("C2 (dsyrk, upper) = [ ");
    for(int i=0;i<N2*N2;i++){ printf(" %f ",C2[i]); }
    printf("]\n");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // To solve triangular system (A*X = B) we use cblas_dtrsm

    //Variables
    double A3[] = {
        2.0, 0.0,
        1.0, 3.0
    }; // lower triangular 2x2
    double B3[] = {
        4.0, 10.0,
        6.0, 12.0
    }; // 2x2 right-hand side

    int M3=2, N3=2;
    double alpha3=1.0;

    // Expected result: solve A3*X = B3 → X = A3^{-1}*B3

    cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower,
                CblasNoTrans, CblasNonUnit,
                M3, N3, alpha3, A3, M3, B3, N3);

    printf("X (dtrsm) = [ ");
    for(int i=0;i<M3*N3;i++){ printf(" %f ",B3[i]); }
    printf("]\n");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    return 0;
}

