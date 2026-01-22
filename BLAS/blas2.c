#include <stdio.h>
#include <mkl.h>
/*

In this file we're a re exploring the Basic Linear Algebra Subprogram Level 2 for operations between vector and matrix 
We'll doing some of current examples of use  and explain how it works in the readme  .
Following are juste tests based on cblas ,the c interface  of the native Blas in Fortran77.

Author: Prosper Affoukou 

Start:january 22, 09:09


*/  
int main(void){

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // To perform matrix-vector multiplication (y = alpha*A*x + beta*y) we use cblas_dgemv

    //Variables
    double A[] = {
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0
};
    double x[] = {1.0, 2.0}; 
    double y[] = {0.0, 0.0,0.0};
    int M=2,N=3;
    int n=4;
   

    double alpha = 1.0;
    double beta = 0.0;
    // We must have as result:
    //y[]={14.00,32.00}

    //Test of cblas_dgemv 

    cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, alpha, A, N, x, 1, beta, y, 1);

    // To verify our output let's do a print 
    printf("y= [%f,%f] ",y[0],y[1]);

    // We've on the output

    //y= [14.000000,32.000000]  ; It's work!

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // To perform matrix-vector outer_product (A= alpha*x*y^T+ A) we use cblas_dger

    //Variables
    double A1[] = {
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0
};
    double x1[] = {1.0, 2.0}; 
    double y1[] = {0.0, 0.0,0.0};
    int M1=2,N1=3;
    int n1=4;
   

    double alpha1 = 2.0;
    // The result must be A again coz the first term must be zero matrix
    // Test of cbals_dger (A =alpha* x * y^T + A)
    cblas_dger(CblasRowMajor, M1, N1,alpha1, x1, 1, y1, 1, A1, N1);
    printf("\n A=[ ");
    for (int i=0 ; i< (M1*N1);i++){
            printf(" %f ",A[i]);
    }
    printf("]\n");

    // We've on the output

    //A=[  1.000000  2.000000  3.000000  4.000000  5.000000  6.000000 ]  ; It's work!


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // To solve triangular T*x=b

    //Variables
    double T[] = {
    1.0, 2.0, 3.0,
    0.0, 5.0, 6.0,
    0.0, 0.0, 3.0
};
    double b[]={1.0,1.0,1.0};

    //The result must be b=[  0.400000  -0.200000  0.333333 ]
    
    //Test of cblas_dtrsv
    double result[]={0.0,0.0,0.0};
    cblas_dtrsv(CblasRowMajor,CblasUpper,CblasNoTrans,CblasNonUnit,3,T,3,b,1);
    
    printf("\n b=[ ");
    for (int i=0 ; i<3 ; i++){
            printf(" %f ",b[i]);
    }
    printf("]\n");

    // We've on the output

    //b=[  0.400000  -0.200000  0.333333 ]  ; It's work!


   return 0;

    printf("By The Bee Guy ! ");
}