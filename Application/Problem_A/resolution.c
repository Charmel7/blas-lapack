#include<stdio.h>
#include <mkl.h>
int main(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// To solve least squares we use LAPACKE_dgels
// Variables
double A[] = {
    1.0, 10.0, 3.0,1.2,
    1.0,12.0,4.0,1.1,
    1.0,8.0,6.0,0.9,
    1.0,15.0,5.0,1.3,
    1.0,9.0,6.0,1.0
}; // Matrix 5*4 in row-major
double b[] = {25.0, 28.0,20.0, 32.0, 22.0}; // Right hand side vector
int m = 5;
int n=4;// Dimension of matrix A
int nrhs = 1; // Number of right hand sides
int lda = 4; // Leading dimension of A
int ldb = 1; // Leading dimension of b

// We must have as result (solving manually):

// Test of LAPACKE_dgesv
lapack_int info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,A,lda,b,ldb);
// To verify our output let's do a print
if (info == 0) {
    printf("Solution: x = [%f, %f, %f, %f]\n", b[0], b[1], b[2],b[3]);
   
    // We've on the output
   
} else {
    printf("Error: info = %d\n", info);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
     return 0;

}