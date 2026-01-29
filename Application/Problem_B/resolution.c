#include<stdio.h>
#include <mkl.h>
int main(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// To solve least squares we use LAPACKE_dgels
// Variables
double  A[]={2.0,1.0,1.3,
            1.0 ,2.0,2.0,     
            3.0,2.0,1.0  }
; // Matrix 5*4 in row-major
double b[]={100.0,90.0,80.0}; // Right hand side vector
double b1[]={0.0,0.0,0.0}; // Right hand side vector
int n=3;// Dimension of matrix A
int nrhs = 1; // Number of right hand sides
int lda = 3; // Leading dimension of A
int ldb = 1; // Leading dimension of b
int ipiv[3];
// We must have as result (solving manually):

// Test of LAPACKE_dgesv
lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,nrhs,A,lda,ipiv,b1,ldb);
// To verify our output let's do a print
if (info == 0) {
    printf("Solution: x = [%f, %f, %f]\n", b[0], b[1], b[2]);
   
    // We've on the output
   
} else {
    printf("Error: info = %d\n", info);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
     return 0;

}