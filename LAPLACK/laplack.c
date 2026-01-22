#include<stdio.h>
#include <mkl.h>
int main(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// To solve linear system Ax = b we use LAPACKE_dgesv
// Variables
double A[] = {
    2.0, 1.0, 1.0,
    1.0, 3.0, 2.0,
    1.0, 2.0, 2.0
}; // Matrix 3x3 in row-major
double b[] = {4.0, 5.0, 6.0}; // Right hand side vector
int n = 3; // Dimension of matrix A
int nrhs = 1; // Number of right hand sides
int lda = 3; // Leading dimension of A
int ldb = 1; // Leading dimension of b
int ipiv[3]; // Pivot indices
// We must have as result (solving manually):
// 2x + y + z = 4
// x + 3y + 2z = 5
// x + 2y + 2z = 6
// Solution: x = [0.67, -1.0, 3.67]
// Test of LAPACKE_dgesv
lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv, b, ldb);
// To verify our output let's do a print
if (info == 0) {
    printf("Solution: x = [%f, %f, %f]\n", b[0], b[1], b[2]);
    // We've on the output
    // Solution: x = [0.666667, -1.000000, 3.666667] ; It's work!
} else {
    printf("Error: info = %d\n", info);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
     return 0;

}