#include <stdio.h>
#include <mkl.h>
/*

In this file we're a re exploring the Basic Linear Algebra Subprogram Level 1 for operations between two vectors 
We'll doing some of current examples of use  and explain how it works in the readme  .
Following are juste tests based on cblas ,the c interface  of the native Blas in Fortran77.

Author: Prosper Affoukou 

Start:january 22, 07:14


*/  
int main(void){

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // To perform linear combination(y = alpha*x+y ) on vector y we use cblas_daxpy

    //Variables
    double x[]={1.0,2.0,3.0,4.0};
    double y[]={10.0,20.0,30.0,40.0};
    int n=4;//The dimension of our vectors
    double alpha=10.0;

   // We must have as result:
   //y[]={20.0,40.0,60.0,80.0}

    //Test of cblas_daxpy

    cblas_daxpy(n,alpha,x,1,y,1);

    // To verify our output let's do a print 
    printf("y= [%f,%f,%f,%f] ",y[0],y[1],y[2],y[3]);

    // We've on the output 
    //y= [20.000000,40.000000,60.000000,80.000000] ; It's work!

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// To compute dot product (scalar product: result = x·y) we use cblas_ddot

// Variables
double a[] = {1.0, 2.0, 3.0, 4.0};
double b[] = {5.0, 6.0, 7.0, 8.0};
int m = 4; // The dimension of our vectors


// We must have as result:
// result = 1*5 + 2*6 + 3*7 + 4*8 = 5 + 12 + 21 + 32 = 70

// Test of cblas_ddot
double result = cblas_ddot(m, a, 1, b, 1);

// To verify our output let's do a print
printf("dot product = %f\n", result);
// We've on the output
// dot product = 70.000000 ; It's work!

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// To compute euclidean norm (result = ||x|| = sqrt(x[0]² + x[1]² + ... + x[n-1]²)) we use cblas_dnrm2

// Variables
double c[] = {3.0, 4.0};
int o = 2; // The dimension of our vector


// We must have as result:
// result = sqrt(3² + 4²) = sqrt(9 + 16) = sqrt(25) = 5.0

// Test of cblas_dnrm2
double norme = cblas_dnrm2(o, c, 1);


// To verify our output let's do a print
printf("norm = %f\n", norme);
// We've on the output
// norm = 5.000000 ; It's work!
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// To perform scalar multiplication (d = alpha*d) on vector d we use cblas_dscal

// Variables
double d[] = {1.0, 2.0, 3.0, 4.0};
int p = 4; // The dimension of our vector
double alpha2 = 3.0;

// We must have as result:
// d[] = {3.0, 6.0, 9.0, 12.0}
// Test of cblas_dscal
cblas_dscal(p, alpha2, d, 1);

// To verify our output let's do a print
printf("x = [%f, %f, %f, %f]\n", d[0], d[1], d[2], d[3]);
// We've on the output
// d = [3.000000, 6.000000, 9.000000, 12.000000] ; It's work!
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// To copy vector e into f (e = f) we use cblas_dcopy

// Variables
double e[] = {1.0, 2.0, 3.0, 4.0};
double f[] = {0.0, 0.0, 0.0, 0.0};
int q = 4; // The dimension of our vectors

// We must have as result:
// f[] = {1.0, 2.0, 3.0, 4.0}
// Test of cblas_dcopy
cblas_dcopy(q, e, 1, f, 1);

// To verify our output let's do a print
printf("y = [%f, %f, %f, %f]\n", f[0], f[1], f[2], f[3]);
// We've on the output
// f = [1.000000, 2.000000, 3.000000, 4.000000] ; It's work!
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// To find index of maximum absolute value we use cblas_idamax

// Variables
double g[] = {-5.0, 2.0, -8.0, 3.0};
int r = 4; // The dimension of our vector

// We must have as result:
// index = 2 (because |-8.0| = 8.0 is the largest absolute value)
// Test of cblas_idamax
int index = cblas_idamax(r, g, 1);

// To verify our output let's do a print
printf("index of max absolute value = %d, value = %f\n", index, g[index]);
// We've on the output
// index of max absolute value = 2, value = -8.000000 ; It's work!
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
printf("What's a wonderful work !");

printf("Yeah I know I'm the Best ");
    return 0;
}