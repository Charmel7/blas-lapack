# BLAS and LAPACK Tutorial with Intel MKL

A practical guide to using BLAS and LAPACK libraries for linear algebra computations in C, using Intel MKL implementation.

## Overview

This repository contains educational materials and code examples for learning BLAS (Basic Linear Algebra Subprograms) and LAPACK (Linear Algebra PACKage) libraries. These are industry-standard libraries for high-performance linear algebra operations.

## Contents

- `README_BLAS.md` - Complete BLAS documentation covering all three levels
- `README_LAPACK.md` - Complete LAPACK documentation covering main problem categories
- `blas_level1.c` - Vector-vector operations examples
- `blas_level2.c` - Matrix-vector operations examples
- `blas_level3.c` - Matrix-matrix operations examples
- `lapack_exemples.c` - LAPACK routines examples
- `Makefile` - Compilation instructions

## BLAS Hierarchy

**Level 1**: Vector-Vector operations (O(n))
- Addition, scalar product, norms

**Level 2**: Matrix-Vector operations (O(n²))
- Matrix-vector multiplication, rank-1 updates

**Level 3**: Matrix-Matrix operations (O(n³))
- Matrix-matrix multiplication (most optimized)

## LAPACK Categories

**Linear Systems**: Solve Ax = b
**Decompositions**: LU, QR, Cholesky, SVD
**Eigenproblems**: Eigenvalues and eigenvectors
**Least Squares**: Minimize ||Ax - b||

## Prerequisites

- GCC compiler
- Intel MKL library
- Linux environment (Ubuntu recommended)

## Installation

Install Intel MKL:

```bash
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt update
sudo apt install intel-oneapi-mkl-devel
```

Load MKL environment:

```bash
source /opt/intel/oneapi/setvars.sh
```

## Compilation

Compile individual files:

```bash
gcc -o blas_level1 blas_level1.c -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
```

Or use the provided Makefile:

```bash
make
```

## Usage Example

```c
#include <stdio.h>
#include "mkl.h"

int main() {
    // BLAS Level 1: Vector addition y = alpha*x + y
    double x[] = {1.0, 2.0, 3.0};
    double y[] = {4.0, 5.0, 6.0};
    int n = 3;
    double alpha = 2.0;
    
    cblas_daxpy(n, alpha, x, 1, y, 1);
    
    printf("Result: [%f, %f, %f]\n", y[0], y[1], y[2]);
    return 0;
}
```

## Key Functions

### BLAS
- `cblas_daxpy` - Vector addition
- `cblas_ddot` - Dot product
- `cblas_dgemv` - Matrix-vector multiplication
- `cblas_dgemm` - Matrix-matrix multiplication

### LAPACK
- `LAPACKE_dgesv` - Solve linear system
- `LAPACKE_dgetrf` - LU decomposition
- `LAPACKE_dsyev` - Eigenvalues/eigenvectors
- `LAPACKE_dgesvd` - Singular value decomposition

## Documentation

Detailed documentation is provided in:
- `README_BLAS.md` for BLAS functions and conventions
- `README_LAPACK.md` for LAPACK routines and usage

## References

- Intel MKL Documentation: https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c
- BLAS Reference: http://www.netlib.org/blas/
- LAPACK Reference: http://www.netlib.org/lapack/

## License

Educational material for learning purposes.