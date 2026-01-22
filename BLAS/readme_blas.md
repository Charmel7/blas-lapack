# BLAS - Basic Linear Algebra Subprograms

## Introduction

BLAS is a standardized library for basic linear algebra operations. It provides optimized implementations for computations on vectors and matrices.

### Why Use BLAS?

- Hardware-optimized performance
- Recognized standard in scientific computing
- Used as a foundation by LAPACK and other libraries
- Multiple implementations (Intel MKL, OpenBLAS, ATLAS)

### CBLAS Interface

In C language, the CBLAS interface is used, which is an adaptation of BLAS (originally in FORTRAN).

**Naming Convention:**
```
cblas_dgemm
  |    | |||
  |    | ||+-- Form (m=matrix, v=vector)
  |    | |+--- Operation (ge=general, sy=symmetric, tr=triangular)
  |    | +---- Matrix type
  |    +------ Data type (d=double, s=float, c=complex, z=double complex)
  +----------- CBLAS prefix
```

---

## BLAS Hierarchy

BLAS organizes operations into 3 levels based on their complexity:

| Level | Operation Type | Complexity | Examples |
|-------|----------------|------------|----------|
| Level 1 | Vector-Vector | O(n) | Addition, dot product |
| Level 2 | Matrix-Vector | O(n²) | Ax, triangular solve |
| Level 3 | Matrix-Matrix | O(n³) | Matrix multiplication |

---

## BLAS Level 1: Vector-Vector Operations
BLAS Level 1 refers to all operations that can be performed between two vectors.

### Complexity
O(n) operations for n elements.

### Main Functions

#### cblas_daxpy
**Operation:** y = alpha * x + y
This function performs a linear combination of two vectors x and y.

**Prototype:**
```c
void cblas_daxpy(const int n, const double alpha, const double *x, 
                 const int incx, double *y, const int incy);
```

**Parameters:**
- `n`: number of elements
- `alpha`: scalar multiplier
- `x`: source vector
- `incx`: increment to traverse x (typically 1)
- `y`: destination vector (modified)
- `incy`: increment to traverse y (typically 1)

**Typical Usage:**
```c
// y = 2*x + y
cblas_daxpy(n, 2.0, x, 1, y, 1);
```

---

#### cblas_ddot
**Operation:** result = x · y (dot product)

**Prototype:**
```c
double cblas_ddot(const int n, const double *x, const int incx,
                  const double *y, const int incy);
```

**Parameters:**
- `n`: number of elements
- `x`: first vector
- `incx`: increment for x
- `y`: second vector
- `incy`: increment for y

**Typical Usage:**
```c
// Compute x·y
double product = cblas_ddot(n, x, 1, y, 1);
```

---

#### cblas_dnrm2
**Operation:** result = ||x|| (Euclidean norm)

**Prototype:**
```c
double cblas_dnrm2(const int n, const double *x, const int incx);
```

**Parameters:**
- `n`: number of elements
- `x`: vector
- `incx`: increment for x

**Typical Usage:**
```c
// Compute the norm of x
double norm = cblas_dnrm2(n, x, 1);
```

---

#### cblas_dscal
**Operation:** x = alpha * x

**Prototype:**
```c
void cblas_dscal(const int n, const double alpha, double *x, const int incx);
```

**Parameters:**
- `n`: number of elements
- `alpha`: scalar multiplier
- `x`: vector (modified)
- `incx`: increment for x

**Typical Usage:**
```c
// Multiply all elements of x by 3
cblas_dscal(n, 3.0, x, 1);
```

---

#### cblas_dcopy
**Operation:** y = x

**Prototype:**
```c
void cblas_dcopy(const int n, const double *x, const int incx,
                 double *y, const int incy);
```

**Typical Usage:**
```c
// Copy x into y
cblas_dcopy(n, x, 1, y, 1);
```

---

#### cblas_idamax
**Operation:** finds the index of the largest element in absolute value

**Prototype:**
```c
CBLAS_INDEX cblas_idamax(const int n, const double *x, const int incx);
```

**Typical Usage:**
```c
// Find the index of the largest element
int index = cblas_idamax(n, x, 1);
```

---

## BLAS Level 2: Matrix-Vector Operations

### Complexity
O(n²) operations for n×n matrices.

### Main Functions

#### cblas_dgemv
**Operation:** y = alpha * A * x + beta * y

**Prototype:**
```c
void cblas_dgemv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha,
                 const double *A, const int lda, const double *x,
                 const int incx, const double beta, double *y, const int incy);
```

**Parameters:**
- `layout`: `CblasRowMajor` (C) or `CblasColMajor` (FORTRAN)
- `TransA`: `CblasNoTrans`, `CblasTrans`, or `CblasConjTrans`
- `M`: number of rows in A
- `N`: number of columns in A
- `alpha`: scalar for A*x
- `A`: M×N matrix
- `lda`: leading dimension of A (typically N for RowMajor)
- `x`: vector of size N
- `incx`: increment for x
- `beta`: scalar for y
- `y`: vector of size M (modified)
- `incy`: increment for y

**Typical Usage:**
```c
// y = A*x (without beta*y term)
cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, 1.0, A, N, x, 1, 0.0, y, 1);
```

---

#### cblas_dger
**Operation:** A = alpha * x * y^T + A (rank-1 update)

**Prototype:**
```c
void cblas_dger(const CBLAS_LAYOUT layout, const int M, const int N,
                const double alpha, const double *x, const int incx,
                const double *y, const int incy, double *A, const int lda);
```

**Typical Usage:**
```c
// A = x * y^T + A
cblas_dger(CblasRowMajor, M, N, 1.0, x, 1, y, 1, A, N);
```

---

#### cblas_dtrsv
**Operation:** solves T * x = b where T is triangular

**Prototype:**
```c
void cblas_dtrsv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda,
                 double *x, const int incx);
```

**Additional Parameters:**
- `Uplo`: `CblasUpper` (upper triangular) or `CblasLower` (lower triangular)
- `Diag`: `CblasNonUnit` or `CblasUnit` (diagonal = 1)

---

## BLAS Level 3: Matrix-Matrix Operations

### Complexity
O(n³) operations for n×n matrices. These are the most optimized operations.

### Main Functions

#### cblas_dgemm
**Operation:** C = alpha * A * B + beta * C

**This is THE most important BLAS function.**

**Prototype:**
```c
void cblas_dgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
```

**Parameters:**
- `layout`: `CblasRowMajor` or `CblasColMajor`
- `TransA`: transposition of A
- `TransB`: transposition of B
- `M`: number of rows of A (and C)
- `N`: number of columns of B (and C)
- `K`: number of columns of A / rows of B
- `alpha`: scalar for A*B
- `A`: M×K matrix
- `lda`: leading dimension of A
- `B`: K×N matrix
- `ldb`: leading dimension of B
- `beta`: scalar for C
- `C`: M×N matrix (modified)
- `ldc`: leading dimension of C

**Typical Usage:**
```c
// C = A * B (simple multiplication)
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
            M, N, K, 1.0, A, K, B, N, 0.0, C, N);
```

---

#### cblas_dsyrk
**Operation:** C = alpha * A * A^T + beta * C (symmetric product)

**Prototype:**
```c
void cblas_dsyrk(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
```

**Typical Usage:**
```c
// C = A * A^T
cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, 
            N, K, 1.0, A, K, 0.0, C, N);
```

---

#### cblas_dtrsm
**Operation:** solves T * X = alpha * B where T is triangular

**Prototype:**
```c
void cblas_dtrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
```

---

## Important Conventions

### Layout (memory organization)

**CblasRowMajor (C standard):**
```
2×3 matrix: [a b c]
            [d e f]
            
Memory: [a, b, c, d, e, f]
```

**CblasColMajor (FORTRAN):**
```
2×3 matrix: [a b c]
            [d e f]
            
Memory: [a, d, b, e, c, f]
```

### Leading Dimension (lda, ldb, ldc)

The leading dimension indicates the number of elements between two consecutive rows/columns in memory.

For `CblasRowMajor`: typically equal to the number of columns.

**Example:**
```c
double A[3][5];  // 3×5 matrix
// lda = 5 (5 columns)
```

### Parameter incx, incy

The increment allows skipping elements in a vector.

**incx = 1:** uses all consecutive elements
**incx = 2:** uses every other element

```c
// x = [1, 2, 3, 4, 5, 6]
// incx = 1 → uses [1, 2, 3, 4, 5, 6]
// incx = 2 → uses [1, 3, 5]
```

---

## Compilation with Intel MKL

### Environment variables
```bash
source /opt/intel/oneapi/setvars.sh
```

### Compilation
```bash
gcc -o program program.c \
    -I${MKLROOT}/include \
    -L${MKLROOT}/lib/intel64 \
    -lmkl_intel_lp64 \
    -lmkl_sequential \
    -lmkl_core \
    -lpthread -lm -ldl
```

### Sample Makefile
```makefile
CC = gcc
CFLAGS = -O2 -I${MKLROOT}/include
LDFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

all: blas_level1 blas_level2 blas_level3

blas_level1: blas_level1.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

blas_level2: blas_level2.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

blas_level3: blas_level3.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f blas_level1 blas_level2 blas_level3
```

---

## Summary of Essential Functions

### Level 1 (most used)
- `cblas_daxpy`: linear combination of vectors
- `cblas_ddot`: dot product
- `cblas_dnrm2`: Euclidean norm
- `cblas_dscal`: scalar multiplication

### Level 2 (most used)
- `cblas_dgemv`: matrix-vector multiplication
- `cblas_dger`: outer product
- `cblas_dtrsv`: triangular system solve

### Level 3 (most used)
- `cblas_dgemm`: matrix-matrix multiplication (MOST IMPORTANT)
- `cblas_dsyrk`: symmetric product
- `cblas_dtrsm`: multiple triangular system solve

---

## References

- Intel MKL Documentation: https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c
- BLAS Reference: http://www.netlib.org/blas/
- CBLAS Interface: http://www.netlib.org/blas/#_cblas