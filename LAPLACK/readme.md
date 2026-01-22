# LAPACK - Linear Algebra PACKage

## Introduction

LAPACK is a library of routines for solving common numerical linear algebra problems. It relies on BLAS for basic operations.

### BLAS-LAPACK Relationship

```
LAPACK (advanced algorithms)
    ↓ uses
BLAS (basic operations)
    ↓ uses
Hardware (CPU, cache, SIMD)
```

LAPACK doesn't reinvent matrix operations; it organizes them intelligently using BLAS.

---

## Areas Covered by LAPACK

| Area | Problems Solved | Complexity |
|------|-----------------|------------|
| Linear systems | Ax = b | O(n³) |
| Decompositions | LU, QR, Cholesky, SVD | O(n³) |
| Eigenvalues | Eigenvalues, eigenvectors | O(n³) |
| Least squares | min ‖Ax - b‖ | O(mn²) |

---

## LAPACK Naming Convention

LAPACK functions follow a strict pattern:

```
dgesv
│││└─ Problem type
││└── Matrix type
│└─── Data type
└──── (nothing = driver routine)
```

### Data Type Prefix

- **d** = double precision
- **s** = single precision
- **c** = complex single
- **z** = complex double

### Matrix Type

- **ge** = GEneral (any matrix)
- **sy** = SYmmetric (symmetric)
- **po** = POsitive definite (positive definite)
- **tr** = TRiangular (triangular)
- **he** = HErmitian (Hermitian, for complex)
- **gb** = General Band (general band)

### Problem Type

- **sv** = SolVe (solve linear system)
- **trf** = TRiangular Factorization (decomposition)
- **trs** = TRiangular Solve (solve with already factorized matrix)
- **ev** = EigenValues (eigenvalues)
- **ls** = Least Squares (least squares)
- **qrf** = QR Factorization
- **potrf** = POsitive definite TRiangular Factorization (Cholesky)

---

## Category 1: Linear System Solving

### Problem: Solve Ax = b

#### LAPACK_dgesv (Driver routine - All-in-one)

**Operation:** Solves Ax = b for general matrix

**Prototype:**
```c
lapack_int LAPACKE_dgesv(int matrix_layout, lapack_int n, lapack_int nrhs,
                         double* A, lapack_int lda, lapack_int* ipiv,
                         double* b, lapack_int ldb);
```

**Parameters:**
- `matrix_layout`: LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
- `n`: dimension of matrix A (n×n)
- `nrhs`: number of right-hand sides (columns of b)
- `A`: system matrix (will be modified - contains LU after)
- `lda`: leading dimension of A
- `ipiv`: pivot vector (size n)
- `b`: right-hand side input, solution output
- `ldb`: leading dimension of b

**Return:** 0 if success, >0 if matrix singular

**What dgesv does:**
1. Decomposes A into LU with pivoting
2. Solves lower triangular system
3. Solves upper triangular system
4. Returns solution in b

---

#### LAPACKE_dgetrf + LAPACKE_dgetrs (Two-step approach)

**Advantage:** Reuse decomposition for multiple right-hand sides

##### LAPACKE_dgetrf - LU Decomposition

**Operation:** A = P*L*U

**Prototype:**
```c
lapack_int LAPACKE_dgetrf(int matrix_layout, lapack_int m, lapack_int n,
                          double* A, lapack_int lda, lapack_int* ipiv);
```

**Parameters:**
- `m`: number of rows of A
- `n`: number of columns of A
- `A`: input matrix, contains L and U in output
- `ipiv`: pivot vector

##### LAPACKE_dgetrs - Solve with LU

**Operation:** Solves A*x = b using LU decomposition

**Prototype:**
```c
lapack_int LAPACKE_dgetrs(int matrix_layout, char trans, lapack_int n,
                          lapack_int nrhs, const double* A, lapack_int lda,
                          const lapack_int* ipiv, double* b, lapack_int ldb);
```

**Parameters:**
- `trans`: 'N' (no transposition), 'T' (transpose), 'C' (conjugate)
- Other parameters identical to dgesv

---

#### LAPACKE_dposv (Symmetric positive definite matrices)

**Operation:** Solves Ax = b with symmetric positive definite A (uses Cholesky)

**Prototype:**
```c
lapack_int LAPACKE_dposv(int matrix_layout, char uplo, lapack_int n,
                         lapack_int nrhs, double* A, lapack_int lda,
                         double* b, lapack_int ldb);
```

**Parameters:**
- `uplo`: 'U' (upper) or 'L' (lower) - which half of A is stored

**Advantage over dgesv:**
- Faster (~2x factor)
- More numerically stable
- But requires A to be positive definite

---

## Category 2: Matrix Decompositions

### LU Decomposition

Already seen with `LAPACKE_dgetrf`

**A = P*L*U**
- P = permutation matrix
- L = lower triangular (diagonal = 1)
- U = upper triangular

---

### QR Decomposition

#### LAPACKE_dgeqrf

**Operation:** A = Q*R

**Prototype:**
```c
lapack_int LAPACKE_dgeqrf(int matrix_layout, lapack_int m, lapack_int n,
                          double* A, lapack_int lda, double* tau);
```

**Parameters:**
- `m`: number of rows
- `n`: number of columns
- `A`: input matrix, contains R and Q factors in output
- `tau`: vector of scalars to reconstruct Q (size min(m,n))

**Utility:**
- Solving overdetermined systems (m > n)
- Computing orthonormal bases
- Least squares

---

### Cholesky Decomposition

#### LAPACKE_dpotrf

**Operation:** A = L*L^T (or A = U^T*U)

**Prototype:**
```c
lapack_int LAPACKE_dpotrf(int matrix_layout, char uplo, lapack_int n,
                          double* A, lapack_int lda);
```

**Parameters:**
- `uplo`: 'U' or 'L'
- `n`: dimension of A
- `A`: symmetric positive definite matrix input, factor L or U output

**Condition:** A must be symmetric positive definite

**Utility:**
- Fast solution of systems with positive definite matrices
- Simulation of correlated Gaussian random variables
- Optimization (Hessian matrices)

---

### Singular Value Decomposition (SVD)

#### LAPACKE_dgesvd

**Operation:** A = U*Σ*V^T

**Prototype:**
```c
lapack_int LAPACKE_dgesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, double* A,
                          lapack_int lda, double* s, double* u,
                          lapack_int ldu, double* vt, lapack_int ldvt,
                          double* superb);
```

**Parameters:**
- `jobu`: 'A' (all columns of U), 'S' (min(m,n) columns), 'N' (none)
- `jobvt`: same for V^T
- `s`: vector of singular values (size min(m,n))
- `u`: matrix U (m×m or m×min(m,n))
- `vt`: matrix V^T (n×n or min(m,n)×n)
- `superb`: workspace

**Utility:**
- Principal Component Analysis (PCA)
- Data compression
- Moore-Penrose pseudo-inverse
- Solving ill-conditioned systems
- Low-rank approximation

---

## Category 3: Eigenvalue Problems

### LAPACKE_dsyev (Symmetric matrices)

**Operation:** Find λ and v such that A*v = λ*v

**Prototype:**
```c
lapack_int LAPACKE_dsyev(int matrix_layout, char jobz, char uplo,
                         lapack_int n, double* A, lapack_int lda,
                         double* w);
```

**Parameters:**
- `jobz`: 'N' (eigenvalues only), 'V' (eigenvalues and eigenvectors)
- `uplo`: 'U' or 'L'
- `A`: symmetric matrix input, eigenvectors output if jobz='V'
- `w`: vector of eigenvalues (size n)

**Guarantee:** For symmetric matrices, eigenvalues are always real

---

### LAPACKE_dgeev (General matrices)

**Operation:** Find λ and v such that A*v = λ*v

**Prototype:**
```c
lapack_int LAPACKE_dgeev(int matrix_layout, char jobvl, char jobvr,
                         lapack_int n, double* A, lapack_int lda,
                         double* wr, double* wi, double* vl,
                         lapack_int ldvl, double* vr, lapack_int ldvr);
```

**Parameters:**
- `jobvl`: compute left eigenvectors
- `jobvr`: compute right eigenvectors
- `wr`: real part of eigenvalues
- `wi`: imaginary part of eigenvalues
- `vl`: left eigenvectors
- `vr`: right eigenvectors

**Warning:** Eigenvalues can be complex even for real matrices

---

## Category 4: Least Squares

### LAPACKE_dgels

**Operation:** Solves min ‖A*x - b‖₂

**Prototype:**
```c
lapack_int LAPACKE_dgels(int matrix_layout, char trans, lapack_int m,
                         lapack_int n, lapack_int nrhs, double* A,
                         lapack_int lda, double* b, lapack_int ldb);
```

**Parameters:**
- `trans`: 'N' (solves min ‖Ax-b‖), 'T' (solves min ‖A^T*x-b‖)
- `m`: number of rows of A
- `n`: number of columns of A
- `b`: vector b input, solution x output (size max(m,n))

**Use cases:**
- Linear regression
- Curve fitting
- Overdetermined systems (more equations than unknowns)

---

## Error Handling

All LAPACKE functions return an integer:

**Return value:**
- **0**: success
- **< 0**: the i-th argument has an illegal value (|return| = i)
- **> 0**: algorithm failure (meaning depends on function)

**Example:**
```c
lapack_int info;
info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv, b, ldb);

if (info == 0) {
    printf("Success\n");
} else if (info < 0) {
    printf("Error: invalid argument %d\n", -info);
} else {
    printf("Singular matrix at element %d\n", info);
}
```

---

## Compilation with Intel MKL

### Required Headers
```c
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
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

---

## Best Practices

### 1. Copy Matrices if Needed

Most LAPACK functions **modify** the input matrix:

```c
// A will be destroyed by dgesv
double A_original[n*n];
double A_copy[n*n];
memcpy(A_copy, A_original, n*n*sizeof(double));
LAPACKE_dgesv(..., A_copy, ...);
// A_original preserved, A_copy contains LU
```

### 2. Check Return Values

```c
lapack_int info = LAPACKE_dgesv(...);
if (info != 0) {
    // Handle error
}
```

### 3. Use the Right Routines

- Symmetric positive definite matrix → `dposv` (Cholesky)
- General matrix → `dgesv` (LU)
- Multiple right-hand sides with same A → `dgetrf` + `dgetrs`

### 4. Prefer Driver Routines

"Driver routines" (like `dgesv`, `dposv`) do all the work in one call. Use separate routines only if necessary.

---

## Summary Table

| Problem | Routine | Matrix Type | Complexity |
|---------|---------|-------------|------------|
| Linear system | `dgesv` | General | O(n³) |
| Linear system | `dposv` | Symmetric positive definite | O(n³/3) |
| LU decomposition | `dgetrf` | General | O(n³) |
| QR decomposition | `dgeqrf` | Any | O(mn²) |
| Cholesky decomposition | `dpotrf` | Sym. positive definite | O(n³/3) |
| SVD | `dgesvd` | Any | O(mn²) |
| Eigenvalues | `dsyev` | Symmetric | O(n³) |
| Eigenvalues | `dgeev` | General | O(n³) |
| Least squares | `dgels` | Any | O(mn²) |

---

## References

- Intel MKL LAPACK Documentation: https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/lapack-routines.html
- LAPACK Users' Guide: http://www.netlib.org/lapack/lug/
- LAPACK Reference: http://www.netlib.org/lapack/