# BLAS - Basic Linear Algebra Subprograms

## Introduction

BLAS est une bibliothèque standardisée pour les opérations d'algèbre linéaire de base. Elle fournit des implémentations optimisées pour les calculs sur vecteurs et matrices.

### Pourquoi utiliser BLAS ?

- Performance optimisée au niveau matériel
- Standard reconnu dans le calcul scientifique
- Utilisé comme base par LAPACK et autres bibliothèques
- Implémentations multiples (Intel MKL, OpenBLAS, ATLAS)

### Interface CBLAS

En langage C, on utilise l'interface CBLAS qui est une adaptation de BLAS (originellement en FORTRAN).

**Convention de nommage :**
```
cblas_dgemm
  |    | |||
  |    | ||+-- Forme (m=matrix, v=vector)
  |    | |+--- Opération (ge=general, sy=symmetric, tr=triangular)
  |    | +---- Type matrice
  |    +------ Type de données (d=double, s=float, c=complex, z=double complex)
  +----------- Préfixe CBLAS
```

---

## Hiérarchie BLAS

BLAS organise les opérations en 3 niveaux selon leur complexité :

| Niveau | Type d'opération | Complexité | Exemples |
|--------|------------------|------------|----------|
| Level 1 | Vecteur-Vecteur | O(n) | Addition, produit scalaire |
| Level 2 | Matrice-Vecteur | O(n²) | Ax, résolution triangulaire |
| Level 3 | Matrice-Matrice | O(n³) | Multiplication matricielle |

---

## BLAS Level 1 : Opérations Vecteur-Vecteur
 BLAS Level 1 dsigne toutes les operations que nous pouvons realiser entre deux vecteurs .

### Complexité
O(n) opérations pour n éléments.

### Principales fonctions

#### cblas_daxpy
**Opération :** y = alpha * x + y
Cette fonction realise une combinaison lineaire de deux vecteur x et y.

**Prototype :**
```c
void cblas_daxpy(const int n, const double alpha, const double *x, 
                 const int incx, double *y, const int incy);
```

**Paramètres :**
- `n` : nombre d'éléments
- `alpha` : scalaire multiplicateur
- `x` : vecteur source
- `incx` : incrément pour parcourir x (généralement 1)
- `y` : vecteur destination (modifié)
- `incy` : incrément pour parcourir y (généralement 1)

**Usage typique :**
```c
// y = 2*x + y
cblas_daxpy(n, 2.0, x, 1, y, 1);
```

---

#### cblas_ddot
**Opération :** résultat = x · y (produit scalaire)

**Prototype :**
```c
double cblas_ddot(const int n, const double *x, const int incx,
                  const double *y, const int incy);
```

**Paramètres :**
- `n` : nombre d'éléments
- `x` : premier vecteur
- `incx` : incrément pour x
- `y` : second vecteur
- `incy` : incrément pour y

**Usage typique :**
```c
// Calculer x·y
double produit = cblas_ddot(n, x, 1, y, 1);
```

---

#### cblas_dnrm2
**Opération :** résultat = ||x|| (norme euclidienne)

**Prototype :**
```c
double cblas_dnrm2(const int n, const double *x, const int incx);
```

**Paramètres :**
- `n` : nombre d'éléments
- `x` : vecteur
- `incx` : incrément pour x

**Usage typique :**
```c
// Calculer la norme de x
double norme = cblas_dnrm2(n, x, 1);
```

---

#### cblas_dscal
**Opération :** x = alpha * x

**Prototype :**
```c
void cblas_dscal(const int n, const double alpha, double *x, const int incx);
```

**Paramètres :**
- `n` : nombre d'éléments
- `alpha` : scalaire multiplicateur
- `x` : vecteur (modifié)
- `incx` : incrément pour x

**Usage typique :**
```c
// Multiplier tous les éléments de x par 3
cblas_dscal(n, 3.0, x, 1);
```

---

#### cblas_dcopy
**Opération :** y = x

**Prototype :**
```c
void cblas_dcopy(const int n, const double *x, const int incx,
                 double *y, const int incy);
```

**Usage typique :**
```c
// Copier x dans y
cblas_dcopy(n, x, 1, y, 1);
```

---

#### cblas_idamax
**Opération :** trouve l'indice du plus grand élément en valeur absolue

**Prototype :**
```c
CBLAS_INDEX cblas_idamax(const int n, const double *x, const int incx);
```

**Usage typique :**
```c
// Trouver l'indice du plus grand élément
int index = cblas_idamax(n, x, 1);
```

---

## BLAS Level 2 : Opérations Matrice-Vecteur

### Complexité
O(n²) opérations pour matrices n×n.

### Principales fonctions

#### cblas_dgemv
**Opération :** y = alpha * A * x + beta * y

**Prototype :**
```c
void cblas_dgemv(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA,
                 const int M, const int N, const double alpha,
                 const double *A, const int lda, const double *x,
                 const int incx, const double beta, double *y, const int incy);
```

**Paramètres :**
- `layout` : `CblasRowMajor` (C) ou `CblasColMajor` (FORTRAN)
- `TransA` : `CblasNoTrans`, `CblasTrans`, ou `CblasConjTrans`
- `M` : nombre de lignes de A
- `N` : nombre de colonnes de A
- `alpha` : scalaire pour A*x
- `A` : matrice M×N
- `lda` : leading dimension de A (généralement N pour RowMajor)
- `x` : vecteur de taille N
- `incx` : incrément pour x
- `beta` : scalaire pour y
- `y` : vecteur de taille M (modifié)
- `incy` : incrément pour y

**Usage typique :**
```c
// y = A*x (sans terme beta*y)
cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, 1.0, A, N, x, 1, 0.0, y, 1);
```

---

#### cblas_dger
**Opération :** A = alpha * x * y^T + A (mise à jour de rang 1)

**Prototype :**
```c
void cblas_dger(const CBLAS_LAYOUT layout, const int M, const int N,
                const double alpha, const double *x, const int incx,
                const double *y, const int incy, double *A, const int lda);
```

**Usage typique :**
```c
// A = x * y^T + A
cblas_dger(CblasRowMajor, M, N, 1.0, x, 1, y, 1, A, N);
```

---

#### cblas_dtrsv
**Opération :** résout T * x = b où T est triangulaire

**Prototype :**
```c
void cblas_dtrsv(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
                 const int N, const double *A, const int lda,
                 double *x, const int incx);
```

**Paramètres supplémentaires :**
- `Uplo` : `CblasUpper` (triangulaire supérieure) ou `CblasLower` (inférieure)
- `Diag` : `CblasNonUnit` ou `CblasUnit` (diagonale = 1)

---

## BLAS Level 3 : Opérations Matrice-Matrice

### Complexité
O(n³) opérations pour matrices n×n. Ce sont les opérations les plus optimisées.

### Principales fonctions

#### cblas_dgemm
**Opération :** C = alpha * A * B + beta * C

**C'est LA fonction la plus importante de BLAS.**

**Prototype :**
```c
void cblas_dgemm(const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc);
```

**Paramètres :**
- `layout` : `CblasRowMajor` ou `CblasColMajor`
- `TransA` : transposition de A
- `TransB` : transposition de B
- `M` : nombre de lignes de A (et C)
- `N` : nombre de colonnes de B (et C)
- `K` : nombre de colonnes de A / lignes de B
- `alpha` : scalaire pour A*B
- `A` : matrice M×K
- `lda` : leading dimension de A
- `B` : matrice K×N
- `ldb` : leading dimension de B
- `beta` : scalaire pour C
- `C` : matrice M×N (modifié)
- `ldc` : leading dimension de C

**Usage typique :**
```c
// C = A * B (multiplication simple)
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
            M, N, K, 1.0, A, K, B, N, 0.0, C, N);
```

---

#### cblas_dsyrk
**Opération :** C = alpha * A * A^T + beta * C (produit symétrique)

**Prototype :**
```c
void cblas_dsyrk(const CBLAS_LAYOUT layout, const CBLAS_UPLO Uplo,
                 const CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);
```

**Usage typique :**
```c
// C = A * A^T
cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, 
            N, K, 1.0, A, K, 0.0, C, N);
```

---

#### cblas_dtrsm
**Opération :** résout T * X = alpha * B où T est triangulaire

**Prototype :**
```c
void cblas_dtrsm(const CBLAS_LAYOUT layout, const CBLAS_SIDE Side,
                 const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_DIAG Diag, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 double *B, const int ldb);
```

---

## Conventions importantes

### Layout (organisation mémoire)

**CblasRowMajor (C standard) :**
```
Matrice 2×3 : [a b c]
              [d e f]
              
Mémoire : [a, b, c, d, e, f]
```

**CblasColMajor (FORTRAN) :**
```
Matrice 2×3 : [a b c]
              [d e f]
              
Mémoire : [a, d, b, e, c, f]
```

### Leading Dimension (lda, ldb, ldc)

Le leading dimension indique le nombre d'éléments entre deux lignes/colonnes consécutives en mémoire.

Pour `CblasRowMajor` : généralement égal au nombre de colonnes.

**Exemple :**
```c
double A[3][5];  // Matrice 3×5
// lda = 5 (5 colonnes)
```

### Paramètre incx, incy

L'incrément permet de sauter des éléments dans un vecteur.

**incx = 1 :** utilise tous les éléments consécutifs
**incx = 2 :** utilise un élément sur deux

```c
// x = [1, 2, 3, 4, 5, 6]
// incx = 1 → utilise [1, 2, 3, 4, 5, 6]
// incx = 2 → utilise [1, 3, 5]
```

---

## Compilation avec Intel MKL

### Variables d'environnement
```bash
source /opt/intel/oneapi/setvars.sh
```

### Compilation
```bash
gcc -o programme programme.c \
    -I${MKLROOT}/include \
    -L${MKLROOT}/lib/intel64 \
    -lmkl_intel_lp64 \
    -lmkl_sequential \
    -lmkl_core \
    -lpthread -lm -ldl
```

### Makefile type
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

## Résumé des fonctions essentielles

### Level 1 (les plus utilisées)
- `cblas_daxpy` : combinaison linéaire de vecteurs
- `cblas_ddot` : produit scalaire
- `cblas_dnrm2` : norme euclidienne
- `cblas_dscal` : multiplication par scalaire

### Level 2 (les plus utilisées)
- `cblas_dgemv` : multiplication matrice-vecteur
- `cblas_dger` : produit extérieur
- `cblas_dtrsv` : résolution système triangulaire

### Level 3 (les plus utilisées)
- `cblas_dgemm` : multiplication matrice-matrice (LA PLUS IMPORTANTE)
- `cblas_dsyrk` : produit symétrique
- `cblas_dtrsm` : résolution système triangulaire multiple

---

## Références

- Documentation Intel MKL : https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c
- BLAS Reference : http://www.netlib.org/blas/
- CBLAS Interface : http://www.netlib.org/blas/#_cblas