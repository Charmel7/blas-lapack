# LAPACK - Linear Algebra PACKage

## Introduction

LAPACK est une bibliothèque de routines pour résoudre des problèmes courants d'algèbre linéaire numérique. Elle s'appuie sur BLAS pour les opérations de base.

### Relation BLAS-LAPACK

```
LAPACK (algorithmes avancés)
    ↓ utilise
BLAS (opérations de base)
    ↓ utilise
Matériel (CPU, cache, SIMD)
```

LAPACK ne réinvente pas les opérations matricielles, elle les organise intelligemment en utilisant BLAS.

---

## Domaines couverts par LAPACK

| Domaine | Problèmes résolus | Complexité |
|---------|-------------------|------------|
| Systèmes linéaires | Ax = b | O(n³) |
| Décompositions | LU, QR, Cholesky, SVD | O(n³) |
| Valeurs propres | Eigenvalues, eigenvectors | O(n³) |
| Moindres carrés | min ‖Ax - b‖ | O(mn²) |

---

## Convention de nommage LAPACK

Les fonctions LAPACK suivent un pattern strict:

```
dgesv
│││└─ Type de problème
││└── Type de matrice
│└─── Type de données
└──── (rien = driver routine)
```

### Préfixe de type de données

- **d** = double précision
- **s** = simple précision
- **c** = complexe simple
- **z** = complexe double

### Type de matrice

- **ge** = GEneral (matrice quelconque)
- **sy** = SYmmetric (symétrique)
- **po** = POsitive definite (définie positive)
- **tr** = TRiangular (triangulaire)
- **he** = HErmitian (hermitienne, pour complexes)
- **gb** = General Band (bande générale)

### Type de problème

- **sv** = SolVe (résolution système linéaire)
- **trf** = TRiangular Factorization (décomposition)
- **trs** = TRiangular Solve (résolution avec matrice déjà factorisée)
- **ev** = EigenValues (valeurs propres)
- **ls** = Least Squares (moindres carrés)
- **qrf** = QR Factorization
- **potrf** = POsitive definite TRiangular Factorization (Cholesky)

---

## Catégorie 1: Résolution de systèmes linéaires

### Problème: Résoudre Ax = b

#### LAPACK_dgesv (Driver routine - Tout en un)

**Opération:** Résout Ax = b pour matrice générale

**Prototype:**
```c
lapack_int LAPACKE_dgesv(int matrix_layout, lapack_int n, lapack_int nrhs,
                         double* A, lapack_int lda, lapack_int* ipiv,
                         double* b, lapack_int ldb);
```

**Paramètres:**
- `matrix_layout` : LAPACK_ROW_MAJOR ou LAPACK_COL_MAJOR
- `n` : dimension de la matrice A (n×n)
- `nrhs` : nombre de seconds membres (colonnes de b)
- `A` : matrice du système (sera modifiée - contient LU après)
- `lda` : leading dimension de A
- `ipiv` : vecteur de pivots (taille n)
- `b` : second membre en entrée, solution en sortie
- `ldb` : leading dimension de b

**Retour:** 0 si succès, >0 si matrice singulière

**Ce que fait dgesv:**
1. Décompose A en LU avec pivotage
2. Résout le système triangulaire inférieur
3. Résout le système triangulaire supérieur
4. Retourne la solution dans b

---

#### LAPACKE_dgetrf + LAPACKE_dgetrs (Approche en deux étapes)

**Avantage:** Réutiliser la décomposition pour plusieurs seconds membres

##### LAPACKE_dgetrf - Décomposition LU

**Opération:** A = P*L*U

**Prototype:**
```c
lapack_int LAPACKE_dgetrf(int matrix_layout, lapack_int m, lapack_int n,
                          double* A, lapack_int lda, lapack_int* ipiv);
```

**Paramètres:**
- `m` : nombre de lignes de A
- `n` : nombre de colonnes de A
- `A` : matrice en entrée, contient L et U en sortie
- `ipiv` : vecteur de pivots

##### LAPACKE_dgetrs - Résolution avec LU

**Opération:** Résout A*x = b en utilisant la décomposition LU

**Prototype:**
```c
lapack_int LAPACKE_dgetrs(int matrix_layout, char trans, lapack_int n,
                          lapack_int nrhs, const double* A, lapack_int lda,
                          const lapack_int* ipiv, double* b, lapack_int ldb);
```

**Paramètres:**
- `trans` : 'N' (pas de transposition), 'T' (transposée), 'C' (conjuguée)
- Les autres paramètres sont identiques à dgesv

---

#### LAPACKE_dposv (Matrices symétriques définies positives)

**Opération:** Résout Ax = b avec A symétrique définie positive (utilise Cholesky)

**Prototype:**
```c
lapack_int LAPACKE_dposv(int matrix_layout, char uplo, lapack_int n,
                         lapack_int nrhs, double* A, lapack_int lda,
                         double* b, lapack_int ldb);
```

**Paramètres:**
- `uplo` : 'U' (upper) ou 'L' (lower) - quelle moitié de A est stockée

**Avantage sur dgesv:**
- Plus rapide (facteur ~2)
- Plus stable numériquement
- Mais nécessite que A soit définie positive

---

## Catégorie 2: Décompositions matricielles

### Décomposition LU

Déjà vu avec `LAPACKE_dgetrf`

**A = P*L*U**
- P = matrice de permutation
- L = triangulaire inférieure (diagonal = 1)
- U = triangulaire supérieure

---

### Décomposition QR

#### LAPACKE_dgeqrf

**Opération:** A = Q*R

**Prototype:**
```c
lapack_int LAPACKE_dgeqrf(int matrix_layout, lapack_int m, lapack_int n,
                          double* A, lapack_int lda, double* tau);
```

**Paramètres:**
- `m` : nombre de lignes
- `n` : nombre de colonnes
- `A` : matrice en entrée, contient R et facteurs de Q en sortie
- `tau` : vecteur de scalaires pour reconstruire Q (taille min(m,n))

**Utilité:**
- Résolution de systèmes surdéterminés (m > n)
- Calcul de bases orthonormées
- Moindres carrés

---

### Décomposition de Cholesky

#### LAPACKE_dpotrf

**Opération:** A = L*L^T (ou A = U^T*U)

**Prototype:**
```c
lapack_int LAPACKE_dpotrf(int matrix_layout, char uplo, lapack_int n,
                          double* A, lapack_int lda);
```

**Paramètres:**
- `uplo` : 'U' ou 'L'
- `n` : dimension de A
- `A` : matrice symétrique définie positive en entrée, facteur L ou U en sortie

**Condition:** A doit être symétrique définie positive

**Utilité:**
- Résolution rapide de systèmes avec matrices définies positives
- Simulation de variables aléatoires gaussiennes corrélées
- Optimisation (matrices hessiennes)

---

### Décomposition en valeurs singulières (SVD)

#### LAPACKE_dgesvd

**Opération:** A = U*Σ*V^T

**Prototype:**
```c
lapack_int LAPACKE_dgesvd(int matrix_layout, char jobu, char jobvt,
                          lapack_int m, lapack_int n, double* A,
                          lapack_int lda, double* s, double* u,
                          lapack_int ldu, double* vt, lapack_int ldvt,
                          double* superb);
```

**Paramètres:**
- `jobu` : 'A' (toutes les colonnes de U), 'S' (min(m,n) colonnes), 'N' (aucune)
- `jobvt` : idem pour V^T
- `s` : vecteur des valeurs singulières (taille min(m,n))
- `u` : matrice U (m×m ou m×min(m,n))
- `vt` : matrice V^T (n×n ou min(m,n)×n)
- `superb` : workspace

**Utilité:**
- Analyse en composantes principales (PCA)
- Compression de données
- Pseudo-inverse de Moore-Penrose
- Résolution de systèmes mal conditionnés
- Approximation de rang faible

---

## Catégorie 3: Problèmes de valeurs propres

### LAPACKE_dsyev (Matrices symétriques)

**Opération:** Trouve λ et v tels que A*v = λ*v

**Prototype:**
```c
lapack_int LAPACKE_dsyev(int matrix_layout, char jobz, char uplo,
                         lapack_int n, double* A, lapack_int lda,
                         double* w);
```

**Paramètres:**
- `jobz` : 'N' (valeurs propres seulement), 'V' (valeurs et vecteurs)
- `uplo` : 'U' ou 'L'
- `A` : matrice symétrique en entrée, vecteurs propres en sortie si jobz='V'
- `w` : vecteur des valeurs propres (taille n)

**Garantie:** Pour matrices symétriques, valeurs propres toujours réelles

---

### LAPACKE_dgeev (Matrices générales)

**Opération:** Trouve λ et v tels que A*v = λ*v

**Prototype:**
```c
lapack_int LAPACKE_dgeev(int matrix_layout, char jobvl, char jobvr,
                         lapack_int n, double* A, lapack_int lda,
                         double* wr, double* wi, double* vl,
                         lapack_int ldvl, double* vr, lapack_int ldvr);
```

**Paramètres:**
- `jobvl` : calculer vecteurs propres à gauche
- `jobvr` : calculer vecteurs propres à droite
- `wr` : partie réelle des valeurs propres
- `wi` : partie imaginaire des valeurs propres
- `vl` : vecteurs propres à gauche
- `vr` : vecteurs propres à droite

**Attention:** Valeurs propres peuvent être complexes même pour matrices réelles

---

## Catégorie 4: Moindres carrés

### LAPACKE_dgels

**Opération:** Résout min ‖A*x - b‖₂

**Prototype:**
```c
lapack_int LAPACKE_dgels(int matrix_layout, char trans, lapack_int m,
                         lapack_int n, lapack_int nrhs, double* A,
                         lapack_int lda, double* b, lapack_int ldb);
```

**Paramètres:**
- `trans` : 'N' (résout min ‖Ax-b‖), 'T' (résout min ‖A^T*x-b‖)
- `m` : nombre de lignes de A
- `n` : nombre de colonnes de A
- `b` : vecteur b en entrée, solution x en sortie (taille max(m,n))

**Cas d'usage:**
- Régression linéaire
- Ajustement de courbes
- Systèmes surdéterminés (plus d'équations que d'inconnues)

---

## Gestion des erreurs

Toutes les fonctions LAPACKE retournent un entier:

**Valeur de retour:**
- **0** : succès
- **< 0** : le i-ème argument a une valeur illégale (|retour| = i)
- **> 0** : échec de l'algorithme (signification dépend de la fonction)

**Exemple:**
```c
lapack_int info;
info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv, b, ldb);

if (info == 0) {
    printf("Succès\n");
} else if (info < 0) {
    printf("Erreur: argument %d invalide\n", -info);
} else {
    printf("Matrice singulière à l'élément %d\n", info);
}
```

---

## Compilation avec Intel MKL

### Headers nécessaires
```c
#include <stdio.h>
#include "mkl.h"
#include "mkl_lapacke.h"
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

---

## Bonnes pratiques

### 1. Copier les matrices si besoin

La plupart des fonctions LAPACK **modifient** la matrice d'entrée:

```c
// A sera détruite par dgesv
double A_original[n*n];
double A_copy[n*n];
memcpy(A_copy, A_original, n*n*sizeof(double));
LAPACKE_dgesv(..., A_copy, ...);
// A_original est préservée, A_copy contient LU
```

### 2. Vérifier les valeurs de retour

```c
lapack_int info = LAPACKE_dgesv(...);
if (info != 0) {
    // Gérer l'erreur
}
```

### 3. Utiliser les bonnes routines

- Matrice symétrique définie positive → `dposv` (Cholesky)
- Matrice générale → `dgesv` (LU)
- Plusieurs seconds membres avec même A → `dgetrf` + `dgetrs`

### 4. Préférer les driver routines

Les "driver routines" (comme `dgesv`, `dposv`) font tout le travail en un appel. Utiliser les routines séparées seulement si nécessaire.

---

## Tableau récapitulatif

| Problème | Routine | Type de matrice | Complexité |
|----------|---------|-----------------|------------|
| Système linéaire | `dgesv` | Générale | O(n³) |
| Système linéaire | `dposv` | Symétrique définie positive | O(n³/3) |
| Décomposition LU | `dgetrf` | Générale | O(n³) |
| Décomposition QR | `dgeqrf` | Quelconque | O(mn²) |
| Décomposition Cholesky | `dpotrf` | Sym. déf. positive | O(n³/3) |
| SVD | `dgesvd` | Quelconque | O(mn²) |
| Valeurs propres | `dsyev` | Symétrique | O(n³) |
| Valeurs propres | `dgeev` | Générale | O(n³) |
| Moindres carrés | `dgels` | Quelconque | O(mn²) |

---

## Références

- Documentation Intel MKL LAPACK: https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/lapack-routines.html
- LAPACK Users' Guide: http://www.netlib.org/lapack/lug/
- LAPACK Reference: http://www.netlib.org/lapack/

Synthese formatee par Claude AI pour plus de clarte.