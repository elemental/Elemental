/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_FUNCS_C_H
#define EL_LAPACK_FUNCS_C_H

#include <El/core/DistMatrix.h>
#include <El/lapack_like/factor.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EL_SIGN_SCALE_NONE,
  EL_SIGN_SCALE_DET,
  EL_SIGN_SCALE_FROB
} ElSignScaling;

typedef struct {
  ElInt maxIts;
  float tol;
  float power;
  ElSignScaling scaling;
  bool progress;
} ElSignCtrl_s;
EL_EXPORT ElError ElSignCtrlDefault_s( ElSignCtrl_s* ctrl );

typedef struct {
  ElInt maxIts;
  double tol;
  double power;
  ElSignScaling scaling;
  bool progress;
} ElSignCtrl_d;
EL_EXPORT ElError ElSignCtrlDefault_d( ElSignCtrl_d* ctrl );

typedef struct {
  ElInt maxIts;
  float tol;
  float power;
  bool progress;
} ElSquareRootCtrl_s;
EL_EXPORT ElError ElSquareRootCtrlDefault_s( ElSquareRootCtrl_s* ctrl );

typedef struct {
  ElInt maxIts;
  double tol;
  double power;
  bool progress;
} ElSquareRootCtrl_d;
EL_EXPORT ElError ElSquareRootCtrlDefault_d( ElSquareRootCtrl_d* ctrl );

/* Hermitian function
   ================== */

/* Real Hermitian function
   ----------------------- */
EL_EXPORT ElError ElRealHermitianFunction_s
( ElUpperOrLower uplo, ElMatrix_s A, float (*func)(float) );
EL_EXPORT ElError ElRealHermitianFunction_d
( ElUpperOrLower uplo, ElMatrix_d A, double (*func)(double) );
EL_EXPORT ElError ElRealHermitianFunction_c
( ElUpperOrLower uplo, ElMatrix_c A, float (*func)(float) );
EL_EXPORT ElError ElRealHermitianFunction_z
( ElUpperOrLower uplo, ElMatrix_z A, double (*func)(double) );

EL_EXPORT ElError ElRealHermitianFunctionDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, float (*func)(float) );
EL_EXPORT ElError ElRealHermitianFunctionDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, double (*func)(double) );
EL_EXPORT ElError ElRealHermitianFunctionDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, float (*func)(float) );
EL_EXPORT ElError ElRealHermitianFunctionDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, double (*func)(double) );

/* Complex Hermitian function
   -------------------------- */
EL_EXPORT ElError ElComplexHermitianFunction_c
( ElUpperOrLower uplo, ElMatrix_c A, complex_float (*func)(float) );
EL_EXPORT ElError ElComplexHermitianFunction_z
( ElUpperOrLower uplo, ElMatrix_z A, complex_double (*func)(double) );

EL_EXPORT ElError ElComplexHermitianFunctionDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, complex_float (*func)(float) );
EL_EXPORT ElError ElComplexHermitianFunctionDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, complex_double (*func)(double) );

/* Inverse
   ======= */

/* General
   ------- */
EL_EXPORT ElError ElInverse_s( ElMatrix_s A );
EL_EXPORT ElError ElInverse_d( ElMatrix_d A );
EL_EXPORT ElError ElInverse_c( ElMatrix_c A );
EL_EXPORT ElError ElInverse_z( ElMatrix_z A );

EL_EXPORT ElError ElInverseDist_s( ElDistMatrix_s A );
EL_EXPORT ElError ElInverseDist_d( ElDistMatrix_d A );
EL_EXPORT ElError ElInverseDist_c( ElDistMatrix_c A );
EL_EXPORT ElError ElInverseDist_z( ElDistMatrix_z A );

/* After LU factorization with partial pivoting
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElInverseAfterLUPartialPiv_s
( ElMatrix_s A, ElConstPermutation P );
EL_EXPORT ElError ElInverseAfterLUPartialPiv_d
( ElMatrix_d A, ElConstPermutation P );
EL_EXPORT ElError ElInverseAfterLUPartialPiv_c
( ElMatrix_c A, ElConstPermutation P );
EL_EXPORT ElError ElInverseAfterLUPartialPiv_z
( ElMatrix_z A, ElConstPermutation P );

EL_EXPORT ElError ElInverseAfterLUPartialPivDist_s
( ElDistMatrix_s A, ElConstDistPermutation P );
EL_EXPORT ElError ElInverseAfterLUPartialPivDist_d
( ElDistMatrix_d A, ElConstDistPermutation P );
EL_EXPORT ElError ElInverseAfterLUPartialPivDist_c
( ElDistMatrix_c A, ElConstDistPermutation P );
EL_EXPORT ElError ElInverseAfterLUPartialPivDist_z
( ElDistMatrix_z A, ElConstDistPermutation P );

/* Hermitian Positive-Definite
   --------------------------- */
EL_EXPORT ElError ElHPDInverse_s( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElHPDInverse_d( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElHPDInverse_c( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElHPDInverse_z( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElHPDInverseDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElHPDInverseDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElHPDInverseDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElHPDInverseDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Hermitian
   --------- */
EL_EXPORT ElError ElHermitianInverse_c( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElHermitianInverse_z( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElHermitianInverseDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElHermitianInverseDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

/* Symmetric
   --------- */
EL_EXPORT ElError ElSymmetricInverse_s( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElSymmetricInverse_d( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElSymmetricInverse_c( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElSymmetricInverse_z( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElSymmetricInverseDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElSymmetricInverseDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElSymmetricInverseDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElSymmetricInverseDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

/* Triangular
   ---------- */
EL_EXPORT ElError ElTriangularInverse_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_s A );
EL_EXPORT ElError ElTriangularInverse_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_d A );
EL_EXPORT ElError ElTriangularInverse_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_c A );
EL_EXPORT ElError ElTriangularInverse_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_z A );

EL_EXPORT ElError ElTriangularInverseDist_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_s A );
EL_EXPORT ElError ElTriangularInverseDist_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_d A );
EL_EXPORT ElError ElTriangularInverseDist_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_c A );
EL_EXPORT ElError ElTriangularInverseDist_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_z A );

/* Pseudoinverse
   ============= */

/* General
   ------- */
EL_EXPORT ElError ElPseudoinverse_s( ElMatrix_s A );
EL_EXPORT ElError ElPseudoinverse_d( ElMatrix_d A );
EL_EXPORT ElError ElPseudoinverse_c( ElMatrix_c A );
EL_EXPORT ElError ElPseudoinverse_z( ElMatrix_z A );

EL_EXPORT ElError ElPseudoinverseDist_s( ElDistMatrix_s A );
EL_EXPORT ElError ElPseudoinverseDist_d( ElDistMatrix_d A );
EL_EXPORT ElError ElPseudoinverseDist_c( ElDistMatrix_c A );
EL_EXPORT ElError ElPseudoinverseDist_z( ElDistMatrix_z A );

/* TODO: Expert version */

/* Hermitian
   --------- */
EL_EXPORT ElError ElHermitianPseudoinverse_s
( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElHermitianPseudoinverse_d
( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElHermitianPseudoinverse_c
( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElHermitianPseudoinverse_z
( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElHermitianPseudoinverseDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElHermitianPseudoinverseDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElHermitianPseudoinverseDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElHermitianPseudoinverseDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Sign
   ==== */

/* General
   ------- */
EL_EXPORT ElError ElSign_s( ElMatrix_s A );
EL_EXPORT ElError ElSign_d( ElMatrix_d A );
EL_EXPORT ElError ElSign_c( ElMatrix_c A );
EL_EXPORT ElError ElSign_z( ElMatrix_z A );

EL_EXPORT ElError ElSignDist_s( ElDistMatrix_s A );
EL_EXPORT ElError ElSignDist_d( ElDistMatrix_d A );
EL_EXPORT ElError ElSignDist_c( ElDistMatrix_c A );
EL_EXPORT ElError ElSignDist_z( ElDistMatrix_z A );

/* TODO: Expert version */

EL_EXPORT ElError ElSignDecomp_s( ElMatrix_s A, ElMatrix_s N );
EL_EXPORT ElError ElSignDecomp_d( ElMatrix_d A, ElMatrix_d N );
EL_EXPORT ElError ElSignDecomp_c( ElMatrix_c A, ElMatrix_c N );
EL_EXPORT ElError ElSignDecomp_z( ElMatrix_z A, ElMatrix_z N );

EL_EXPORT ElError ElSignDecompDist_s( ElDistMatrix_s A, ElDistMatrix_s N );
EL_EXPORT ElError ElSignDecompDist_d( ElDistMatrix_d A, ElDistMatrix_d N );
EL_EXPORT ElError ElSignDecompDist_c( ElDistMatrix_c A, ElDistMatrix_c N );
EL_EXPORT ElError ElSignDecompDist_z( ElDistMatrix_z A, ElDistMatrix_z N );

/* TODO: Expert version */

/* Hermitian
   --------- */
EL_EXPORT ElError ElHermitianSign_s( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElHermitianSign_d( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElHermitianSign_c( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElHermitianSign_z( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElHermitianSignDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElHermitianSignDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElHermitianSignDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElHermitianSignDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

EL_EXPORT ElError ElHermitianSignDecomp_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s N );
EL_EXPORT ElError ElHermitianSignDecomp_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d N );
EL_EXPORT ElError ElHermitianSignDecomp_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_c N );
EL_EXPORT ElError ElHermitianSignDecomp_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_z N );

EL_EXPORT ElError ElHermitianSignDecompDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s N );
EL_EXPORT ElError ElHermitianSignDecompDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d N );
EL_EXPORT ElError ElHermitianSignDecompDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c N );
EL_EXPORT ElError ElHermitianSignDecompDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z N );

/* TODO: Expert version */

/* Square-root
   =========== */

/* General
   ------- */
EL_EXPORT ElError ElSquareRoot_s( ElMatrix_s A );
EL_EXPORT ElError ElSquareRoot_d( ElMatrix_d A );
EL_EXPORT ElError ElSquareRoot_c( ElMatrix_c A );
EL_EXPORT ElError ElSquareRoot_z( ElMatrix_z A );

EL_EXPORT ElError ElSquareRootDist_s( ElDistMatrix_s A );
EL_EXPORT ElError ElSquareRootDist_d( ElDistMatrix_d A );
EL_EXPORT ElError ElSquareRootDist_c( ElDistMatrix_c A );
EL_EXPORT ElError ElSquareRootDist_z( ElDistMatrix_z A );

/* TODO: Expert version */

/* Hermitian Positive Semi-Definite
   -------------------------------- */
EL_EXPORT ElError ElHPSDSquareRoot_s( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElHPSDSquareRoot_d( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElHPSDSquareRoot_c( ElUpperOrLower uplo, ElMatrix_c A );
EL_EXPORT ElError ElHPSDSquareRoot_z( ElUpperOrLower uplo, ElMatrix_z A );

EL_EXPORT ElError ElHPSDSquareRootDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElHPSDSquareRootDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElHPSDSquareRootDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A );
EL_EXPORT ElError ElHPSDSquareRootDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_FUNCS_C_H */
