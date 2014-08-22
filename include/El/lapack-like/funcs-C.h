/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_FUNCS_C_H
#define EL_LAPACK_FUNCS_C_H

#include "El/core/DistMatrix-C.h"

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
ElError ElSignCtrlDefault_s( ElSignCtrl_s* ctrl );

typedef struct {
  ElInt maxIts;
  double tol;
  double power;
  ElSignScaling scaling;
  bool progress;
} ElSignCtrl_d;
ElError ElSignCtrlDefault_d( ElSignCtrl_d* ctrl );

typedef struct {
  ElInt maxIts;
  float tol;
  float power;
  bool progress;
} ElSquareRootCtrl_s;
ElError ElSquareRootCtrlDefault_s( ElSquareRootCtrl_s* ctrl );

typedef struct {
  ElInt maxIts;
  double tol;
  double power;
  bool progress;
} ElSquareRootCtrl_d;
ElError ElSquareRootCtrlDefault_d( ElSquareRootCtrl_d* ctrl );

/* Hermitian function
   ================== */

/* Real Hermitian function
   ----------------------- */
ElError ElRealHermitianFunction_s
( ElUpperOrLower uplo, ElMatrix_s A, float (*func)(float) );
ElError ElRealHermitianFunction_d
( ElUpperOrLower uplo, ElMatrix_d A, double (*func)(double) );
ElError ElRealHermitianFunction_c
( ElUpperOrLower uplo, ElMatrix_c A, float (*func)(float) );
ElError ElRealHermitianFunction_z
( ElUpperOrLower uplo, ElMatrix_z A, double (*func)(double) );

ElError ElRealHermitianFunctionDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, float (*func)(float) );
ElError ElRealHermitianFunctionDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, double (*func)(double) );
ElError ElRealHermitianFunctionDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, float (*func)(float) );
ElError ElRealHermitianFunctionDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, double (*func)(double) );

/* Complex Hermitian function
   -------------------------- */
ElError ElComplexHermitianFunction_c
( ElUpperOrLower uplo, ElMatrix_c A, complex_float (*func)(float) );
ElError ElComplexHermitianFunction_z
( ElUpperOrLower uplo, ElMatrix_z A, complex_double (*func)(double) );

ElError ElComplexHermitianFunctionDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, complex_float (*func)(float) );
ElError ElComplexHermitianFunctionDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, complex_double (*func)(double) );

/* Inverse
   ======= */

/* General
   ------- */
ElError ElInverse_s( ElMatrix_s A );
ElError ElInverse_d( ElMatrix_d A );
ElError ElInverse_c( ElMatrix_c A );
ElError ElInverse_z( ElMatrix_z A );

ElError ElInverseDist_s( ElDistMatrix_s A );
ElError ElInverseDist_d( ElDistMatrix_d A );
ElError ElInverseDist_c( ElDistMatrix_c A );
ElError ElInverseDist_z( ElDistMatrix_z A );

/* After LU factorization with partial pivoting
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
ElError ElInverseAfterLUPartialPiv_s( ElMatrix_s A, ElConstMatrix_i p );
ElError ElInverseAfterLUPartialPiv_d( ElMatrix_d A, ElConstMatrix_i p );
ElError ElInverseAfterLUPartialPiv_c( ElMatrix_c A, ElConstMatrix_i p );
ElError ElInverseAfterLUPartialPiv_z( ElMatrix_z A, ElConstMatrix_i p );

ElError ElInverseAfterLUPartialPivDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p );
ElError ElInverseAfterLUPartialPivDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p );
ElError ElInverseAfterLUPartialPivDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p );
ElError ElInverseAfterLUPartialPivDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p );

/* Hermitian Positive-Definite
   --------------------------- */
ElError ElHPDInverse_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElHPDInverse_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElHPDInverse_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElHPDInverse_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElHPDInverseDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElHPDInverseDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElHPDInverseDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElHPDInverseDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Hermitian
   --------- */
ElError ElHermitianInverse_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElHermitianInverse_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElHermitianInverseDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElHermitianInverseDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

/* Symmetric
   --------- */
ElError ElSymmetricInverse_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElSymmetricInverse_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElSymmetricInverse_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElSymmetricInverse_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElSymmetricInverseDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElSymmetricInverseDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElSymmetricInverseDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElSymmetricInverseDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

/* Triangular
   ---------- */
ElError ElTriangularInverse_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_s A );
ElError ElTriangularInverse_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_d A );
ElError ElTriangularInverse_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_c A );
ElError ElTriangularInverse_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_z A );

ElError ElTriangularInverseDist_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_s A );
ElError ElTriangularInverseDist_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_d A );
ElError ElTriangularInverseDist_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_c A );
ElError ElTriangularInverseDist_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_z A );

/* Pseudoinverse
   ============= */

/* General
   ------- */
ElError ElPseudoinverse_s( ElMatrix_s A );
ElError ElPseudoinverse_d( ElMatrix_d A );
ElError ElPseudoinverse_c( ElMatrix_c A );
ElError ElPseudoinverse_z( ElMatrix_z A );

ElError ElPseudoinverseDist_s( ElDistMatrix_s A );
ElError ElPseudoinverseDist_d( ElDistMatrix_d A );
ElError ElPseudoinverseDist_c( ElDistMatrix_c A );
ElError ElPseudoinverseDist_z( ElDistMatrix_z A );

/* TODO: Expert version */

/* Hermitian
   --------- */
ElError ElHermitianPseudoinverse_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElHermitianPseudoinverse_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElHermitianPseudoinverse_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElHermitianPseudoinverse_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElHermitianPseudoinverseDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElHermitianPseudoinverseDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElHermitianPseudoinverseDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElHermitianPseudoinverseDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* Sign
   ==== */

/* General
   ------- */
ElError ElSign_s( ElMatrix_s A );
ElError ElSign_d( ElMatrix_d A );
ElError ElSign_c( ElMatrix_c A );
ElError ElSign_z( ElMatrix_z A );

ElError ElSignDist_s( ElDistMatrix_s A );
ElError ElSignDist_d( ElDistMatrix_d A );
ElError ElSignDist_c( ElDistMatrix_c A );
ElError ElSignDist_z( ElDistMatrix_z A );

/* TODO: Expert version */

ElError ElSignDecomp_s( ElMatrix_s A, ElMatrix_s N );
ElError ElSignDecomp_d( ElMatrix_d A, ElMatrix_d N );
ElError ElSignDecomp_c( ElMatrix_c A, ElMatrix_c N );
ElError ElSignDecomp_z( ElMatrix_z A, ElMatrix_z N );

ElError ElSignDecompDist_s( ElDistMatrix_s A, ElDistMatrix_s N );
ElError ElSignDecompDist_d( ElDistMatrix_d A, ElDistMatrix_d N );
ElError ElSignDecompDist_c( ElDistMatrix_c A, ElDistMatrix_c N );
ElError ElSignDecompDist_z( ElDistMatrix_z A, ElDistMatrix_z N );

/* TODO: Expert version */

/* Hermitian
   --------- */
ElError ElHermitianSign_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElHermitianSign_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElHermitianSign_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElHermitianSign_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElHermitianSignDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElHermitianSignDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElHermitianSignDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElHermitianSignDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

ElError ElHermitianSignDecomp_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s N );
ElError ElHermitianSignDecomp_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d N );
ElError ElHermitianSignDecomp_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_c N );
ElError ElHermitianSignDecomp_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_z N );

ElError ElHermitianSignDecompDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s N );
ElError ElHermitianSignDecompDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d N );
ElError ElHermitianSignDecompDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_c N );
ElError ElHermitianSignDecompDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_z N );

/* TODO: Expert version */

/* Square-root
   =========== */

/* General
   ------- */
ElError ElSquareRoot_s( ElMatrix_s A );
ElError ElSquareRoot_d( ElMatrix_d A );
ElError ElSquareRoot_c( ElMatrix_c A );
ElError ElSquareRoot_z( ElMatrix_z A );

ElError ElSquareRootDist_s( ElDistMatrix_s A );
ElError ElSquareRootDist_d( ElDistMatrix_d A );
ElError ElSquareRootDist_c( ElDistMatrix_c A );
ElError ElSquareRootDist_z( ElDistMatrix_z A );

/* TODO: Expert version */

/* Hermitian Positive Semi-Definite
   -------------------------------- */
ElError ElHPSDSquareRoot_s( ElUpperOrLower uplo, ElMatrix_s A );
ElError ElHPSDSquareRoot_d( ElUpperOrLower uplo, ElMatrix_d A );
ElError ElHPSDSquareRoot_c( ElUpperOrLower uplo, ElMatrix_c A );
ElError ElHPSDSquareRoot_z( ElUpperOrLower uplo, ElMatrix_z A );

ElError ElHPSDSquareRootDist_s( ElUpperOrLower uplo, ElDistMatrix_s A );
ElError ElHPSDSquareRootDist_d( ElUpperOrLower uplo, ElDistMatrix_d A );
ElError ElHPSDSquareRootDist_c( ElUpperOrLower uplo, ElDistMatrix_c A );
ElError ElHPSDSquareRootDist_z( ElUpperOrLower uplo, ElDistMatrix_z A );

/* TODO: Expert version */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_FUNCS_C_H */
