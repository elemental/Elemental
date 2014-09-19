/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_SOLVE_C_H
#define EL_LAPACK_SOLVE_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Gaussian Elimination
   ==================== */
ElError ElGaussianElimination_s( ElMatrix_s A, ElMatrix_s B );
ElError ElGaussianElimination_d( ElMatrix_d A, ElMatrix_d B );
ElError ElGaussianElimination_c( ElMatrix_c A, ElMatrix_c B );
ElError ElGaussianElimination_z( ElMatrix_z A, ElMatrix_z B );

ElError ElGaussianEliminationDist_s( ElDistMatrix_s A, ElDistMatrix_s B );
ElError ElGaussianEliminationDist_d( ElDistMatrix_d A, ElDistMatrix_d B );
ElError ElGaussianEliminationDist_c( ElDistMatrix_c A, ElDistMatrix_c B );
ElError ElGaussianEliminationDist_z( ElDistMatrix_z A, ElDistMatrix_z B );

/* General Linear Model
   ==================== */
/* Solve min_{X,Y} || Y ||_F subject to D = A X + B Y */

ElError ElGLM_s( ElMatrix_s A, ElMatrix_s B, ElMatrix_s D, ElMatrix_s Y );
ElError ElGLM_d( ElMatrix_d A, ElMatrix_d B, ElMatrix_d D, ElMatrix_d Y );
ElError ElGLM_c( ElMatrix_c A, ElMatrix_c B, ElMatrix_c D, ElMatrix_c Y );
ElError ElGLM_z( ElMatrix_z A, ElMatrix_z B, ElMatrix_z D, ElMatrix_z Y );

ElError ElGLMDist_s
( ElDistMatrix_s A, ElDistMatrix_s B, ElDistMatrix_s D, ElDistMatrix_s Y );
ElError ElGLMDist_d
( ElDistMatrix_d A, ElDistMatrix_d B, ElDistMatrix_d D, ElDistMatrix_d Y );
ElError ElGLMDist_c
( ElDistMatrix_c A, ElDistMatrix_c B, ElDistMatrix_c D, ElDistMatrix_c Y );
ElError ElGLMDist_z
( ElDistMatrix_z A, ElDistMatrix_z B, ElDistMatrix_z D, ElDistMatrix_z Y );

/* Hermitian solve
   =============== */
ElError ElHermitianSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_c A, ElMatrix_c B );
ElError ElHermitianSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_z A, ElMatrix_z B );

ElError ElHermitianSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_c A, ElDistMatrix_c B );
ElError ElHermitianSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_z A, ElDistMatrix_z B );

/* TODO: Expert version for choosing pivot strategy */

/* HPD solve
   ========= */
ElError ElHPDSolve_s
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_s A, ElMatrix_s B );
ElError ElHPDSolve_d
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_d A, ElMatrix_d B );
ElError ElHPDSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_c A, ElMatrix_c B );
ElError ElHPDSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_z A, ElMatrix_z B );

ElError ElHPDSolveDist_s
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_s A, ElDistMatrix_s B );
ElError ElHPDSolveDist_d
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_d A, ElDistMatrix_d B );
ElError ElHPDSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_c A, ElDistMatrix_c B );
ElError ElHPDSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_z A, ElDistMatrix_z B );

/* Least squares
   ============= */
/* Solves min_X || A X - B ||_F */

ElError ElLeastSquares_s
( ElOrientation orientation, ElMatrix_s A, ElConstMatrix_s B, ElMatrix_s X );
ElError ElLeastSquares_d
( ElOrientation orientation, ElMatrix_d A, ElConstMatrix_d B, ElMatrix_d X );
ElError ElLeastSquares_c
( ElOrientation orientation, ElMatrix_c A, ElConstMatrix_c B, ElMatrix_c X );
ElError ElLeastSquares_z
( ElOrientation orientation, ElMatrix_z A, ElConstMatrix_z B, ElMatrix_z X );

ElError ElLeastSquaresDist_s
( ElOrientation orientation, 
  ElDistMatrix_s A, ElConstDistMatrix_s B, ElDistMatrix_s X );
ElError ElLeastSquaresDist_d
( ElOrientation orientation, 
  ElDistMatrix_d A, ElConstDistMatrix_d B, ElDistMatrix_d X );
ElError ElLeastSquaresDist_c
( ElOrientation orientation, 
  ElDistMatrix_c A, ElConstDistMatrix_c B, ElDistMatrix_c X );
ElError ElLeastSquaresDist_z
( ElOrientation orientation, 
  ElDistMatrix_z A, ElConstDistMatrix_z B, ElDistMatrix_z X );

/* Equality-constrained Least Squares
   ================================== */
/* Solves min_X || A X - C ||_F subject to B X = D */

ElError ElLSE_s
( ElMatrix_s A, ElMatrix_s B, ElMatrix_s C, ElMatrix_s D, ElMatrix_s X );
ElError ElLSE_d
( ElMatrix_d A, ElMatrix_d B, ElMatrix_d C, ElMatrix_d D, ElMatrix_d X );
ElError ElLSE_c
( ElMatrix_c A, ElMatrix_c B, ElMatrix_c C, ElMatrix_c D, ElMatrix_c X );
ElError ElLSE_z
( ElMatrix_z A, ElMatrix_z B, ElMatrix_z C, ElMatrix_z D, ElMatrix_z X );

ElError ElLSEDist_s
( ElDistMatrix_s A, ElDistMatrix_s B, ElDistMatrix_s C, ElDistMatrix_s D, 
  ElDistMatrix_s X );
ElError ElLSEDist_d
( ElDistMatrix_d A, ElDistMatrix_d B, ElDistMatrix_d C, ElDistMatrix_d D, 
  ElDistMatrix_d X );
ElError ElLSEDist_c
( ElDistMatrix_c A, ElDistMatrix_c B, ElDistMatrix_c C, ElDistMatrix_c D, 
  ElDistMatrix_c X );
ElError ElLSEDist_z
( ElDistMatrix_z A, ElDistMatrix_z B, ElDistMatrix_z C, ElDistMatrix_z D, 
  ElDistMatrix_z X );

/* TODO: Expert version which also returns the residual */

/* TODO: Hessenberg solve */

/* Multi-shift Hessenberg solve
   ============================ */
ElError ElMultiShiftHessSolve_s
( ElUpperOrLower uplo, ElOrientation orientation, float alpha,
  ElConstMatrix_s H, ElConstMatrix_s shifts, ElMatrix_s X );
ElError ElMultiShiftHessSolve_d
( ElUpperOrLower uplo, ElOrientation orientation, double alpha,
  ElConstMatrix_d H, ElConstMatrix_d shifts, ElMatrix_d X );
ElError ElMultiShiftHessSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, complex_float alpha,
  ElConstMatrix_c H, ElConstMatrix_c shifts, ElMatrix_c X );
ElError ElMultiShiftHessSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, complex_double alpha,
  ElConstMatrix_z H, ElConstMatrix_z shifts, ElMatrix_z X );

ElError ElMultiShiftHessSolveDist_s
( ElUpperOrLower uplo, ElOrientation orientation, float alpha,
  ElConstDistMatrix_s H, ElConstDistMatrix_s shifts, ElDistMatrix_s X );
ElError ElMultiShiftHessSolveDist_d
( ElUpperOrLower uplo, ElOrientation orientation, double alpha,
  ElConstDistMatrix_d H, ElConstDistMatrix_d shifts, ElDistMatrix_d X );
ElError ElMultiShiftHessSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, complex_float alpha,
  ElConstDistMatrix_c H, ElConstDistMatrix_c shifts, ElDistMatrix_c X );
ElError ElMultiShiftHessSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, complex_double alpha,
  ElConstDistMatrix_z H, ElConstDistMatrix_z shifts, ElDistMatrix_z X );

/* Ridge regression
   ================ */
/* NOTE: This is simply Tikhonov regularization with Gamma = alpha I */

typedef enum {
  EL_RIDGE_CHOLESKY,
  EL_RIDGE_QR,
  EL_RIDGE_SVD
} ElRidgeAlg;

ElError ElRidge_s
( ElConstMatrix_s A, ElConstMatrix_s B, float alpha, ElMatrix_s X,
  ElRidgeAlg alg );
ElError ElRidge_d
( ElConstMatrix_d A, ElConstMatrix_d B, double alpha, ElMatrix_d X,
  ElRidgeAlg alg );
ElError ElRidge_c
( ElConstMatrix_c A, ElConstMatrix_c B, float alpha, ElMatrix_c X,
  ElRidgeAlg alg );
ElError ElRidge_z
( ElConstMatrix_z A, ElConstMatrix_z B, double alpha, ElMatrix_z X,
  ElRidgeAlg alg );

ElError ElRidgeDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float alpha, 
  ElDistMatrix_s X, ElRidgeAlg alg );
ElError ElRidgeDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double alpha, 
  ElDistMatrix_d X, ElRidgeAlg alg );
ElError ElRidgeDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, float alpha, 
  ElDistMatrix_c X, ElRidgeAlg alg );
ElError ElRidgeDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, double alpha, 
  ElDistMatrix_z X, ElRidgeAlg alg );

/* Symmetric solve
   =============== */
ElError ElSymmetricSolve_s
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_s A, ElMatrix_s B );
ElError ElSymmetricSolve_d
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_d A, ElMatrix_d B );
ElError ElSymmetricSolve_c
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_c A, ElMatrix_c B );
ElError ElSymmetricSolve_z
( ElUpperOrLower uplo, ElOrientation orientation, ElMatrix_z A, ElMatrix_z B );

ElError ElSymmetricSolveDist_s
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_s A, ElDistMatrix_s B );
ElError ElSymmetricSolveDist_d
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_d A, ElDistMatrix_d B );
ElError ElSymmetricSolveDist_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_c A, ElDistMatrix_c B );
ElError ElSymmetricSolveDist_z
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElDistMatrix_z A, ElDistMatrix_z B );

/* TODO: Expert version for choosing pivot strategy */

/* Tikhonov regularization
   ======================= */
/* Solve min_X || op(A) X - B ||_2^2 + || Gamma X ||_2^2 */

typedef enum {
  EL_TIKHONOV_CHOLESKY,
  EL_TIKHONOV_QR
} ElTikhonovAlg;

ElError ElTikhonov_s
( ElConstMatrix_s A, ElConstMatrix_s B, ElConstMatrix_s Gamma, ElMatrix_s X, 
  ElTikhonovAlg alg );
ElError ElTikhonov_d
( ElConstMatrix_d A, ElConstMatrix_d B, ElConstMatrix_d Gamma, ElMatrix_d X, 
  ElTikhonovAlg alg );
ElError ElTikhonov_c
( ElConstMatrix_c A, ElConstMatrix_c B, ElConstMatrix_c Gamma, ElMatrix_c X, 
  ElTikhonovAlg alg );
ElError ElTikhonov_z
( ElConstMatrix_z A, ElConstMatrix_z B, ElConstMatrix_z Gamma, ElMatrix_z X, 
  ElTikhonovAlg alg );

ElError ElTikhonovDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, ElConstDistMatrix_s Gamma, 
  ElDistMatrix_s X, ElTikhonovAlg alg );
ElError ElTikhonovDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, ElConstDistMatrix_d Gamma, 
  ElDistMatrix_d X, ElTikhonovAlg alg );
ElError ElTikhonovDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, ElConstDistMatrix_c Gamma, 
  ElDistMatrix_c X, ElTikhonovAlg alg );
ElError ElTikhonovDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, ElConstDistMatrix_z Gamma, 
  ElDistMatrix_z X, ElTikhonovAlg alg );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_SOLVE_C_H */
