/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_EUCLIDEANMIN_C_H
#define EL_LAPACK_EUCLIDEANMIN_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Euclidean minimization
   ====================== */

/* General Linear Model
   -------------------- */
/* Solve min_{X,Y} || Y ||_F subject to D = A X + B Y */

EL_EXPORT ElError ElGLM_s
( ElMatrix_s A, ElMatrix_s B, ElMatrix_s D, ElMatrix_s Y );
EL_EXPORT ElError ElGLM_d
( ElMatrix_d A, ElMatrix_d B, ElMatrix_d D, ElMatrix_d Y );
EL_EXPORT ElError ElGLM_c
( ElMatrix_c A, ElMatrix_c B, ElMatrix_c D, ElMatrix_c Y );
EL_EXPORT ElError ElGLM_z
( ElMatrix_z A, ElMatrix_z B, ElMatrix_z D, ElMatrix_z Y );

EL_EXPORT ElError ElGLMDist_s
( ElDistMatrix_s A, ElDistMatrix_s B, ElDistMatrix_s D, ElDistMatrix_s Y );
EL_EXPORT ElError ElGLMDist_d
( ElDistMatrix_d A, ElDistMatrix_d B, ElDistMatrix_d D, ElDistMatrix_d Y );
EL_EXPORT ElError ElGLMDist_c
( ElDistMatrix_c A, ElDistMatrix_c B, ElDistMatrix_c D, ElDistMatrix_c Y );
EL_EXPORT ElError ElGLMDist_z
( ElDistMatrix_z A, ElDistMatrix_z B, ElDistMatrix_z D, ElDistMatrix_z Y );

/* Least squares
   ------------- */
/* When height(A) >= width(A), solve

     min_X || A X - B ||_F,

   otherwise, solve

     min_X || X ||_F s.t. A X = B. */

EL_EXPORT ElError ElLeastSquares_s
( ElOrientation orientation, ElMatrix_s A, ElConstMatrix_s B, ElMatrix_s X );
EL_EXPORT ElError ElLeastSquares_d
( ElOrientation orientation, ElMatrix_d A, ElConstMatrix_d B, ElMatrix_d X );
EL_EXPORT ElError ElLeastSquares_c
( ElOrientation orientation, ElMatrix_c A, ElConstMatrix_c B, ElMatrix_c X );
EL_EXPORT ElError ElLeastSquares_z
( ElOrientation orientation, ElMatrix_z A, ElConstMatrix_z B, ElMatrix_z X );

EL_EXPORT ElError ElLeastSquaresDist_s
( ElOrientation orientation, 
  ElDistMatrix_s A, ElConstDistMatrix_s B, ElDistMatrix_s X );
EL_EXPORT ElError ElLeastSquaresDist_d
( ElOrientation orientation, 
  ElDistMatrix_d A, ElConstDistMatrix_d B, ElDistMatrix_d X );
EL_EXPORT ElError ElLeastSquaresDist_c
( ElOrientation orientation, 
  ElDistMatrix_c A, ElConstDistMatrix_c B, ElDistMatrix_c X );
EL_EXPORT ElError ElLeastSquaresDist_z
( ElOrientation orientation, 
  ElDistMatrix_z A, ElConstDistMatrix_z B, ElDistMatrix_z X );

EL_EXPORT ElError ElLeastSquaresSparse_s
( ElOrientation orientation, ElSparseMatrix_s A, 
  ElConstMatrix_s B, ElMatrix_s X );
EL_EXPORT ElError ElLeastSquaresSparse_d
( ElOrientation orientation, ElSparseMatrix_d A, 
  ElConstMatrix_d B, ElMatrix_d X );
EL_EXPORT ElError ElLeastSquaresSparse_c
( ElOrientation orientation, ElSparseMatrix_c A, 
  ElConstMatrix_c B, ElMatrix_c X );
EL_EXPORT ElError ElLeastSquaresSparse_z
( ElOrientation orientation, ElSparseMatrix_z A, 
  ElConstMatrix_z B, ElMatrix_z X );

EL_EXPORT ElError ElLeastSquaresDistSparse_s
( ElOrientation orientation,
  ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s B, ElDistMultiVec_s X );
EL_EXPORT ElError ElLeastSquaresDistSparse_d
( ElOrientation orientation,
  ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d B, ElDistMultiVec_d X );
EL_EXPORT ElError ElLeastSquaresDistSparse_c
( ElOrientation orientation,
  ElConstDistSparseMatrix_c A, ElConstDistMultiVec_c B, ElDistMultiVec_c X );
EL_EXPORT ElError ElLeastSquaresDistSparse_z
( ElOrientation orientation,
  ElConstDistSparseMatrix_z A, ElConstDistMultiVec_z B, ElDistMultiVec_z X );

/* Equality-constrained Least Squares
   ---------------------------------- */
/* Solves min_X || A X - C ||_F subject to B X = D */

EL_EXPORT ElError ElLSE_s
( ElMatrix_s A, ElMatrix_s B, ElMatrix_s C, ElMatrix_s D, ElMatrix_s X );
EL_EXPORT ElError ElLSE_d
( ElMatrix_d A, ElMatrix_d B, ElMatrix_d C, ElMatrix_d D, ElMatrix_d X );
EL_EXPORT ElError ElLSE_c
( ElMatrix_c A, ElMatrix_c B, ElMatrix_c C, ElMatrix_c D, ElMatrix_c X );
EL_EXPORT ElError ElLSE_z
( ElMatrix_z A, ElMatrix_z B, ElMatrix_z C, ElMatrix_z D, ElMatrix_z X );

EL_EXPORT ElError ElLSEDist_s
( ElDistMatrix_s A, ElDistMatrix_s B, ElDistMatrix_s C, ElDistMatrix_s D, 
  ElDistMatrix_s X );
EL_EXPORT ElError ElLSEDist_d
( ElDistMatrix_d A, ElDistMatrix_d B, ElDistMatrix_d C, ElDistMatrix_d D, 
  ElDistMatrix_d X );
EL_EXPORT ElError ElLSEDist_c
( ElDistMatrix_c A, ElDistMatrix_c B, ElDistMatrix_c C, ElDistMatrix_c D, 
  ElDistMatrix_c X );
EL_EXPORT ElError ElLSEDist_z
( ElDistMatrix_z A, ElDistMatrix_z B, ElDistMatrix_z C, ElDistMatrix_z D, 
  ElDistMatrix_z X );

/* TODO: Expert version which also returns the residual */

/* Ridge regression
   ---------------- */
/* NOTE: This is simply Tikhonov regularization with Gamma = alpha I */

typedef enum {
  EL_RIDGE_CHOLESKY,
  EL_RIDGE_QR,
  EL_RIDGE_SVD
} ElRidgeAlg;

EL_EXPORT ElError ElRidge_s
( ElConstMatrix_s A, ElConstMatrix_s B, float alpha, ElMatrix_s X,
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidge_d
( ElConstMatrix_d A, ElConstMatrix_d B, double alpha, ElMatrix_d X,
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidge_c
( ElConstMatrix_c A, ElConstMatrix_c B, float alpha, ElMatrix_c X,
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidge_z
( ElConstMatrix_z A, ElConstMatrix_z B, double alpha, ElMatrix_z X,
  ElRidgeAlg alg );

EL_EXPORT ElError ElRidgeDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, float alpha, 
  ElDistMatrix_s X, ElRidgeAlg alg );
EL_EXPORT ElError ElRidgeDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, double alpha, 
  ElDistMatrix_d X, ElRidgeAlg alg );
EL_EXPORT ElError ElRidgeDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, float alpha, 
  ElDistMatrix_c X, ElRidgeAlg alg );
EL_EXPORT ElError ElRidgeDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, double alpha, 
  ElDistMatrix_z X, ElRidgeAlg alg );

EL_EXPORT ElError ElRidgeDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s Y,
  float alpha, ElDistMultiVec_s X );
EL_EXPORT ElError ElRidgeDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d Y,
  double alpha, ElDistMultiVec_d X );
EL_EXPORT ElError ElRidgeDistSparse_c
( ElConstDistSparseMatrix_c A, ElConstDistMultiVec_c Y,
  float alpha, ElDistMultiVec_c X );
EL_EXPORT ElError ElRidgeDistSparse_z
( ElConstDistSparseMatrix_z A, ElConstDistMultiVec_z Y,
  double alpha, ElDistMultiVec_z X );

/* Tikhonov regularization
   ----------------------- */
/* Solve min_X || op(A) X - B ||_2^2 + || Gamma X ||_2^2 */

typedef enum {
  EL_TIKHONOV_CHOLESKY,
  EL_TIKHONOV_QR
} ElTikhonovAlg;

EL_EXPORT ElError ElTikhonov_s
( ElConstMatrix_s A, ElConstMatrix_s B, ElConstMatrix_s Gamma, ElMatrix_s X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonov_d
( ElConstMatrix_d A, ElConstMatrix_d B, ElConstMatrix_d Gamma, ElMatrix_d X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonov_c
( ElConstMatrix_c A, ElConstMatrix_c B, ElConstMatrix_c Gamma, ElMatrix_c X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonov_z
( ElConstMatrix_z A, ElConstMatrix_z B, ElConstMatrix_z Gamma, ElMatrix_z X, 
  ElTikhonovAlg alg );

EL_EXPORT ElError ElTikhonovDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, ElConstDistMatrix_s Gamma, 
  ElDistMatrix_s X, ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonovDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, ElConstDistMatrix_d Gamma, 
  ElDistMatrix_d X, ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonovDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, ElConstDistMatrix_c Gamma, 
  ElDistMatrix_c X, ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonovDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, ElConstDistMatrix_z Gamma, 
  ElDistMatrix_z X, ElTikhonovAlg alg );

EL_EXPORT ElError ElTikhonovDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s Y,
  ElConstDistSparseMatrix_s Gamma, ElDistMultiVec_s X );
EL_EXPORT ElError ElTikhonovDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d Y,
  ElConstDistSparseMatrix_d Gamma, ElDistMultiVec_d X );
EL_EXPORT ElError ElTikhonovDistSparse_c
( ElConstDistSparseMatrix_c A, ElConstDistMultiVec_c Y,
  ElConstDistSparseMatrix_c Gamma, ElDistMultiVec_c X );
EL_EXPORT ElError ElTikhonovDistSparse_z
( ElConstDistSparseMatrix_z A, ElConstDistMultiVec_z Y,
  ElConstDistSparseMatrix_z Gamma, ElDistMultiVec_z X );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_EUCLIDEANMIN_C_H */
