/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_EUCLIDEANMIN_C_H
#define EL_LAPACK_EUCLIDEANMIN_C_H

#include <El/core/DistMatrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Least squares
   ============= */
/* When height(A) >= width(A), solve

     min_X || A X - B ||_F,

   otherwise, solve

     min_X || X ||_F s.t. A X = B. */

EL_EXPORT ElError ElLeastSquares_s
( ElOrientation orientation, ElConstMatrix_s A, 
  ElConstMatrix_s B, ElMatrix_s X );
EL_EXPORT ElError ElLeastSquares_d
( ElOrientation orientation, ElConstMatrix_d A, 
  ElConstMatrix_d B, ElMatrix_d X );
EL_EXPORT ElError ElLeastSquares_c
( ElOrientation orientation, ElConstMatrix_c A, 
  ElConstMatrix_c B, ElMatrix_c X );
EL_EXPORT ElError ElLeastSquares_z
( ElOrientation orientation, ElConstMatrix_z A, 
  ElConstMatrix_z B, ElMatrix_z X );

EL_EXPORT ElError ElLeastSquaresDist_s
( ElOrientation orientation, ElConstDistMatrix_s A, 
  ElConstDistMatrix_s B, ElDistMatrix_s X );
EL_EXPORT ElError ElLeastSquaresDist_d
( ElOrientation orientation, ElConstDistMatrix_d A, 
  ElConstDistMatrix_d B, ElDistMatrix_d X );
EL_EXPORT ElError ElLeastSquaresDist_c
( ElOrientation orientation, ElConstDistMatrix_c A, 
  ElConstDistMatrix_c B, ElDistMatrix_c X );
EL_EXPORT ElError ElLeastSquaresDist_z
( ElOrientation orientation, ElConstDistMatrix_z A, 
  ElConstDistMatrix_z B, ElDistMatrix_z X );

EL_EXPORT ElError ElLeastSquaresSparse_s
( ElOrientation orientation, ElConstSparseMatrix_s A, 
  ElConstMatrix_s B, ElMatrix_s X );
EL_EXPORT ElError ElLeastSquaresSparse_d
( ElOrientation orientation, ElConstSparseMatrix_d A, 
  ElConstMatrix_d B, ElMatrix_d X );
EL_EXPORT ElError ElLeastSquaresSparse_c
( ElOrientation orientation, ElConstSparseMatrix_c A, 
  ElConstMatrix_c B, ElMatrix_c X );
EL_EXPORT ElError ElLeastSquaresSparse_z
( ElOrientation orientation, ElConstSparseMatrix_z A, 
  ElConstMatrix_z B, ElMatrix_z X );

EL_EXPORT ElError ElLeastSquaresDistSparse_s
( ElOrientation orientation, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s B, ElDistMultiVec_s X );
EL_EXPORT ElError ElLeastSquaresDistSparse_d
( ElOrientation orientation, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d B, ElDistMultiVec_d X );
EL_EXPORT ElError ElLeastSquaresDistSparse_c
( ElOrientation orientation, ElConstDistSparseMatrix_c A, 
  ElConstDistMultiVec_c B, ElDistMultiVec_c X );
EL_EXPORT ElError ElLeastSquaresDistSparse_z
( ElOrientation orientation, ElConstDistSparseMatrix_z A, 
  ElConstDistMultiVec_z B, ElDistMultiVec_z X );

/* Expert versions
   --------------- */
typedef struct {
  bool scaleTwoNorm;
  ElInt basisSize;
  float alpha;
  float reg0Tmp;
  float reg0Perm;
  float reg1Tmp;
  float reg1Perm;
  ElRegSolveCtrl_s solveCtrl;
  bool equilibrate;
  bool progress;
  bool time;
} ElLeastSquaresCtrl_s;
EL_EXPORT ElError ElLeastSquaresCtrlDefault_s( ElLeastSquaresCtrl_s* ctrl );

typedef struct {
  bool scaleTwoNorm;
  ElInt basisSize;
  double alpha;
  double reg0Tmp;
  double reg0Perm;
  double reg1Tmp;
  double reg1Perm;
  ElRegSolveCtrl_d solveCtrl;
  bool equilibrate;
  bool progress;
  bool time;
} ElLeastSquaresCtrl_d;
EL_EXPORT ElError ElLeastSquaresCtrlDefault_d( ElLeastSquaresCtrl_d* ctrl );

EL_EXPORT ElError ElLeastSquaresXSparse_s
( ElOrientation orientation, ElConstSparseMatrix_s A, 
  ElConstMatrix_s B, ElMatrix_s X, ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLeastSquaresXSparse_d
( ElOrientation orientation, ElConstSparseMatrix_d A, 
  ElConstMatrix_d B, ElMatrix_d X, ElLeastSquaresCtrl_d ctrl );
EL_EXPORT ElError ElLeastSquaresXSparse_c
( ElOrientation orientation, ElConstSparseMatrix_c A, 
  ElConstMatrix_c B, ElMatrix_c X, ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLeastSquaresXSparse_z
( ElOrientation orientation, ElConstSparseMatrix_z A, 
  ElConstMatrix_z B, ElMatrix_z X, ElLeastSquaresCtrl_d ctrl );

EL_EXPORT ElError ElLeastSquaresXDistSparse_s
( ElOrientation orientation, ElConstDistSparseMatrix_s A, 
  ElConstDistMultiVec_s B, ElDistMultiVec_s X,
  ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLeastSquaresXDistSparse_d
( ElOrientation orientation, ElConstDistSparseMatrix_d A, 
  ElConstDistMultiVec_d B, ElDistMultiVec_d X,
  ElLeastSquaresCtrl_d ctrl );
EL_EXPORT ElError ElLeastSquaresXDistSparse_c
( ElOrientation orientation, ElConstDistSparseMatrix_c A, 
  ElConstDistMultiVec_c B, ElDistMultiVec_c X,
  ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLeastSquaresXDistSparse_z
( ElOrientation orientation, ElConstDistSparseMatrix_z A, 
  ElConstDistMultiVec_z B, ElDistMultiVec_z X,
  ElLeastSquaresCtrl_d ctrl );

/* Ridge regression
   ================ */
/* Ridge regression is a special case of Tikhonov regularization with 
   the regularization matrix equal to gamma I */

typedef enum {
  EL_RIDGE_CHOLESKY,
  EL_RIDGE_QR,
  EL_RIDGE_SVD
} ElRidgeAlg;

EL_EXPORT ElError ElRidge_s
( ElOrientation orientation,
  ElConstMatrix_s A, ElConstMatrix_s B, 
  float gamma,       ElMatrix_s X,
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidge_d
( ElOrientation orientation,
  ElConstMatrix_d A, ElConstMatrix_d B, 
  double gamma,      ElMatrix_d X,
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidge_c
( ElOrientation orientation,
  ElConstMatrix_c A, ElConstMatrix_c B, 
  float gamma,       ElMatrix_c X,
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidge_z
( ElOrientation orientation,
  ElConstMatrix_z A, ElConstMatrix_z B, 
  double gamma,      ElMatrix_z X,
  ElRidgeAlg alg );

EL_EXPORT ElError ElRidgeDist_s
( ElOrientation orientation,
  ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  float gamma,           ElDistMatrix_s X, 
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidgeDist_d
( ElOrientation orientation,
  ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  double gamma,          ElDistMatrix_d X, 
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidgeDist_c
( ElOrientation orientation,
  ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  float gamma,           ElDistMatrix_c X, 
  ElRidgeAlg alg );
EL_EXPORT ElError ElRidgeDist_z
( ElOrientation orientation,
  ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  double gamma,          ElDistMatrix_z X, 
  ElRidgeAlg alg );

EL_EXPORT ElError ElRidgeSparse_s
( ElOrientation orientation,
  ElConstSparseMatrix_s A, ElConstMatrix_s B,
  float gamma,             ElMatrix_s X );
EL_EXPORT ElError ElRidgeSparse_d
( ElOrientation orientation,
  ElConstSparseMatrix_d A, ElConstMatrix_d B,
  double gamma,            ElMatrix_d X );
EL_EXPORT ElError ElRidgeSparse_c
( ElOrientation orientation,
  ElConstSparseMatrix_c A, ElConstMatrix_c B,
  float gamma,             ElMatrix_c X );
EL_EXPORT ElError ElRidgeSparse_z
( ElOrientation orientation,
  ElConstSparseMatrix_z A, ElConstMatrix_z B,
  double gamma,            ElMatrix_z X );

EL_EXPORT ElError ElRidgeDistSparse_s
( ElOrientation orientation,
  ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s B,
  float gamma,                 ElDistMultiVec_s X );
EL_EXPORT ElError ElRidgeDistSparse_d
( ElOrientation orientation,
  ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d B,
  double gamma,                ElDistMultiVec_d X );
EL_EXPORT ElError ElRidgeDistSparse_c
( ElOrientation orientation,
  ElConstDistSparseMatrix_c A, ElConstDistMultiVec_c B,
  float gamma,                 ElDistMultiVec_c X );
EL_EXPORT ElError ElRidgeDistSparse_z
( ElOrientation orientation,
  ElConstDistSparseMatrix_z A, ElConstDistMultiVec_z B,
  double gamma,                ElDistMultiVec_z X );

/* Tikhonov regularization
   ======================= */
/* Defining W = op(A), where op(A) is either A, A^T, or A^H, Tikhonov
   regularization involves the solution of either

   Regularized Least Squares:

     min_X || [W; G] X - [B; 0] ||_F^2

   or

   Regularized Minimum Length:

     min_{X,S} || [X, S] ||_F^
     s.t.      [W, G] [X; S] = B.
*/

typedef enum {
  EL_TIKHONOV_CHOLESKY,
  EL_TIKHONOV_QR
} ElTikhonovAlg;

EL_EXPORT ElError ElTikhonov_s
( ElOrientation orientation,
  ElConstMatrix_s A, ElConstMatrix_s B, 
  ElConstMatrix_s G, ElMatrix_s X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonov_d
( ElOrientation orientation,
  ElConstMatrix_d A, ElConstMatrix_d B, 
  ElConstMatrix_d G, ElMatrix_d X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonov_c
( ElOrientation orientation,
  ElConstMatrix_c A, ElConstMatrix_c B, 
  ElConstMatrix_c G, ElMatrix_c X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonov_z
( ElOrientation orientation,
  ElConstMatrix_z A, ElConstMatrix_z B, 
  ElConstMatrix_z G, ElMatrix_z X, 
  ElTikhonovAlg alg );

EL_EXPORT ElError ElTikhonovDist_s
( ElOrientation orientation,
  ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  ElConstDistMatrix_s G, ElDistMatrix_s X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonovDist_d
( ElOrientation orientation,
  ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  ElConstDistMatrix_d G, ElDistMatrix_d X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonovDist_c
( ElOrientation orientation,
  ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  ElConstDistMatrix_c G, ElDistMatrix_c X, 
  ElTikhonovAlg alg );
EL_EXPORT ElError ElTikhonovDist_z
( ElOrientation orientation,
  ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  ElConstDistMatrix_z G, ElDistMatrix_z X, 
  ElTikhonovAlg alg );

EL_EXPORT ElError ElTikhonovSparse_s
( ElOrientation orientation,
  ElConstSparseMatrix_s A, ElConstMatrix_s B,
  ElConstSparseMatrix_s G, ElMatrix_s X );
EL_EXPORT ElError ElTikhonovSparse_d
( ElOrientation orientation,
  ElConstSparseMatrix_d A, ElConstMatrix_d B,
  ElConstSparseMatrix_d G, ElMatrix_d X );
EL_EXPORT ElError ElTikhonovSparse_c
( ElOrientation orientation,
  ElConstSparseMatrix_c A, ElConstMatrix_c B,
  ElConstSparseMatrix_c G, ElMatrix_c X );
EL_EXPORT ElError ElTikhonovSparse_z
( ElOrientation orientation,
  ElConstSparseMatrix_z A, ElConstMatrix_z B,
  ElConstSparseMatrix_z G, ElMatrix_z X );

EL_EXPORT ElError ElTikhonovDistSparse_s
( ElOrientation orientation,
  ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s B,
  ElConstDistSparseMatrix_s G, ElDistMultiVec_s X );
EL_EXPORT ElError ElTikhonovDistSparse_d
( ElOrientation orientation,
  ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d B,
  ElConstDistSparseMatrix_d G, ElDistMultiVec_d X );
EL_EXPORT ElError ElTikhonovDistSparse_c
( ElOrientation orientation,
  ElConstDistSparseMatrix_c A, ElConstDistMultiVec_c B,
  ElConstDistSparseMatrix_c G, ElDistMultiVec_c X );
EL_EXPORT ElError ElTikhonovDistSparse_z
( ElOrientation orientation,
  ElConstDistSparseMatrix_z A, ElConstDistMultiVec_z B,
  ElConstDistSparseMatrix_z G, ElDistMultiVec_z X );

/* Equality-constrained Least Squares
   ================================== */
/* Solves min_X || A X - C ||_F subject to B X = D */

EL_EXPORT ElError ElLSE_s
( ElConstMatrix_s A, ElConstMatrix_s B, 
  ElConstMatrix_s C, ElConstMatrix_s D, 
  ElMatrix_s X );
EL_EXPORT ElError ElLSE_d
( ElConstMatrix_d A, ElConstMatrix_d B, 
  ElConstMatrix_d C, ElConstMatrix_d D, 
  ElMatrix_d X );
EL_EXPORT ElError ElLSE_c
( ElConstMatrix_c A, ElConstMatrix_c B, 
  ElConstMatrix_c C, ElConstMatrix_c D, 
  ElMatrix_c X );
EL_EXPORT ElError ElLSE_z
( ElConstMatrix_z A, ElConstMatrix_z B, 
  ElConstMatrix_z C, ElConstMatrix_z D, 
  ElMatrix_z X );

EL_EXPORT ElError ElLSEDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  ElConstDistMatrix_s C, ElConstDistMatrix_s D, 
  ElDistMatrix_s X );
EL_EXPORT ElError ElLSEDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  ElConstDistMatrix_d C, ElConstDistMatrix_d D, 
  ElDistMatrix_d X );
EL_EXPORT ElError ElLSEDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  ElConstDistMatrix_c C, ElConstDistMatrix_c D, 
  ElDistMatrix_c X );
EL_EXPORT ElError ElLSEDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  ElConstDistMatrix_z C, ElConstDistMatrix_z D, 
  ElDistMatrix_z X );

EL_EXPORT ElError ElLSESparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s B,
  ElConstMatrix_s C,       ElConstMatrix_s D,
  ElMatrix_s X );
EL_EXPORT ElError ElLSESparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d B,
  ElConstMatrix_d C,       ElConstMatrix_d D,
  ElMatrix_d X );
EL_EXPORT ElError ElLSESparse_c
( ElConstSparseMatrix_c A, ElConstSparseMatrix_c B,
  ElConstMatrix_c C,       ElConstMatrix_c D,
  ElMatrix_c X );
EL_EXPORT ElError ElLSESparse_z
( ElConstSparseMatrix_z A, ElConstSparseMatrix_z B,
  ElConstMatrix_z C,       ElConstMatrix_z D,
  ElMatrix_z X );

EL_EXPORT ElError ElLSEDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s B,
  ElConstDistMultiVec_s C,     ElConstDistMultiVec_s D,
  ElDistMultiVec_s X );
EL_EXPORT ElError ElLSEDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d B,
  ElConstDistMultiVec_d C,     ElConstDistMultiVec_d D,
  ElDistMultiVec_d X );
EL_EXPORT ElError ElLSEDistSparse_c
( ElConstDistSparseMatrix_c A, ElConstDistSparseMatrix_c B,
  ElConstDistMultiVec_c C,     ElConstDistMultiVec_c D,
  ElDistMultiVec_c X );
EL_EXPORT ElError ElLSEDistSparse_z
( ElConstDistSparseMatrix_z A, ElConstDistSparseMatrix_z B,
  ElConstDistMultiVec_z C,     ElConstDistMultiVec_z D,
  ElDistMultiVec_z X );

/* Expert versions
   --------------- */
/* TODO: Dense expert version which also return the residual */
EL_EXPORT ElError ElLSEXSparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s B,
  ElConstMatrix_s C,       ElConstMatrix_s D,
  ElMatrix_s X, ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLSEXSparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d B,
  ElConstMatrix_d C,       ElConstMatrix_d D,
  ElMatrix_d X, ElLeastSquaresCtrl_d ctrl );
EL_EXPORT ElError ElLSEXSparse_c
( ElConstSparseMatrix_c A, ElConstSparseMatrix_c B,
  ElConstMatrix_c C,       ElConstMatrix_c D,
  ElMatrix_c X, ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLSEXSparse_z
( ElConstSparseMatrix_z A, ElConstSparseMatrix_z B,
  ElConstMatrix_z C,       ElConstMatrix_z D,
  ElMatrix_z X, ElLeastSquaresCtrl_d ctrl );

EL_EXPORT ElError ElLSEXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s B,
  ElConstDistMultiVec_s C,     ElConstDistMultiVec_s D,
  ElDistMultiVec_s X, ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLSEXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d B,
  ElConstDistMultiVec_d C,     ElConstDistMultiVec_d D,
  ElDistMultiVec_d X, ElLeastSquaresCtrl_d ctrl );
EL_EXPORT ElError ElLSEXDistSparse_c
( ElConstDistSparseMatrix_c A, ElConstDistSparseMatrix_c B,
  ElConstDistMultiVec_c C,     ElConstDistMultiVec_c D,
  ElDistMultiVec_c X, ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElLSEXDistSparse_z
( ElConstDistSparseMatrix_z A, ElConstDistSparseMatrix_z B,
  ElConstDistMultiVec_z C,     ElConstDistMultiVec_z D,
  ElDistMultiVec_z X, ElLeastSquaresCtrl_d ctrl );

/* General Linear Model
   ==================== */
/* Solve min_{X,Y} || Y ||_F subject to A X + B Y = D */

EL_EXPORT ElError ElGLM_s
( ElConstMatrix_s A, ElConstMatrix_s B, 
  ElConstMatrix_s D, 
  ElMatrix_s X,      ElMatrix_s Y );
EL_EXPORT ElError ElGLM_d
( ElConstMatrix_d A, ElConstMatrix_d B, 
  ElConstMatrix_d D, 
  ElMatrix_d X,      ElMatrix_d Y );
EL_EXPORT ElError ElGLM_c
( ElConstMatrix_c A, ElConstMatrix_c B, 
  ElConstMatrix_c D, 
  ElMatrix_c X,      ElMatrix_c Y );
EL_EXPORT ElError ElGLM_z
( ElConstMatrix_z A, ElConstMatrix_z B, 
  ElConstMatrix_z D, 
  ElMatrix_z X,      ElMatrix_z Y );

EL_EXPORT ElError ElGLMDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  ElConstDistMatrix_s D, 
  ElDistMatrix_s X,      ElDistMatrix_s Y );
EL_EXPORT ElError ElGLMDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  ElConstDistMatrix_d D, 
  ElDistMatrix_d X,      ElDistMatrix_d Y );
EL_EXPORT ElError ElGLMDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  ElConstDistMatrix_c D, 
  ElDistMatrix_c X,      ElDistMatrix_c Y );
EL_EXPORT ElError ElGLMDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  ElConstDistMatrix_z D, 
  ElDistMatrix_z X,      ElDistMatrix_z Y );

EL_EXPORT ElError ElGLMSparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s B,
  ElConstMatrix_s D,
  ElMatrix_s X,            ElMatrix_s Y );
EL_EXPORT ElError ElGLMSparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d B,
  ElConstMatrix_d D,
  ElMatrix_d X,            ElMatrix_d Y );
EL_EXPORT ElError ElGLMSparse_c
( ElConstSparseMatrix_c A, ElConstSparseMatrix_c B,
  ElConstMatrix_c D,
  ElMatrix_c X,            ElMatrix_c Y );
EL_EXPORT ElError ElGLMSparse_z
( ElConstSparseMatrix_z A, ElConstSparseMatrix_z B,
  ElConstMatrix_z D,
  ElMatrix_z X,            ElMatrix_z Y );

EL_EXPORT ElError ElGLMDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s B,
  ElConstDistMultiVec_s D,
  ElDistMultiVec_s X,          ElDistMultiVec_s Y );
EL_EXPORT ElError ElGLMDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d B,
  ElConstDistMultiVec_d D,
  ElDistMultiVec_d X,          ElDistMultiVec_d Y );
EL_EXPORT ElError ElGLMDistSparse_c
( ElConstDistSparseMatrix_c A, ElConstDistSparseMatrix_c B,
  ElConstDistMultiVec_c D,
  ElDistMultiVec_c X,          ElDistMultiVec_c Y );
EL_EXPORT ElError ElGLMDistSparse_z
( ElConstDistSparseMatrix_z A, ElConstDistSparseMatrix_z B,
  ElConstDistMultiVec_z D,
  ElDistMultiVec_z X,          ElDistMultiVec_z Y );

/* Expert versions
   --------------- */
EL_EXPORT ElError ElGLMXSparse_s
( ElConstSparseMatrix_s A, ElConstSparseMatrix_s B,
  ElConstMatrix_s D,
  ElMatrix_s X,            ElMatrix_s Y,
  ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElGLMXSparse_d
( ElConstSparseMatrix_d A, ElConstSparseMatrix_d B,
  ElConstMatrix_d D,
  ElMatrix_d X,            ElMatrix_d Y,
  ElLeastSquaresCtrl_d ctrl );
EL_EXPORT ElError ElGLMXSparse_c
( ElConstSparseMatrix_c A, ElConstSparseMatrix_c B,
  ElConstMatrix_c D,
  ElMatrix_c X,            ElMatrix_c Y,
  ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElGLMXSparse_z
( ElConstSparseMatrix_z A, ElConstSparseMatrix_z B,
  ElConstMatrix_z D,
  ElMatrix_z X,            ElMatrix_z Y,
  ElLeastSquaresCtrl_d ctrl );

EL_EXPORT ElError ElGLMXDistSparse_s
( ElConstDistSparseMatrix_s A, ElConstDistSparseMatrix_s B,
  ElConstDistMultiVec_s D,
  ElDistMultiVec_s X,          ElDistMultiVec_s Y,
  ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElGLMXDistSparse_d
( ElConstDistSparseMatrix_d A, ElConstDistSparseMatrix_d B,
  ElConstDistMultiVec_d D,
  ElDistMultiVec_d X,          ElDistMultiVec_d Y,
  ElLeastSquaresCtrl_d ctrl );
EL_EXPORT ElError ElGLMXDistSparse_c
( ElConstDistSparseMatrix_c A, ElConstDistSparseMatrix_c B,
  ElConstDistMultiVec_c D,
  ElDistMultiVec_c X,          ElDistMultiVec_c Y,
  ElLeastSquaresCtrl_s ctrl );
EL_EXPORT ElError ElGLMXDistSparse_z
( ElConstDistSparseMatrix_z A, ElConstDistSparseMatrix_z B,
  ElConstDistMultiVec_z D,
  ElDistMultiVec_z X,          ElDistMultiVec_z Y,
  ElLeastSquaresCtrl_d ctrl );

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus
#include <El/lapack_like/euclidean_min/CReflect.hpp>
#endif

#endif /* ifndef EL_LAPACK_EUCLIDEANMIN_C_H */
