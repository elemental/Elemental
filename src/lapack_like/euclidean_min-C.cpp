/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

ElError ElLeastSquaresCtrlDefault_s( ElLeastSquaresCtrl_s* ctrl )
{
    const float eps = limits::Epsilon<float>();
    ctrl->scaleTwoNorm = true;
    ctrl->basisSize = 15;
    ctrl->alpha = Pow(eps,float(0.25));
    ElRegSolveCtrlDefault_s( &ctrl->solveCtrl );
    ctrl->equilibrate = false;
    ctrl->progress = false;
    ctrl->time = false;
    return EL_SUCCESS;
}

ElError ElLeastSquaresCtrlDefault_d( ElLeastSquaresCtrl_d* ctrl )
{
    const double eps = limits::Epsilon<double>();
    ctrl->scaleTwoNorm = true;
    ctrl->basisSize = 15;
    ctrl->alpha = Pow(eps,double(0.25));
    ElRegSolveCtrlDefault_d( &ctrl->solveCtrl );
    ctrl->equilibrate = false;
    ctrl->progress = false;
    ctrl->time = false;
    return EL_SUCCESS;
}

#define C_PROTO(SIG,SIGBASE,F) \
  /* Least squares
     ------------- */ \
  ElError ElLeastSquares_ ## SIG \
  ( ElOrientation orientation, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG X ) \
  { EL_TRY( LeastSquares( CReflect(orientation), *CReflect(A), \
                          *CReflect(B), *CReflect(X) ) ) } \
  ElError ElLeastSquaresDist_ ## SIG \
  ( ElOrientation orientation, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG B, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( LeastSquares( CReflect(orientation), *CReflect(A), \
                          *CReflect(B), *CReflect(X) ) ) } \
  ElError ElLeastSquaresSparse_ ## SIG \
  ( ElOrientation orientation, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG X ) \
  { EL_TRY( LeastSquares( CReflect(orientation), *CReflect(A), \
                          *CReflect(B), *CReflect(X) ) ) } \
  ElError ElLeastSquaresDistSparse_ ## SIG \
  ( ElOrientation orientation, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG B, ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( LeastSquares( CReflect(orientation), \
      *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  /* Expert versions
     ^^^^^^^^^^^^^^^ */ \
  ElError ElLeastSquaresXSparse_ ## SIG \
  ( ElOrientation orientation, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG B,   ElMatrix_ ## SIG X, \
    ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LeastSquares( CReflect(orientation), *CReflect(A), \
                          *CReflect(B), *CReflect(X), CReflect(ctrl) ) ) } \
  ElError ElLeastSquaresXDistSparse_ ## SIG \
  ( ElOrientation orientation, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG B, ElDistMultiVec_ ## SIG X, \
    ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LeastSquares( CReflect(orientation), \
      *CReflect(A), *CReflect(B), *CReflect(X), CReflect(ctrl) ) ) } \
  /* Ridge regression
     ---------------- */ \
  ElError ElRidge_ ## SIG \
  ( ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    Base<F> gamma,           ElMatrix_ ## SIG X, \
    ElRidgeAlg alg ) \
  { EL_TRY( Ridge( CReflect(orientation), \
                   *CReflect(A), *CReflect(B), \
                   gamma,        *CReflect(X), CReflect(alg) ) ) } \
  ElError ElRidgeDist_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    Base<F> gamma,               ElDistMatrix_ ## SIG X, \
    ElRidgeAlg alg ) \
  { EL_TRY( Ridge( CReflect(orientation), \
                   *CReflect(A), *CReflect(B), \
                   gamma,        *CReflect(X), CReflect(alg) ) ) } \
  ElError ElRidgeSparse_ ## SIG \
  ( ElOrientation orientation, \
    ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    Base<F> gamma,                 ElMatrix_ ## SIG X ) \
  { EL_TRY( Ridge( CReflect(orientation), \
                   *CReflect(A), *CReflect(B), \
                   gamma,        *CReflect(X) ) ) } \
  ElError ElRidgeDistSparse_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    Base<F> gamma,                     ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( Ridge( CReflect(orientation), \
                   *CReflect(A), *CReflect(B), \
                   gamma,        *CReflect(X) ) ) } \
  /* Tikhonov regularization
     ----------------------- */ \
  ElError ElTikhonov_ ## SIG \
  ( ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG G, ElMatrix_ ## SIG X, \
    ElTikhonovAlg alg ) \
  { EL_TRY( Tikhonov( CReflect(orientation), \
                      *CReflect(A), *CReflect(B), \
                      *CReflect(G), *CReflect(X), CReflect(alg) ) ) } \
  ElError ElTikhonovDist_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElConstDistMatrix_ ## SIG G, ElDistMatrix_ ## SIG X, \
    ElTikhonovAlg alg ) \
  { EL_TRY( Tikhonov( CReflect(orientation), \
                      *CReflect(A), *CReflect(B), \
                      *CReflect(G), *CReflect(X), CReflect(alg) ) ) } \
  ElError ElTikhonovSparse_ ## SIG \
  ( ElOrientation orientation, \
    ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElConstSparseMatrix_ ## SIG G, ElMatrix_ ## SIG X ) \
  { EL_TRY( Tikhonov( CReflect(orientation), \
                      *CReflect(A), *CReflect(B), \
                      *CReflect(G), *CReflect(X) ) ) } \
  ElError ElTikhonovDistSparse_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElConstDistSparseMatrix_ ## SIG G, ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( Tikhonov( CReflect(orientation), \
                      *CReflect(A), *CReflect(B), \
                      *CReflect(G), *CReflect(X) ) ) } \
  /* Equality-constrained Least Squares
     ---------------------------------- */ \
  ElError ElLSE_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG C, ElConstMatrix_ ## SIG D, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( LSE( *CReflect(A), *CReflect(B), \
                 *CReflect(C), *CReflect(D), \
                 *CReflect(X) ) ) } \
  ElError ElLSEDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElConstDistMatrix_ ## SIG C, ElConstDistMatrix_ ## SIG D, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( LSE( *CReflect(A), *CReflect(B), \
                 *CReflect(C), *CReflect(D), \
                  *CReflect(X) ) ) } \
  ElError ElLSESparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG C,       ElConstMatrix_ ## SIG D, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( LSE( *CReflect(A), *CReflect(B), \
                 *CReflect(C), *CReflect(D), \
                  *CReflect(X) ) ) } \
  ElError ElLSEDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG B, \
    ElConstDistMultiVec_ ## SIG C,     ElConstDistMultiVec_ ## SIG D, \
    ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( LSE( *CReflect(A), *CReflect(B), \
                 *CReflect(C), *CReflect(D), \
                 *CReflect(X) ) ) } \
  /* Expert versions
     ^^^^^^^^^^^^^^^ */ \
  ElError ElLSEXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG C,       ElConstMatrix_ ## SIG D, \
    ElMatrix_ ## SIG X, ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LSE( *CReflect(A), *CReflect(B), \
                 *CReflect(C), *CReflect(D), \
                 *CReflect(X), CReflect(ctrl) ) ) } \
  ElError ElLSEXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG B, \
    ElConstDistMultiVec_ ## SIG C,     ElConstDistMultiVec_ ## SIG D, \
    ElDistMultiVec_ ## SIG X, ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LSE( *CReflect(A), *CReflect(B), \
                 *CReflect(C), *CReflect(D), \
                 *CReflect(X), CReflect(ctrl) ) ) } \
  /* General Linear Model
     -------------------- */ \
  ElError ElGLM_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG D, \
    ElMatrix_ ## SIG X,      ElMatrix_ ## SIG Y ) \
  { EL_TRY( GLM( *CReflect(A), *CReflect(B), \
                 *CReflect(D), \
                 *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElGLMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElConstDistMatrix_ ## SIG D, \
    ElDistMatrix_ ## SIG X,      ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( GLM( *CReflect(A), *CReflect(B), \
                 *CReflect(D), \
                 *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElGLMSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG D, \
    ElMatrix_ ## SIG X,            ElMatrix_ ## SIG Y ) \
  { EL_TRY( GLM( *CReflect(A), *CReflect(B), \
                 *CReflect(D), \
                 *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElGLMDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG B, \
    ElConstDistMultiVec_ ## SIG D, \
    ElDistMultiVec_ ## SIG X,          ElDistMultiVec_ ## SIG Y ) \
  { EL_TRY( GLM( *CReflect(A), *CReflect(B), \
                 *CReflect(D), \
                 *CReflect(X), *CReflect(Y) ) ) } \
  /* Expert versions
     ^^^^^^^^^^^^^^^ */ \
  ElError ElGLMXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG D, \
    ElMatrix_ ## SIG X,            ElMatrix_ ## SIG Y, \
    ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( GLM( *CReflect(A), *CReflect(B), \
                 *CReflect(D), \
                 *CReflect(X), *CReflect(Y), CReflect(ctrl) ) ) } \
  ElError ElGLMXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG B, \
    ElConstDistMultiVec_ ## SIG D, \
    ElDistMultiVec_ ## SIG X,          ElDistMultiVec_ ## SIG Y, \
    ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( GLM( *CReflect(A), *CReflect(B), \
                 *CReflect(D), \
                 *CReflect(X), *CReflect(Y), CReflect(ctrl) ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
