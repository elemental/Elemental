/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

#define C_PROTO(SIG,SIGBASE,F) \
  ElError ElSymmetricSolveDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( SymmetricSolve( *CReflect(A), *CReflect(X) ) ) } \
  ElError ElLeastSquaresDistSparse_ ## SIG \
  ( ElOrientation orientation, \
    ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG X, \
    ElDistMultiVec_ ## SIG Y ) \
  { EL_TRY( LeastSquares( CReflect(orientation), \
      *CReflect(A), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElRidgeDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG X, \
    Base<F> alpha, ElDistMultiVec_ ## SIG Y ) \
  { EL_TRY( Ridge( *CReflect(A), *CReflect(X), alpha, *CReflect(Y) ) ) } \
  ElError ElTikhonovDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG X, \
    ElConstDistSparseMatrix_ ## SIG Gamma, ElDistMultiVec_ ## SIG Y ) \
  { EL_TRY( Tikhonov( *CReflect(A), *CReflect(X), *CReflect(Gamma), \
                      *CReflect(Y) ) ) }

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO(SIG,SIGBASE,F) \
  ElError ElHermitianSolveDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( HermitianSolve( *CReflect(A), *CReflect(X) ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
