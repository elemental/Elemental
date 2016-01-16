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

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Linear solve
     ------------ */ \
  ElError ElLinearSolve_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( LinearSolve( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElLinearSolveDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( LinearSolve( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElLinearSolveSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( LinearSolve( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElLinearSolveDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG B ) \
  { EL_TRY( LinearSolve( *CReflect(A), *CReflect(B) ) ) } \
  /* Expert versions
     ^^^^^^^^^^^^^^^ */ \
  ElError ElLinearSolveXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG B, \
    ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LinearSolve( *CReflect(A), *CReflect(B), CReflect(ctrl) ) ) } \
  ElError ElLinearSolveXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG B, \
    ElLeastSquaresCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( LinearSolve( *CReflect(A), *CReflect(B), CReflect(ctrl) ) ) } \
  /* HPD solve
     --------- */ \
  ElError ElHPDSolve_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( HPDSolve( CReflect(uplo), CReflect(orientation), \
                      *CReflect(A), *CReflect(B) ) ) } \
  ElError ElHPDSolveDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( HPDSolve( CReflect(uplo), CReflect(orientation), \
                      *CReflect(A), *CReflect(B) ) ) } \
  ElError ElHPDSolveSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( HPDSolve( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElHPDSolveDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG B ) \
  { EL_TRY( HPDSolve( *CReflect(A), *CReflect(B) ) ) } \
  /* Multi-shift Hessenberg solve
     ---------------------------- */ \
  ElError ElMultiShiftHessSolve_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, CREFLECT(F) alpha, \
    ElConstMatrix_ ## SIG H, ElConstMatrix_ ## SIG shifts, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( MultiShiftHessSolve( \
      CReflect(uplo), CReflect(orientation), CReflect(alpha), \
      *CReflect(H), *CReflect(shifts), *CReflect(X) ) ) } \
  ElError ElMultiShiftHessSolveDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, CREFLECT(F) alpha, \
    ElConstDistMatrix_ ## SIG H, ElConstDistMatrix_ ## SIG shifts, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( MultiShiftHessSolve( \
      CReflect(uplo), CReflect(orientation), CReflect(alpha), \
      *CReflect(H), *CReflect(shifts), *CReflect(X) ) ) } \
  /* Symmetric solve
     --------------- */ \
  ElError ElSymmetricSolve_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( SymmetricSolve( CReflect(uplo), CReflect(orientation), \
                            *CReflect(A), *CReflect(B) ) ) } \
  ElError ElSymmetricSolveDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( SymmetricSolve( CReflect(uplo), CReflect(orientation), \
                            *CReflect(A), *CReflect(B) ) ) } \
  ElError ElSymmetricSolveSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG B, bool tryLDL ) \
  { EL_TRY( SymmetricSolve( *CReflect(A), *CReflect(B), tryLDL ) ) } \
  ElError ElSymmetricSolveDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG B, bool tryLDL ) \
  { EL_TRY( SymmetricSolve( *CReflect(A), *CReflect(B), tryLDL ) ) }

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Hermitian solve
     --------------- */ \
  ElError ElHermitianSolve_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( HermitianSolve( CReflect(uplo), CReflect(orientation), \
                            *CReflect(A), *CReflect(B) ) ) } \
  ElError ElHermitianSolveDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( HermitianSolve( CReflect(uplo), CReflect(orientation), \
                            *CReflect(A), *CReflect(B) ) ) } \
  ElError ElHermitianSolveSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG B, \
    bool tryLDL ) \
  { EL_TRY( HermitianSolve( *CReflect(A), *CReflect(B), tryLDL ) ) } \
  ElError ElHermitianSolveDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG B, \
    bool tryLDL ) \
  { EL_TRY( HermitianSolve( *CReflect(A), *CReflect(B), tryLDL ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
