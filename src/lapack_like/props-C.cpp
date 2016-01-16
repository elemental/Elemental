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

#define C_PROTO_BASE(SIG,SIGBASE,T) \
  /* Trace
     ===== */ \
  ElError ElTrace_ ## SIG \
  ( ElConstMatrix_ ## SIG A, CREFLECT(T)* trace ) \
  { EL_TRY( *trace = CReflect(Trace(*CReflect(A))) ) } \
  ElError ElTraceDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, CREFLECT(T)* trace ) \
  { EL_TRY( *trace = CReflect(Trace(*CReflect(A))) ) } \
  /* Norm
     ==== */ \
  /* Max norm
     -------- */ \
  ElError ElMaxNorm_ ## SIG ( ElConstMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = MaxNorm( *CReflect(A) ) ) } \
  ElError ElMaxNormDist_ ## SIG ( ElConstDistMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = MaxNorm( *CReflect(A) ) ) } \
  ElError ElMaxNormSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = MaxNorm( *CReflect(A) ) ) } \
  ElError ElMaxNormDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = MaxNorm( *CReflect(A) ) ) } \
  ElError ElMaxNormDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = MaxNorm( *CReflect(A) ) ) } \
  ElError ElSymmetricMaxNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = SymmetricMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricMaxNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = SymmetricMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricMaxNormSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstSparseMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = SymmetricMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricMaxNormDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistSparseMatrix_ ## SIG A, Base<T>* norm ) \
  { EL_TRY( *norm = SymmetricMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Zero "norm"
     ----------- */ \
  ElError ElZeroNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<T> tol, ElInt* numNonzero ) \
  { EL_TRY( *numNonzero = ZeroNorm( *CReflect(A), tol ) ) } \
  ElError ElZeroNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<T> tol, ElInt* numNonzero) \
  { EL_TRY( *numNonzero = ZeroNorm( *CReflect(A), tol ) ) } \
  ElError ElZeroNormSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, Base<T> tol, ElInt* numNonzero ) \
  { EL_TRY( *numNonzero = ZeroNorm( *CReflect(A), tol ) ) } \
  ElError ElZeroNormDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, Base<T> tol, ElInt* numNonzero ) \
  { EL_TRY( *numNonzero = ZeroNorm( *CReflect(A), tol ) ) }

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  C_PROTO_BASE(SIG,SIGBASE,F) \
  /* Condition number
     ================ */ \
  ElError ElCondition_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElNormType type, Base<F>* cond ) \
  { EL_TRY( *cond = Condition( *CReflect(A), CReflect(type) ) ) } \
  ElError ElConditionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElNormType type, Base<F>* cond ) \
  { EL_TRY( *cond = Condition( *CReflect(A), CReflect(type) ) ) } \
  /* Frobenius
     --------- */ \
  ElError ElFrobeniusCondition_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = FrobeniusCondition( *CReflect(A) ) ) } \
  ElError ElFrobeniusConditionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = FrobeniusCondition( *CReflect(A) ) ) } \
  /* Infinity
     -------- */ \
  ElError ElInfinityCondition_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = InfinityCondition( *CReflect(A) ) ) } \
  ElError ElInfinityConditionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = InfinityCondition( *CReflect(A) ) ) } \
  /* Max
     --- */ \
  ElError ElMaxCondition_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = MaxCondition( *CReflect(A) ) ) } \
  ElError ElMaxConditionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = MaxCondition( *CReflect(A) ) ) } \
  /* One
     --- */ \
  ElError ElOneCondition_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = OneCondition( *CReflect(A) ) ) } \
  ElError ElOneConditionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = OneCondition( *CReflect(A) ) ) } \
  /* Two
     --- */ \
  ElError ElTwoCondition_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = TwoCondition( *CReflect(A) ) ) } \
  ElError ElTwoConditionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* cond ) \
  { EL_TRY( *cond = TwoCondition( *CReflect(A) ) ) } \
  /* Determinant
     =========== */ \
  /* Return the result in a safer, expanded manner
     --------------------------------------------- */ \
  ElError ElSafeDeterminant_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElSafeProduct_ ## SIG* prod ) \
  { EL_TRY( *prod = CReflect(SafeDeterminant(*CReflect(A))) ) } \
  ElError ElSafeDeterminantDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElSafeProduct_ ## SIG* prod ) \
  { EL_TRY( *prod = CReflect(SafeDeterminant(*CReflect(A))) ) } \
  ElError ElSafeHPDDeterminant_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElSafeProduct_ ## SIGBASE* prod ) \
  { EL_TRY( *prod = CReflect(\
      SafeHPDDeterminant(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElSafeHPDDeterminantDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, \
    ElSafeProduct_ ## SIGBASE* prod ) \
  { EL_TRY( *prod = CReflect(\
      SafeHPDDeterminant(CReflect(uplo),*CReflect(A))) ) } \
  /* Return the direct result
     ------------------------ */ \
  ElError ElDeterminant_ ## SIG \
  ( ElConstMatrix_ ## SIG A, CREFLECT(F)* prod ) \
  { EL_TRY( *prod = CReflect(Determinant(*CReflect(A))) ) } \
  ElError ElDeterminantDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, CREFLECT(F)* prod ) \
  { EL_TRY( *prod = CReflect(Determinant(*CReflect(A))) ) } \
  ElError ElHPDDeterminant_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* prod ) \
  { EL_TRY( *prod = HPDDeterminant(CReflect(uplo),*CReflect(A)) ) } \
  ElError ElHPDDeterminantDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* prod ) \
  { EL_TRY( *prod = HPDDeterminant(CReflect(uplo),*CReflect(A)) ) } \
  /* Inertia
     ======= */ \
  ElError ElInertia_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInertiaType* inertia ) \
  { EL_TRY( *inertia = CReflect(Inertia(CReflect(uplo), *CReflect(A))) ) } \
  ElError ElInertiaDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInertiaType* inertia ) \
  { EL_TRY( *inertia = CReflect(Inertia(CReflect(uplo), *CReflect(A))) ) } \
  /* TODO: Expert versions */ \
  /* Norm
     ==== */ \
  ElError ElNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElNormType normType, Base<F>* norm ) \
  { EL_TRY( *norm = Norm( *CReflect(A), CReflect(normType) ) ) } \
  ElError ElNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElNormType normType, Base<F>* norm ) \
  { EL_TRY( *norm = Norm( *CReflect(A), CReflect(normType) ) ) } \
  ElError ElSymmetricNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElNormType normType, \
    Base<F>* norm ) \
  { EL_TRY( *norm = \
      SymmetricNorm( CReflect(uplo), *CReflect(A), CReflect(normType) ) ) } \
  ElError ElSymmetricNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElNormType normType, \
    Base<F>* norm ) \
  { EL_TRY( *norm = \
      SymmetricNorm( CReflect(uplo), *CReflect(A), CReflect(normType) ) ) } \
  /* Entrywise norm
     -------------- */ \
  ElError ElEntrywiseNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = EntrywiseNorm( *CReflect(A), p ) ) } \
  ElError ElEntrywiseNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = EntrywiseNorm( *CReflect(A), p ) ) } \
  ElError ElEntrywiseNormSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = EntrywiseNorm( *CReflect(A), p ) ) } \
  ElError ElEntrywiseNormDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = EntrywiseNorm( *CReflect(A), p ) ) } \
  ElError ElEntrywiseNormDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = EntrywiseNorm( *CReflect(A), p ) ) } \
  ElError ElSymmetricEntrywiseNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElSymmetricEntrywiseNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElSymmetricEntrywiseNormSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstSparseMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElSymmetricEntrywiseNormDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistSparseMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  /* Frobenius norm
     -------------- */ \
  ElError ElFrobeniusNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = FrobeniusNorm( *CReflect(A) ) ) } \
  ElError ElFrobeniusNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = FrobeniusNorm( *CReflect(A) ) ) } \
  ElError ElFrobeniusNormSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = FrobeniusNorm( *CReflect(A) ) ) } \
  ElError ElFrobeniusNormDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = FrobeniusNorm( *CReflect(A) ) ) } \
  ElError ElFrobeniusNormDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = FrobeniusNorm( *CReflect(A) ) ) } \
  ElError ElSymmetricFrobeniusNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricFrobeniusNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricFrobeniusNormSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricFrobeniusNormDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Infinity norm
     ------------- */ \
  ElError ElInfinityNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = InfinityNorm( *CReflect(A) ) ) } \
  ElError ElInfinityNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = InfinityNorm( *CReflect(A) ) ) } \
  ElError ElInfinityNormSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = InfinityNorm( *CReflect(A) ) ) } \
  ElError ElInfinityNormDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = InfinityNorm( *CReflect(A) ) ) } \
  ElError ElSymmetricInfinityNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricInfinityNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricInfinityNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricInfinityNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Ky-Fan norm
     ----------- */ \
  ElError ElKyFanNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt k, Base<F>* norm ) \
  { EL_TRY( *norm = KyFanNorm( *CReflect(A), k ) ) } \
  ElError ElKyFanNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt k, Base<F>* norm ) \
  { EL_TRY( *norm = KyFanNorm( *CReflect(A), k ) ) } \
  ElError ElSymmetricKyFanNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElInt k, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricKyFanNorm( CReflect(uplo), *CReflect(A), k ) ) } \
  ElError ElSymmetricKyFanNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElInt k, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricKyFanNorm( CReflect(uplo), *CReflect(A), k ) ) } \
  /* Ky-Fan-Schatten norm
     -------------------- */ \
  ElError ElKyFanSchattenNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElInt k, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = KyFanSchattenNorm( *CReflect(A), k, p ) ) } \
  ElError ElKyFanSchattenNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElInt k, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = KyFanSchattenNorm( *CReflect(A), k, p ) ) } \
  ElError ElSymmetricKyFanSchattenNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElInt k, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricKyFanSchattenNorm( \
      CReflect(uplo), *CReflect(A), k, p ) ) } \
  ElError ElSymmetricKyFanSchattenNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElInt k, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricKyFanSchattenNorm( \
      CReflect(uplo), *CReflect(A), k, p ) ) } \
  /* Nuclear norm
     ------------ */ \
  ElError ElNuclearNorm_ ## SIG ( ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = NuclearNorm( *CReflect(A) ) ) } \
  ElError ElNuclearNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = NuclearNorm( *CReflect(A) ) ) } \
  ElError ElSymmetricNuclearNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricNuclearNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricNuclearNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricNuclearNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* One norm
     -------- */ \
  ElError ElOneNorm_ ## SIG ( ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = OneNorm( *CReflect(A) ) ) } \
  ElError ElOneNormDist_ ## SIG ( ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = OneNorm( *CReflect(A) ) ) } \
  ElError ElOneNormSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = OneNorm( *CReflect(A) ) ) } \
  ElError ElOneNormDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = OneNorm( *CReflect(A) ) ) } \
  ElError ElSymmetricOneNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricOneNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricOneNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricOneNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Schatten norm
     ------------- */ \
  ElError ElSchattenNorm_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = SchattenNorm( *CReflect(A), p ) ) } \
  ElError ElSchattenNormDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = SchattenNorm( *CReflect(A), p ) ) } \
  ElError ElSymmetricSchattenNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricSchattenNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElSymmetricSchattenNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricSchattenNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  /* Two norm
     -------- */ \
  ElError ElTwoNorm_ ## SIG ( ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = TwoNorm( *CReflect(A) ) ) } \
  ElError ElTwoNormDist_ ## SIG ( ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = TwoNorm( *CReflect(A) ) ) } \
  ElError ElSymmetricTwoNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricTwoNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricTwoNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = SymmetricTwoNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Two norm estimate
     ----------------- */ \
  ElError ElTwoNormEstimate_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F> tol, ElInt maxIts, Base<F>* normEst ) \
  { EL_TRY( *normEst = TwoNormEstimate( *CReflect(A), tol, maxIts ) ) } \
  ElError ElTwoNormEstimateDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F> tol, ElInt maxIts, Base<F>* normEst ) \
  { EL_TRY( *normEst = TwoNormEstimate( *CReflect(A), tol, maxIts ) ) } \
  ElError ElTwoNormEstimateSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt basisSize, Base<F>* normEst ) \
  { EL_TRY( *normEst = TwoNormEstimate( *CReflect(A), basisSize ) ) } \
  ElError ElTwoNormEstimateDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt basisSize, Base<F>* normEst ) \
  { EL_TRY( *normEst = TwoNormEstimate( *CReflect(A), basisSize ) ) } \
  ElError ElSymmetricTwoNormEstimate_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    Base<F> tol, ElInt maxIts, Base<F>* normEst ) \
  { EL_TRY( *normEst = SymmetricTwoNormEstimate( \
      CReflect(uplo), *CReflect(A), tol, maxIts ) ) } \
  ElError ElSymmetricTwoNormEstimateDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, \
    Base<F> tol, ElInt maxIts, Base<F>* normEst ) \
  { EL_TRY( *normEst = SymmetricTwoNormEstimate( \
      CReflect(uplo), *CReflect(A), tol, maxIts ) ) }

#define C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F) \
  /* Norm
     ==== */ \
  ElError ElHermitianNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElNormType normType, \
    Base<F>* norm ) \
  { EL_TRY( *norm = \
      HermitianNorm( CReflect(uplo), *CReflect(A), CReflect(normType) ) ) } \
  ElError ElHermitianNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElNormType normType, \
    Base<F>* norm ) \
  { EL_TRY( *norm = \
      HermitianNorm( CReflect(uplo), *CReflect(A), CReflect(normType) ) ) } \
  /* Entrywise norm
     -------------- */ \
  ElError ElHermitianEntrywiseNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElHermitianEntrywiseNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = HermitianEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElHermitianEntrywiseNormSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstSparseMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = HermitianEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElHermitianEntrywiseNormDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistSparseMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = HermitianEntrywiseNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  /* Frobenius norm
     -------------- */ \
  ElError ElHermitianFrobeniusNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianFrobeniusNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianFrobeniusNormSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianFrobeniusNormDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianFrobeniusNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Infinity norm
     ------------- */ \
  ElError ElHermitianInfinityNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianInfinityNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianInfinityNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianInfinityNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Ky-Fan norm
     ----------- */ \
  ElError ElHermitianKyFanNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElInt k, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianKyFanNorm( CReflect(uplo), *CReflect(A), k ) ) } \
  ElError ElHermitianKyFanNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElInt k, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianKyFanNorm( CReflect(uplo), *CReflect(A), k ) ) } \
  /* Ky-Fan-Schatten norm
     -------------------- */ \
  ElError ElHermitianKyFanSchattenNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElInt k, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = HermitianKyFanSchattenNorm( \
      CReflect(uplo), *CReflect(A), k, p ) ) } \
  ElError ElHermitianKyFanSchattenNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElInt k, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = HermitianKyFanSchattenNorm( \
      CReflect(uplo), *CReflect(A), k, p ) ) } \
  /* Max norm
     -------- */ \
  ElError ElHermitianMaxNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianMaxNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianMaxNormSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianMaxNormDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistSparseMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianMaxNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Nuclear norm
     ------------ */ \
  ElError ElHermitianNuclearNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianNuclearNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianNuclearNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianNuclearNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* One norm
     -------- */ \
  ElError ElHermitianOneNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianOneNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianOneNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianOneNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Schatten norm
     ------------- */ \
  ElError ElHermitianSchattenNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F> p, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianSchattenNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  ElError ElHermitianSchattenNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F> p, \
    Base<F>* norm ) \
  { EL_TRY( *norm = HermitianSchattenNorm( \
      CReflect(uplo), *CReflect(A), p ) ) } \
  /* Two norm
     -------- */ \
  ElError ElHermitianTwoNorm_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianTwoNorm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianTwoNormDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = HermitianTwoNorm( CReflect(uplo), *CReflect(A) ) ) } \
  /* Two norm estimate
     ----------------- */ \
  ElError ElHermitianTwoNormEstimate_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    Base<F> tol, ElInt maxIts, Base<F>* normEst ) \
  { EL_TRY( *normEst = HermitianTwoNormEstimate( \
      CReflect(uplo), *CReflect(A), tol, maxIts ) ) } \
  ElError ElHermitianTwoNormEstimateDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, \
    Base<F> tol, ElInt maxIts, Base<F>* normEst ) \
  { EL_TRY( *normEst = HermitianTwoNormEstimate( \
      CReflect(uplo), *CReflect(A), tol, maxIts ) ) }

#define C_PROTO_INT(SIG,T) C_PROTO_BASE(SIG,SIG,T)

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F)

#include "El/macros/CInstantiate.h"

} // extern "C"
