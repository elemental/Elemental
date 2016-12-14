/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <El.h>
using namespace El;

extern "C" {

#define C_PROTO_FIELD(SIG,SIGBASE,T) \
  /* Hyperbolic reflector
     ==================== */ \
  /* Left application
     ---------------- */ \
  ElError ElLeftHyperbolicReflector_ ## SIG \
  ( CREFLECT(T)* chi, ElMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = \
      LeftHyperbolicReflector( *CReflect(chi), *CReflect(x) ) ) } \
  ElError ElLeftHyperbolicReflectorDist_ ## SIG \
  ( CREFLECT(T)* chi, ElDistMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = \
      LeftHyperbolicReflector( *CReflect(chi), *CReflect(x) ) ) } \
  /* Right application
     ----------------- */ \
  ElError ElRightHyperbolicReflector_ ## SIG \
  ( CREFLECT(T)* chi, ElMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = \
      RightHyperbolicReflector( *CReflect(chi), *CReflect(x) ) ) } \
  ElError ElRightHyperbolicReflectorDist_ ## SIG \
  ( CREFLECT(T)* chi, ElDistMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = \
      RightHyperbolicReflector( *CReflect(chi), *CReflect(x) ) ) } \
  /* Householder reflector
     ===================== */ \
  /* Left application
     ---------------- */ \
  ElError ElLeftReflector_ ## SIG \
  ( CREFLECT(T)* chi, ElMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = LeftReflector( *CReflect(chi), *CReflect(x) ) ) } \
  ElError ElLeftReflectorDist_ ## SIG \
  ( CREFLECT(T)* chi, ElDistMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = LeftReflector( *CReflect(chi), *CReflect(x) ) ) } \
  /* Right application
     ----------------- */ \
  ElError ElRightReflector_ ## SIG \
  ( CREFLECT(T)* chi, ElMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = \
      RightReflector( *CReflect(chi), *CReflect(x) ) ) } \
  ElError ElRightReflectorDist_ ## SIG \
  ( CREFLECT(T)* chi, ElDistMatrix_ ## SIG x, CREFLECT(T)* tau ) \
  { EL_TRY( *CReflect(tau) = \
      RightReflector( *CReflect(chi), *CReflect(x) ) ) }

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* Apply packed reflectors
     ======================= */ \
  ElError ElApplyPackedReflectors_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElVerticalOrHorizontal dir, ElForwardOrBackward order, \
    ElInt offset, ElConstMatrix_ ## SIG H, ElConstMatrix_ ## SIG t, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( ApplyPackedReflectors( CReflect(side), CReflect(uplo), \
      CReflect(dir), CReflect(order), UNCONJUGATED, offset, \
      *CReflect(H), *CReflect(t), *CReflect(A) ) ) } \
  ElError ElApplyPackedReflectorsDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElVerticalOrHorizontal dir, ElForwardOrBackward order, \
    ElInt offset, ElConstDistMatrix_ ## SIG H, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG A ) \
  { EL_TRY( ApplyPackedReflectors( CReflect(side), CReflect(uplo), \
      CReflect(dir), CReflect(order), UNCONJUGATED, offset, \
      *CReflect(H), *CReflect(t), *CReflect(A) ) ) } \
  /* Expand packed reflectors
     ======================== */ \
  ElError ElExpandPackedReflectors_ ## SIG \
  ( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, \
    ElInt offset, ElMatrix_ ## SIG H, ElConstMatrix_ ## SIG t ) \
  { EL_TRY( ExpandPackedReflectors( CReflect(uplo), CReflect(dir), \
      UNCONJUGATED, offset, *CReflect(H), *CReflect(t) ) ) } \
  ElError ElExpandPackedReflectorsDist_ ## SIG \
  ( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, \
    ElInt offset, ElDistMatrix_ ## SIG H, ElConstDistMatrix_ ## SIG t ) \
  { EL_TRY( ExpandPackedReflectors( CReflect(uplo), CReflect(dir), \
      UNCONJUGATED, offset, *CReflect(H), *CReflect(t) ) ) }

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Apply packed reflectors
     ======================= */ \
  ElError ElApplyPackedReflectors_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElVerticalOrHorizontal dir, ElForwardOrBackward order, \
    ElConjugation conjugate, ElInt offset, \
    ElConstMatrix_ ## SIG H, ElConstMatrix_ ## SIG t, ElMatrix_ ## SIG A ) \
  { EL_TRY( ApplyPackedReflectors( CReflect(side), CReflect(uplo), \
      CReflect(dir), CReflect(order), CReflect(conjugate), offset, \
      *CReflect(H), *CReflect(t), *CReflect(A) ) ) } \
  ElError ElApplyPackedReflectorsDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElVerticalOrHorizontal dir, ElForwardOrBackward order, \
    ElConjugation conjugate, ElInt offset, \
    ElConstDistMatrix_ ## SIG H, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG A ) \
  { EL_TRY( ApplyPackedReflectors( CReflect(side), CReflect(uplo), \
      CReflect(dir), CReflect(order), CReflect(conjugate), offset, \
      *CReflect(H), *CReflect(t), *CReflect(A) ) ) } \
  /* Expand packed reflectors
     ======================== */ \
  ElError ElExpandPackedReflectors_ ## SIG \
  ( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugate, \
    ElInt offset, ElMatrix_ ## SIG H, ElConstMatrix_ ## SIG t ) \
  { EL_TRY( ExpandPackedReflectors( CReflect(uplo), CReflect(dir), \
      CReflect(conjugate), offset, *CReflect(H), *CReflect(t) ) ) } \
  ElError ElExpandPackedReflectorsDist_ ## SIG \
  ( ElUpperOrLower uplo, ElVerticalOrHorizontal dir, ElConjugation conjugate, \
    ElInt offset, ElDistMatrix_ ## SIG H, ElConstDistMatrix_ ## SIG t ) \
  { EL_TRY( ExpandPackedReflectors( CReflect(uplo), CReflect(dir), \
      CReflect(conjugate), offset, *CReflect(H), *CReflect(t) ) ) }

#define EL_NO_INT_PROTO
#include <El/macros/CInstantiate.h>

} // extern "C"
