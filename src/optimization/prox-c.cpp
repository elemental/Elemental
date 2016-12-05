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

#define C_PROTO_FIELD(SIG,SIGBASE,Field) \
  /* Frobenius-norm prox
     ------------------- */ \
  ElError ElFrobeniusProx_ ## SIG ( ElMatrix_ ## SIG A, Base<Field> rho ) \
  { EL_TRY( FrobeniusProx( *CReflect(A), rho ) ) } \
  ElError ElFrobeniusProxDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, Base<Field> rho ) \
  { EL_TRY( FrobeniusProx( *CReflect(A), rho ) ) } \
  /* Singular-value soft-thresholding
     -------------------------------- */ \
  ElError ElSVT_ ## SIG \
  ( ElMatrix_ ## SIG A, Base<Field> rho, bool relative ) \
  { EL_TRY( SVT( *CReflect(A), rho, relative ) ) } \
  ElError ElSVTDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, Base<Field> rho, bool relative ) \
  { EL_TRY( SVT( *CReflect(A), rho, relative ) ) } \
  /* Soft-thresholding
     ----------------- */ \
  ElError ElSoftThreshold_ ## SIG \
  ( ElMatrix_ ## SIG A, Base<Field> rho, bool relative ) \
  { EL_TRY( SoftThreshold( *CReflect(A), rho, relative ) ) } \
  ElError ElSoftThresholdDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, Base<Field> rho, bool relative ) \
  { EL_TRY( SoftThreshold( *CReflect(A), rho, relative ) ) }

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* Clip
     ---- */ \
  ElError ElLowerClip_ ## SIG ( ElMatrix_ ## SIG X, Real lowerBound ) \
  { EL_TRY( LowerClip( *CReflect(X), lowerBound ) ) } \
  ElError ElLowerClipDist_ ## SIG ( ElDistMatrix_ ## SIG X, Real lowerBound ) \
  { EL_TRY( LowerClip( *CReflect(X), lowerBound ) ) } \
  ElError ElUpperClip_ ## SIG ( ElMatrix_ ## SIG X, Real upperBound ) \
  { EL_TRY( UpperClip( *CReflect(X), upperBound ) ) } \
  ElError ElUpperClipDist_ ## SIG ( ElDistMatrix_ ## SIG X, Real upperBound ) \
  { EL_TRY( UpperClip( *CReflect(X), upperBound ) ) } \
  ElError ElClip_ ## SIG \
  ( ElMatrix_ ## SIG X, Real lowerBound, Real upperBound ) \
  { EL_TRY( Clip( *CReflect(X), lowerBound, upperBound ) ) } \
  ElError ElClipDist_ ## SIG \
  ( ElDistMatrix_ ## SIG X, Real lowerBound, Real upperBound ) \
  { EL_TRY( Clip( *CReflect(X), lowerBound, upperBound ) ) } \
  /* Hinge-loss prox 
     --------------- */ \
  ElError ElHingeLossProx_ ## SIG ( ElMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( HingeLossProx( *CReflect(A), rho ) ) } \
  ElError ElHingeLossProxDist_ ## SIG ( ElDistMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( HingeLossProx( *CReflect(A), rho ) ) } \
  /* Logistic prox
     ------------- */ \
  ElError ElLogisticProx_ ## SIG ( ElMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( LogisticProx( *CReflect(A), rho ) ) } \
  ElError ElLogisticProxDist_ ## SIG ( ElDistMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( LogisticProx( *CReflect(A), rho ) ) }

#define C_PROTO_COMPLEX(SIG,SIGBASE,Field) \
  C_PROTO_FIELD(SIG,SIGBASE,Field) \

#define EL_NO_INT_PROTO
#include <El/macros/CInstantiate.h>

} // extern "C"
