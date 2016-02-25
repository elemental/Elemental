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

ElError ElSignCtrlDefault_s( ElSignCtrl_s* ctrl )
{
    ctrl->maxIts = 100;
    ctrl->tol = 0;
    ctrl->power = 1;
    ctrl->scaling = EL_SIGN_SCALE_FROB;
    ctrl->progress = false;
    return EL_SUCCESS;
}
ElError ElSignCtrlDefault_d( ElSignCtrl_d* ctrl )
{
    ctrl->maxIts = 100;
    ctrl->tol = 0;
    ctrl->power = 1;
    ctrl->scaling = EL_SIGN_SCALE_FROB;
    ctrl->progress = false;
    return EL_SUCCESS;
}

ElError ElSquareRootCtrlDefault_s( ElSquareRootCtrl_s* ctrl )
{
    ctrl->maxIts = 100;
    ctrl->tol = 0;
    ctrl->power = 1;
    ctrl->progress = false;
    return EL_SUCCESS;
}
ElError ElSquareRootCtrlDefault_d( ElSquareRootCtrl_d* ctrl )
{
    ctrl->maxIts = 100;
    ctrl->tol = 0;
    ctrl->power = 1;
    ctrl->progress = false;
    return EL_SUCCESS;
}

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* HermitianFunction [Real] 
     ------------------------ */ \
  ElError ElRealHermitianFunction_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, Base<F> (*funcC)(Base<F>) ) \
  { try { \
      function<Base<F>(Base<F>)> func( funcC ); \
      HermitianFunction( CReflect(uplo), *CReflect(A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElRealHermitianFunctionDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, Base<F> (*funcC)(Base<F>) ) \
  { try { \
      function<Base<F>(Base<F>)> func( funcC ); \
      HermitianFunction( CReflect(uplo), *CReflect(A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Inverse
     ------- */ \
  /* General */ \
  ElError ElInverse_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Inverse( *CReflect(A) ) ) } \
  ElError ElInverseDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Inverse( *CReflect(A) ) ) } \
  ElError ElInverseAfterLUPartialPiv_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstPermutation P ) \
  { EL_TRY( inverse::AfterLUPartialPiv( *CReflect(A), *CReflect(P) ) ) } \
  ElError ElInverseAfterLUPartialPivDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistPermutation P ) \
  { EL_TRY( inverse::AfterLUPartialPiv( *CReflect(A), *CReflect(P) ) ) } \
  /* HPD */ \
  ElError ElHPDInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HPDInverse( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHPDInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HPDInverse( CReflect(uplo), *CReflect(A) ) ) } \
  /* Symmetric */ \
  ElError ElSymmetricInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( SymmetricInverse( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElSymmetricInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( SymmetricInverse( CReflect(uplo), *CReflect(A) ) ) } \
  /* Triangular */ \
  ElError ElTriangularInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_ ## SIG A ) \
  { EL_TRY( TriangularInverse \
            ( CReflect(uplo), CReflect(diag), *CReflect(A) ) ) } \
  ElError ElTriangularInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( TriangularInverse \
            ( CReflect(uplo), CReflect(diag), *CReflect(A) ) ) } \
  /* Pseudoinverse
     ------------- */ \
  /* General */ \
  ElError ElPseudoinverse_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Pseudoinverse( *CReflect(A) ) ) } \
  ElError ElPseudoinverseDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Pseudoinverse( *CReflect(A) ) ) } \
  /* Hermitian */ \
  ElError ElHermitianPseudoinverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPseudoinverse( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianPseudoinverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPseudoinverse( CReflect(uplo), *CReflect(A) ) ) } \
  /* Sign
     ---- */ \
  /* General */ \
  ElError ElSign_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Sign( *CReflect(A) ) ) } \
  ElError ElSignDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Sign( *CReflect(A) ) ) } \
  ElError ElSignDecomp_ ## SIG ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG N ) \
  { EL_TRY( Sign( *CReflect(A), *CReflect(N) ) ) } \
  ElError ElSignDecompDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG N ) \
  { EL_TRY( Sign( *CReflect(A), *CReflect(N) ) ) } \
  /* Hermitian */ \
  ElError ElHermitianSign_ ## SIG ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianSign( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianSignDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianSign( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianSignDecomp_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIG N ) \
  { EL_TRY( HermitianSign( \
      CReflect(uplo), *CReflect(A), *CReflect(N) ) ) } \
  ElError ElHermitianSignDecompDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG N ) \
  { EL_TRY( HermitianSign( CReflect(uplo), *CReflect(A), *CReflect(N) ) ) } \
  /* Square-root
     ----------- */ \
  /* General */ \
  ElError ElSquareRoot_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( SquareRoot( *CReflect(A) ) ) } \
  ElError ElSquareRootDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( SquareRoot( *CReflect(A) ) ) } \
  /* HPSD */ \
  ElError ElHPSDSquareRoot_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HPSDSquareRoot( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHPSDSquareRootDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HPSDSquareRoot( CReflect(uplo), *CReflect(A) ) ) }

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* HermitianFunction [Complex]
     --------------------------- */ \
  ElError ElComplexHermitianFunction_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, CREFLECT(F) (*funcC)(Base<F>) ) \
  { try { \
      auto funcLambda = \
        [&]( Base<F> alpha ) { return CReflect(funcC(alpha)); }; \
      function<F(Base<F>)> func( funcLambda ); \
      HermitianFunction( CReflect(uplo), *CReflect(A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElComplexHermitianFunctionDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    CREFLECT(F) (*funcC)(Base<F>) ) \
  { try { \
      auto funcLambda = \
        [&]( Base<F> alpha ) { return CReflect(funcC(alpha)); }; \
      function<F(Base<F>)> func( funcLambda ); \
      HermitianFunction( CReflect(uplo), *CReflect(A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Hermitian */ \
  ElError ElHermitianInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianInverse( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianInverse( CReflect(uplo), *CReflect(A) ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
