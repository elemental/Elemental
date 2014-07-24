/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"
using namespace El;

#define DM_CAST(T,A) dynamic_cast<DistMatrix<T>&>(*Reinterpret(A))
#define DM_CAST_CONST(T,A) dynamic_cast<const DistMatrix<T>&>(*Reinterpret(A))

#define DM_STAR_VR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,STAR,VR>&>(*Reinterpret(A))
#define DM_STAR_VR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,STAR,VR>&>(*Reinterpret(A))

#define DM_VC_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VC,STAR>&>(*Reinterpret(A))
#define DM_VC_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VC,STAR>&>(*Reinterpret(A))

#define DM_VR_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VR,STAR>&>(*Reinterpret(A))
#define DM_VR_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VR,STAR>&>(*Reinterpret(A))

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
      std::function<Base<F>(Base<F>)> func( funcC ); \
      HermitianFunction( Reinterpret(uplo), *Reinterpret(A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElRealHermitianFunctionDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, Base<F> (*funcC)(Base<F>) ) \
  { try { \
      std::function<Base<F>(Base<F>)> func( funcC ); \
      HermitianFunction( Reinterpret(uplo), DM_CAST(F,A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Inverse
     ------- */ \
  /* General */ \
  ElError ElInverse_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Inverse( *Reinterpret(A) ) ) } \
  ElError ElInverseDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Inverse( DM_CAST(F,A) ) ) } \
  ElError ElInverseAfterLUPartialPiv_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i p ) \
  { EL_TRY( inverse::AfterLUPartialPiv( *Reinterpret(A), *Reinterpret(p) ) ) } \
  ElError ElInverseAfterLUPartialPivDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p ) \
  { EL_TRY( inverse::AfterLUPartialPiv \
            ( DM_CAST(F,A), DM_VC_STAR_CAST_CONST(Int,p) ) ) } \
  /* HPD */ \
  ElError ElHPDInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HPDInverse( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElHPDInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HPDInverse( Reinterpret(uplo), DM_CAST(F,A) ) ) } \
  /* Symmetric */ \
  ElError ElSymmetricInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( SymmetricInverse( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElSymmetricInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( SymmetricInverse( Reinterpret(uplo), DM_CAST(F,A) ) ) } \
  /* Triangular */ \
  ElError ElTriangularInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_ ## SIG A ) \
  { EL_TRY( TriangularInverse \
            ( Reinterpret(uplo), Reinterpret(diag), *Reinterpret(A) ) ) } \
  ElError ElTriangularInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( TriangularInverse \
            ( Reinterpret(uplo), Reinterpret(diag), DM_CAST(F,A) ) ) } \
  /* Pseudoinverse
     ------------- */ \
  /* General */ \
  ElError ElPseudoinverse_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Pseudoinverse( *Reinterpret(A) ) ) } \
  ElError ElPseudoinverseDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Pseudoinverse( DM_CAST(F,A) ) ) } \
  /* Hermitian */ \
  ElError ElHermitianPseudoinverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPseudoinverse( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElHermitianPseudoinverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPseudoinverse( Reinterpret(uplo), DM_CAST(F,A) ) ) } \
  /* Sign
     ---- */ \
  /* General */ \
  ElError ElSign_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Sign( *Reinterpret(A) ) ) } \
  ElError ElSignDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Sign( DM_CAST(F,A) ) ) } \
  ElError ElSignDecomp_ ## SIG ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG N ) \
  { EL_TRY( Sign( *Reinterpret(A), *Reinterpret(N) ) ) } \
  ElError ElSignDecompDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG N ) \
  { EL_TRY( Sign( DM_CAST(F,A), DM_CAST(F,N) ) ) } \
  /* Hermitian */ \
  ElError ElHermitianSign_ ## SIG ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianSign( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElHermitianSignDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianSign( Reinterpret(uplo), DM_CAST(F,A) ) ) } \
  ElError ElHermitianSignDecomp_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIG N ) \
  { EL_TRY( HermitianSign( \
      Reinterpret(uplo), *Reinterpret(A), *Reinterpret(N) ) ) } \
  ElError ElHermitianSignDecompDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG N ) \
  { EL_TRY( HermitianSign( Reinterpret(uplo), DM_CAST(F,A), DM_CAST(F,N) ) ) } \
  /* Square-root
     ----------- */ \
  /* General */ \
  ElError ElSquareRoot_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( SquareRoot( *Reinterpret(A) ) ) } \
  ElError ElSquareRootDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( SquareRoot( DM_CAST(F,A) ) ) } \
  /* HPSD */ \
  ElError ElHPSDSquareRoot_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HPSDSquareRoot( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElSHPSDquareRootDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HPSDSquareRoot( Reinterpret(uplo), DM_CAST(F,A) ) ) }

#define C_PROTO_REAL(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* HermitianFunction [Complex]
     --------------------------- */ \
  ElError ElComplexHermitianFunction_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, CREFLECT(F) (*funcC)(Base<F>) ) \
  { try { \
      auto funcLambda = \
        [&]( Base<F> alpha ) { return Reinterpret(funcC(alpha)); }; \
      std::function<F(Base<F>)> func( funcLambda ); \
      HermitianFunction( Reinterpret(uplo), *Reinterpret(A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElComplexHermitianFunctionDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    CREFLECT(F) (*funcC)(Base<F>) ) \
  { try { \
      auto funcLambda = \
        [&]( Base<F> alpha ) { return Reinterpret(funcC(alpha)); }; \
      std::function<F(Base<F>)> func( funcLambda ); \
      HermitianFunction( Reinterpret(uplo), DM_CAST(F,A), func ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Hermitian */ \
  ElError ElHermitianInverse_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianInverse( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElHermitianInverseDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianInverse( Reinterpret(uplo), DM_CAST(F,A) ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
