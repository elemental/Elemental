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

ElError ElHermitianTridiagCtrlDefault_s( ElHermitianTridiagCtrl* ctrl )
{
    ctrl->approach = EL_HERMITIAN_TRIDIAG_DEFAULT;
    ctrl->order = EL_ROW_MAJOR;
    ElSymvCtrlDefault_s( &ctrl->symvCtrl );
    return EL_SUCCESS;
}
ElError ElHermitianTridiagCtrlDefault_d( ElHermitianTridiagCtrl* ctrl )
{
    ctrl->approach = EL_HERMITIAN_TRIDIAG_DEFAULT;
    ctrl->order = EL_ROW_MAJOR;
    ElSymvCtrlDefault_d( &ctrl->symvCtrl );
    return EL_SUCCESS;
}
ElError ElHermitianTridiagCtrlDefault_c( ElHermitianTridiagCtrl* ctrl )
{
    ctrl->approach = EL_HERMITIAN_TRIDIAG_DEFAULT;
    ctrl->order = EL_ROW_MAJOR;
    ElSymvCtrlDefault_c( &ctrl->symvCtrl );
    return EL_SUCCESS;
}
ElError ElHermitianTridiagCtrlDefault_z( ElHermitianTridiagCtrl* ctrl )
{
    ctrl->approach = EL_HERMITIAN_TRIDIAG_DEFAULT;
    ctrl->order = EL_ROW_MAJOR;
    ElSymvCtrlDefault_z( &ctrl->symvCtrl );
    return EL_SUCCESS;
}

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Bidiag
     ====== */ \
  /* Return the packed reduction to bidiagonal form, B := Q^H A P */ \
  ElError ElBidiag_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG tP, ElMatrix_ ## SIG tQ ) \
  { EL_TRY( Bidiag( *CReflect(A), *CReflect(tP), *CReflect(tQ) ) ) } \
  ElError ElBidiagDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG tP, ElDistMatrix_ ## SIG tQ ) \
  { EL_TRY( Bidiag( *CReflect(A), *CReflect(tP), *CReflect(tQ) ) ) } \
  /* Overwrite A with B = Q^H A P and also return P and Q */ \
  ElError ElBidiagExplicit_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG P, ElMatrix_ ## SIG Q ) \
  { EL_TRY( bidiag::Explicit( *CReflect(A), *CReflect(P), *CReflect(Q) ) ) } \
  ElError ElBidiagExplicitDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG P, ElDistMatrix_ ## SIG Q ) \
  { EL_TRY( bidiag::Explicit( *CReflect(A), *CReflect(P), *CReflect(Q) ) ) } \
  /* Only return the condensed bidiagonal matrix, B := Q^H A P */ \
  ElError ElBidiagOnly_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( bidiag::ExplicitCondensed( *CReflect(A) ) ) } \
  ElError ElBidiagOnlyDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( bidiag::ExplicitCondensed( *CReflect(A) ) ) } \
  /* Apply Q from B := Q^H A P to a set of vectors */ \
  ElError ElApplyQAfterBidiag_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, ElMatrix_ ## SIG B ) \
  { EL_TRY( bidiag::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(B) ) ) } \
  ElError ElApplyQAfterBidiagDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( bidiag::ApplyQ( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(B) ) ) } \
  /* Apply P from B := Q^H A P to a set of vectors */ \
  ElError ElApplyPAfterBidiag_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, ElMatrix_ ## SIG B ) \
  { EL_TRY( bidiag::ApplyP( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(B) ) ) } \
  ElError ElApplyPAfterBidiagDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( bidiag::ApplyP( \
      CReflect(side), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(B) ) ) } \
  /* HermitianTridiag
     ================ */ \
  /* Return the packed reduction to condensed form */ \
  ElError ElHermitianTridiag_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIG t ) \
  { EL_TRY( HermitianTridiag( \
      CReflect(uplo), *CReflect(A), *CReflect(t) ) ) } \
  ElError ElHermitianTridiagDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t ) \
  { EL_TRY( HermitianTridiag( \
      CReflect(uplo), *CReflect(A), *CReflect(t) ) ) } \
  /* Return only the condensed form */ \
  ElError ElHermitianTridiagOnly_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( herm_tridiag::ExplicitCondensed \
      ( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianTridiagOnlyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( herm_tridiag::ExplicitCondensed \
      ( CReflect(uplo), *CReflect(A) ) ) } \
  /* ApplyQ after HermitianTridiag */ \
  ElError ElApplyQAfterHermitianTridiag_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, ElMatrix_ ## SIG B ) \
  { EL_TRY( herm_tridiag::ApplyQ( \
    CReflect(side), CReflect(uplo), CReflect(orientation), \
    *CReflect(A), *CReflect(t), *CReflect(B) ) ) } \
  ElError ElApplyQAfterHermitianTridiagDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( herm_tridiag::ApplyQ( \
    CReflect(side), CReflect(uplo), CReflect(orientation), \
    *CReflect(A), *CReflect(t), *CReflect(B) ) ) } \
  /* Hessenberg
     ========== */ \
  /* Packed reduction to Hessenberg form, H := Q^H A Q */ \
  ElError ElHessenberg_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIG t ) \
  { EL_TRY( Hessenberg( \
      CReflect(uplo), *CReflect(A), *CReflect(t) ) ) } \
  ElError ElHessenbergDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t ) \
  { EL_TRY( Hessenberg( \
      CReflect(uplo), *CReflect(A), *CReflect(t) ) ) } \
  /* Only return the similar Hessenberg matrix */ \
  ElError ElHessenbergOnly_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( hessenberg::ExplicitCondensed( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHessenbergOnlyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( hessenberg::ExplicitCondensed( CReflect(uplo), *CReflect(A) ) ) } \
  /* Apply Q from a Hessenberg decomposition, H := Q^H A Q */ \
  ElError ElApplyQAfterHessenberg_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, \
    ElMatrix_ ## SIG B ) \
  { EL_TRY( hessenberg::ApplyQ( \
      CReflect(side), CReflect(uplo), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(B) ) ) } \
  ElError ElApplyQAfterHessenbergDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( hessenberg::ApplyQ( \
      CReflect(side), CReflect(uplo), CReflect(orientation), \
      *CReflect(A), *CReflect(t), *CReflect(B) ) ) }

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F)

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
