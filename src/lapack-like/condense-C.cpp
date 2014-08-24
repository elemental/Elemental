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

extern "C" {

ElError ElHermitianTridiagCtrlDefault( ElHermitianTridiagCtrl* ctrl )
{
    ctrl->approach = EL_HERMITIAN_TRIDIAG_DEFAULT;
    ctrl->order = EL_ROW_MAJOR;
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
  /* Only return the condensed bidiagonal matrix, B := Q^H A P */ \
  ElError ElBidiagOnly_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Bidiag( *CReflect(A) ) ) } \
  ElError ElBidiagOnlyDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Bidiag( *CReflect(A) ) ) } \
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
  ElError ElHermitianTridiagXDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t, \
    ElHermitianTridiagCtrl ctrl ) \
  { EL_TRY( HermitianTridiag( \
      CReflect(uplo), *CReflect(A), *CReflect(t), CReflect(ctrl) ) ) } \
  /* Return only the condensed form */ \
  ElError ElHermitianTridiagOnly_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianTridiag( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianTridiagOnlyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianTridiag( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianTridiagOnlyXDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElHermitianTridiagCtrl ctrl ) \
  { EL_TRY( HermitianTridiag( CReflect(uplo), *CReflect(A), \
      CReflect(ctrl) ) ) } \
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
  { EL_TRY( Hessenberg( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHessenbergOnlyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Hessenberg( CReflect(uplo), *CReflect(A) ) ) } \
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
