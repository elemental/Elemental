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

#define DM_CAST(T,A) dynamic_cast<DistMatrix<T>&>(*CReflect(A))
#define DM_CAST_CONST(T,A) dynamic_cast<const DistMatrix<T>&>(*CReflect(A))

#define DM_MD_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,MD,STAR>&>(*CReflect(A))
#define DM_MD_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,MD,STAR>&>(*CReflect(A))

#define DM_STAR_VR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,STAR,VR>&>(*CReflect(A))
#define DM_STAR_VR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,STAR,VR>&>(*CReflect(A))

#define DM_STAR_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,STAR,STAR>&>(*CReflect(A))
#define DM_STAR_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,STAR,STAR>&>(*CReflect(A))

#define DM_VC_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VC,STAR>&>(*CReflect(A))
#define DM_VC_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VC,STAR>&>(*CReflect(A))

#define DM_VR_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VR,STAR>&>(*CReflect(A))
#define DM_VR_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VR,STAR>&>(*CReflect(A))

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
  { EL_TRY( Bidiag( \
      DM_CAST(F,A), DM_STAR_STAR_CAST(F,tP), DM_STAR_STAR_CAST(F,tQ) ) ) } \
  /* Only return the condensed bidiagonal matrix, B := Q^H A P */ \
  ElError ElBidiagOnly_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Bidiag( *CReflect(A) ) ) } \
  ElError ElBidiagOnlyDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Bidiag( DM_CAST(F,A) ) ) } \
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
  { EL_TRY( \
      if( CReflect(t)->DistData().colDist == MD ) \
        bidiag::ApplyQ \
        ( CReflect(side), CReflect(orientation), \
          DM_CAST_CONST(F,A), DM_MD_STAR_CAST_CONST(F,t), DM_CAST(F,B) ); \
      else \
        bidiag::ApplyQ \
        ( CReflect(side), CReflect(orientation), \
          DM_CAST_CONST(F,A), DM_STAR_STAR_CAST_CONST(F,t), DM_CAST(F,B) ) ) } \
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
  { EL_TRY( \
      if( CReflect(t)->DistData().colDist == MD ) \
        bidiag::ApplyP \
        ( CReflect(side), CReflect(orientation), \
          DM_CAST_CONST(F,A), DM_MD_STAR_CAST_CONST(F,t), DM_CAST(F,B) ); \
      else \
        bidiag::ApplyP \
        ( CReflect(side), CReflect(orientation), \
          DM_CAST_CONST(F,A), DM_STAR_STAR_CAST_CONST(F,t), DM_CAST(F,B) ) ) } \
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
      CReflect(uplo), DM_CAST(F,A), DM_STAR_STAR_CAST(F,t) ) ) } \
  ElError ElHermitianTridiagXDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t, \
    ElHermitianTridiagCtrl ctrl ) \
  { EL_TRY( HermitianTridiag( \
      CReflect(uplo), DM_CAST(F,A), DM_STAR_STAR_CAST(F,t), \
      CReflect(ctrl) ) ) } \
  /* Return only the condensed form */ \
  ElError ElHermitianTridiagOnly_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianTridiag( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianTridiagOnlyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianTridiag( CReflect(uplo), DM_CAST(F,A) ) ) } \
  ElError ElHermitianTridiagOnlyXDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElHermitianTridiagCtrl ctrl ) \
  { EL_TRY( HermitianTridiag( \
      CReflect(uplo), DM_CAST(F,A), CReflect(ctrl) ) ) } \
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
  { EL_TRY( \
    if( CReflect(t)->DistData().colDist == MD ) \
      herm_tridiag::ApplyQ \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        DM_CAST_CONST(F,A), DM_MD_STAR_CAST_CONST(F,t), DM_CAST(F,B) ); \
    else \
      herm_tridiag::ApplyQ \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        DM_CAST_CONST(F,A), DM_STAR_STAR_CAST_CONST(F,t), DM_CAST(F,B) ) ) } \
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
      CReflect(uplo), DM_CAST(F,A), DM_STAR_STAR_CAST(F,t) ) ) } \
  /* Only return the similar Hessenberg matrix */ \
  ElError ElHessenbergOnly_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( Hessenberg( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHessenbergOnlyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Hessenberg( CReflect(uplo), DM_CAST(F,A) ) ) } \
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
  { EL_TRY( \
    if( CReflect(t)->DistData().colDist == MD ) \
      hessenberg::ApplyQ \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        DM_CAST_CONST(F,A), DM_MD_STAR_CAST_CONST(F,t), DM_CAST(F,B) ); \
    else \
      hessenberg::ApplyQ \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        DM_CAST_CONST(F,A), DM_STAR_STAR_CAST_CONST(F,t), DM_CAST(F,B) ) ) }

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F)

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
