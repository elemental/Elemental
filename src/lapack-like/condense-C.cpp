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

#define DM_MD_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,MD,STAR>&>(*Reinterpret(A))
#define DM_MD_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,MD,STAR>&>(*Reinterpret(A))

#define DM_STAR_VR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,STAR,VR>&>(*Reinterpret(A))
#define DM_STAR_VR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,STAR,VR>&>(*Reinterpret(A))

#define DM_STAR_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,STAR,STAR>&>(*Reinterpret(A))
#define DM_STAR_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,STAR,STAR>&>(*Reinterpret(A))

#define DM_VC_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VC,STAR>&>(*Reinterpret(A))
#define DM_VC_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VC,STAR>&>(*Reinterpret(A))

#define DM_VR_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VR,STAR>&>(*Reinterpret(A))
#define DM_VR_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VR,STAR>&>(*Reinterpret(A))

extern "C" {

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Bidiag
     ====== */ \
  /* Return the packed reduction to bidiagonal form, B := Q^H A P */ \
  ElError ElBidiag_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG tP, ElMatrix_ ## SIG tQ ) \
  { EL_TRY( Bidiag( *Reinterpret(A), *Reinterpret(tP), *Reinterpret(tQ) ) ) } \
  ElError ElBidiagDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG tP, ElDistMatrix_ ## SIG tQ ) \
  { EL_TRY( Bidiag( \
      DM_CAST(F,A), DM_STAR_STAR_CAST(F,tP), DM_STAR_STAR_CAST(F,tQ) ) ) } \
  /* Only return the condensed bidiagonal matrix, B := Q^H A P */ \
  ElError ElBidiagOnly_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Bidiag( *Reinterpret(A) ) ) } \
  ElError ElBidiagOnlyDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Bidiag( DM_CAST(F,A) ) ) } \
  /* Apply Q from B := Q^H A P to a set of vectors */ \
  ElError ElApplyQAfterBidiag_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, ElMatrix_ ## SIG B ) \
  { EL_TRY( bidiag::ApplyQ( \
      Reinterpret(side), Reinterpret(orientation), \
      *Reinterpret(A), *Reinterpret(t), *Reinterpret(B) ) ) } \
  ElError ElApplyQAfterBidiagDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      if( Reinterpret(t)->DistData().colDist == MD ) \
        bidiag::ApplyQ \
        ( Reinterpret(side), Reinterpret(orientation), \
          DM_CAST_CONST(F,A), DM_MD_STAR_CAST_CONST(F,t), DM_CAST(F,B) ); \
      else \
        bidiag::ApplyQ \
        ( Reinterpret(side), Reinterpret(orientation), \
          DM_CAST_CONST(F,A), DM_STAR_STAR_CAST_CONST(F,t), DM_CAST(F,B) ) ) } \
  /* Apply P from B := Q^H A P to a set of vectors */ \
  ElError ElApplyPAfterBidiag_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, ElMatrix_ ## SIG B ) \
  { EL_TRY( bidiag::ApplyP( \
      Reinterpret(side), Reinterpret(orientation), \
      *Reinterpret(A), *Reinterpret(t), *Reinterpret(B) ) ) } \
  ElError ElApplyPAfterBidiagDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      if( Reinterpret(t)->DistData().colDist == MD ) \
        bidiag::ApplyP \
        ( Reinterpret(side), Reinterpret(orientation), \
          DM_CAST_CONST(F,A), DM_MD_STAR_CAST_CONST(F,t), DM_CAST(F,B) ); \
      else \
        bidiag::ApplyP \
        ( Reinterpret(side), Reinterpret(orientation), \
          DM_CAST_CONST(F,A), DM_STAR_STAR_CAST_CONST(F,t), DM_CAST(F,B) ) ) } \
  /* HermitianTridiag
     ================ */ \
  /* Return the packed reduction to condensed form */ \
  ElError ElHermitianTridiag_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIG t ) \
  { EL_TRY( HermitianTridiag( \
      Reinterpret(uplo), *Reinterpret(A), *Reinterpret(t) ) ) } \
  ElError ElHermitianTridiagDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t ) \
  { EL_TRY( HermitianTridiag( \
      Reinterpret(uplo), DM_CAST(F,A), DM_STAR_STAR_CAST(F,t) ) ) } \
  ElError ElHermitianTridiagXDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG t, \
    ElHermitianTridiagCtrl ctrl ) \
  { EL_TRY( HermitianTridiag( \
      Reinterpret(uplo), DM_CAST(F,A), DM_STAR_STAR_CAST(F,t), \
      Reinterpret(ctrl) ) ) } \
  /* Return only the condensed form */ \
  ElError ElHermitianTridiagOnly_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianTridiag( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElHermitianTridiagOnlyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianTridiag( Reinterpret(uplo), DM_CAST(F,A) ) ) } \
  ElError ElHermitianTridiagOnlyXDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElHermitianTridiagCtrl ctrl ) \
  { EL_TRY( HermitianTridiag( \
      Reinterpret(uplo), DM_CAST(F,A), Reinterpret(ctrl) ) ) } \
  /* ApplyQ after HermitianTridiag */ \
  ElError ElApplyQAfterHermitianTridiag_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG t, ElMatrix_ ## SIG B ) \
  { EL_TRY( herm_tridiag::ApplyQ( \
    Reinterpret(side), Reinterpret(uplo), Reinterpret(orientation), \
    *Reinterpret(A), *Reinterpret(t), *Reinterpret(B) ) ) } \
  ElError ElApplyQAfterHermitianTridiagDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG t, \
    ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
    if( Reinterpret(t)->DistData().colDist == MD ) \
      herm_tridiag::ApplyQ \
      ( Reinterpret(side), Reinterpret(uplo), Reinterpret(orientation), \
        DM_CAST_CONST(F,A), DM_MD_STAR_CAST_CONST(F,t), DM_CAST(F,B) ); \
    else \
      herm_tridiag::ApplyQ \
      ( Reinterpret(side), Reinterpret(uplo), Reinterpret(orientation), \
        DM_CAST_CONST(F,A), DM_STAR_STAR_CAST_CONST(F,t), DM_CAST(F,B) ) ) }

#define C_PROTO_REAL(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* TODO */

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* TODO */

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
