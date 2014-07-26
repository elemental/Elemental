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

#define DM_VR_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VR,STAR>&>(*CReflect(A))
#define DM_VR_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VR,STAR>&>(*CReflect(A))

extern "C" {

#define C_PROTO_BASE(SIG,SIGBASE,T) \
  /* Gemm */ \
  ElError ElGemm_ ## SIG \
  ( ElOrientation orientationOfA, ElOrientation orientationOfB, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Gemm( CReflect(orientationOfA), CReflect(orientationOfB), \
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C)) ) } \
  ElError ElGemmDist_ ## SIG \
  ( ElOrientation orientationOfA, ElOrientation orientationOfB, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Gemm( CReflect(orientationOfA), CReflect(orientationOfB), \
            CReflect(alpha), DM_CAST_CONST(T,A), DM_CAST_CONST(T,B), \
            CReflect(beta), DM_CAST(T,C) ) ) } \
  ElError ElGemmXDist_ ## SIG \
  ( ElOrientation orientationOfA, ElOrientation orientationOfB, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElDistMatrix_ ## SIG C, ElGemmAlgorithm alg ) \
  { EL_TRY( \
      Gemm( CReflect(orientationOfA), CReflect(orientationOfB), \
            CReflect(alpha), DM_CAST_CONST(T,A), DM_CAST_CONST(T,B), \
            CReflect(beta), DM_CAST(T,C), CReflect(alg) ) ) } 

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* MultiShiftQuasiTrsm */ \
  ElError ElMultiShiftQuasiTrsm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG shifts, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      MultiShiftQuasiTrsm( \
        CReflect(side), CReflect(uplo), CReflect(orientation), \
        CReflect(alpha),*CReflect(A), *CReflect(shifts), \
        *CReflect(B) ) ) } \
  ElError ElMultiShiftQuasiTrsmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG shifts, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      MultiShiftQuasiTrsm( \
        CReflect(side), CReflect(uplo), CReflect(orientation), \
        CReflect(alpha), DM_CAST_CONST(F,A), \
        DM_VR_STAR_CAST_CONST(F,shifts), DM_CAST(F,B) ) ) } \
  /* MultiShiftTrsm */ \
  ElError ElMultiShiftTrsm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG shifts, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      MultiShiftTrsm( \
        CReflect(side), CReflect(uplo), CReflect(orientation), \
        CReflect(alpha), *CReflect(A), *CReflect(shifts), \
        *CReflect(B) ) ) } \
  ElError ElMultiShiftTrsmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG shifts, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      MultiShiftTrsm( \
        CReflect(side), CReflect(uplo), CReflect(orientation), \
        CReflect(alpha), DM_CAST_CONST(F,A), \
        DM_VR_STAR_CAST_CONST(F,shifts), DM_CAST(F,B) ) ) } \
  /* QuasiTrsm */ \
  ElError ElQuasiTrsm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      QuasiTrsm( \
        CReflect(side), CReflect(uplo), CReflect(orientation), \
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
  ElError ElQuasiTrsmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      QuasiTrsm( \
        CReflect(side), CReflect(uplo), CReflect(orientation), \
        CReflect(alpha), DM_CAST_CONST(F,A), DM_CAST(F,B) ) ) } \
  /* Symm */ \
  ElError ElSymm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    CREFLECT(F) beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Symm( CReflect(side), CReflect(uplo), \
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElSymmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(F) beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Symm( CReflect(side), CReflect(uplo), \
            CReflect(alpha), DM_CAST_CONST(F,A), DM_CAST_CONST(F,B), \
            CReflect(beta), DM_CAST(F,C) ) ) } \
  /* Syrk */ \
  ElError ElSyrk_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, \
    CREFLECT(F) beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Syrk( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), *CReflect(A), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElSyrkDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
    CREFLECT(F) beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Syrk( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), DM_CAST_CONST(F,A), \
            CReflect(beta), DM_CAST(F,C) ) ) } \
  /* Syr2k */ \
  ElError ElSyr2k_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    CREFLECT(F) beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Syr2k( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElSyr2kDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(F) beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Syr2k( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), DM_CAST_CONST(F,A), DM_CAST_CONST(F,B), \
            CReflect(beta), DM_CAST(F,C) ) ) } \
  /* Trdtrmm */ \
  ElError ElTrdtrmm_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trdtrmm( CReflect(uplo), *CReflect(A), conjugate ) ) } \
  ElError ElTrdtrmmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trdtrmm( CReflect(uplo), DM_CAST(F,A), conjugate ) ) } \
  /* TrdtrmmQuasi */ \
  ElError ElTrdtrmmQuasi_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG dOff, bool conjugate ) \
  { EL_TRY( \
      Trdtrmm( \
        CReflect(uplo), *CReflect(A), *CReflect(dOff), \
        conjugate ) ) } \
  ElError ElTrdtrmmQuasiDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG dOff, bool conjugate ) \
  { EL_TRY( \
      Trdtrmm( \
        CReflect(uplo), DM_CAST(F,A), DM_MD_STAR_CAST_CONST(F,dOff), \
        conjugate ) ) } \
  /* Trmm */ \
  ElError ElTrmm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElOrientation orientation, ElUnitOrNonUnit diag, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      Trmm( \
        CReflect(side), CReflect(uplo), \
        CReflect(orientation), CReflect(diag), \
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
  ElError ElTrmmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElOrientation orientation, ElUnitOrNonUnit diag, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      Trmm( \
        CReflect(side), CReflect(uplo), \
        CReflect(orientation), CReflect(diag), \
        CReflect(alpha), DM_CAST_CONST(F,A), DM_CAST(F,B) ) ) } \
  /* Trsm */ \
  ElError ElTrsm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElOrientation orientation, ElUnitOrNonUnit diag, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      Trsm( \
        CReflect(side), CReflect(uplo), \
        CReflect(orientation), CReflect(diag), \
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
  ElError ElTrsmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElOrientation orientation, ElUnitOrNonUnit diag, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      Trsm( \
        CReflect(side), CReflect(uplo), \
        CReflect(orientation), CReflect(diag), \
        CReflect(alpha), DM_CAST_CONST(F,A), DM_CAST(F,B) ) ) } \
  /* Trstrm */ \
  ElError ElTrstrm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElOrientation orientation, ElUnitOrNonUnit diag, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( \
      Trstrm( \
        CReflect(side), CReflect(uplo), \
        CReflect(orientation), CReflect(diag), \
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
  ElError ElTrstrmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElOrientation orientation, ElUnitOrNonUnit diag, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      Trstrm( \
        CReflect(side), CReflect(uplo), \
        CReflect(orientation), CReflect(diag), \
        CReflect(alpha), DM_CAST_CONST(F,A), DM_CAST(F,B) ) ) } \
  /* Trtrmm */ \
  ElError ElTrtrmm_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trtrmm( CReflect(uplo), *CReflect(A), conjugate ) ) } \
  ElError ElTrtrmmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trtrmm( CReflect(uplo), DM_CAST(F,A), conjugate ) ) } \
  /* TwoSidedTrmm */ \
  ElError ElTwoSidedTrmm_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, \
    ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG B ) \
  { EL_TRY( \
      TwoSidedTrmm( \
        CReflect(uplo), CReflect(diag), \
        *CReflect(A), *CReflect(B) ) ) } \
  ElError ElTwoSidedTrmmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, \
    ElDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      TwoSidedTrmm( \
        CReflect(uplo), CReflect(diag), \
        DM_CAST(F,A), DM_CAST_CONST(F,B) ) ) } \
  /* TwoSidedTrsm */ \
  ElError ElTwoSidedTrsm_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, \
    ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG B ) \
  { EL_TRY( \
      TwoSidedTrsm( \
        CReflect(uplo), CReflect(diag), \
        *CReflect(A), *CReflect(B) ) ) } \
  ElError ElTwoSidedTrsmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElUnitOrNonUnit diag, \
    ElDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B ) \
  { EL_TRY( \
      TwoSidedTrsm( \
        CReflect(uplo), CReflect(diag), \
        DM_CAST(F,A), DM_CAST_CONST(F,B) ) ) }

#define C_PROTO_INT(SIG,T) C_PROTO_BASE(SIG,SIG,T)

#define C_PROTO_REAL(SIG,T) \
  C_PROTO_BASE(SIG,SIG,T) \
  C_PROTO_FIELD(SIG,SIG,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO_BASE(SIG,SIGBASE,T) \
  C_PROTO_FIELD(SIG,SIGBASE,T) \
  /* Hemm */ \
  ElError ElHemm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Hemm( CReflect(side), CReflect(uplo), \
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElHemmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Hemm( CReflect(side), CReflect(uplo), \
            CReflect(alpha), DM_CAST_CONST(T,A), DM_CAST_CONST(T,B), \
            CReflect(beta), DM_CAST(T,C) ) ) } \
  /* Herk */ \
  ElError ElHerk_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, \
    CREFLECT(T) beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Herk( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), *CReflect(A), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElHerkDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
    CREFLECT(T) beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Herk( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), DM_CAST_CONST(T,A), \
            CReflect(beta), DM_CAST(T,C) ) ) } \
  /* Her2k */ \
  ElError ElHer2k_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Her2k( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElHer2kDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Her2k( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), DM_CAST_CONST(T,A), DM_CAST_CONST(T,B), \
            CReflect(beta), DM_CAST(T,C) ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
