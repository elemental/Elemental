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
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElGemmXDist_ ## SIG \
  ( ElOrientation orientationOfA, ElOrientation orientationOfB, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T) beta, ElDistMatrix_ ## SIG C, ElGemmAlgorithm alg ) \
  { EL_TRY( \
      Gemm( CReflect(orientationOfA), CReflect(orientationOfB), \
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C), CReflect(alg) ) ) } \
  ElError ElMultiply_ ## SIG \
  ( ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG X, CREFLECT(T) beta, ElMatrix_ ## SIG Y ) \
  { EL_TRY( Multiply( CReflect(orientation), \
                      CReflect(alpha), *CReflect(A), *CReflect(X), \
                      CReflect(beta), *CReflect(Y) ) ) } \
  ElError ElMultiplyDist_ ## SIG \
  ( ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG X, CREFLECT(T) beta, \
    ElDistMultiVec_ ## SIG Y ) \
  { EL_TRY( Multiply( CReflect(orientation), \
                      CReflect(alpha), *CReflect(A), *CReflect(X), \
                      CReflect(beta), *CReflect(Y) ) ) }

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
        CReflect(alpha), *CReflect(A), \
        *CReflect(shifts), *CReflect(B) ) ) } \
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
        CReflect(alpha), *CReflect(A), \
        *CReflect(shifts), *CReflect(B) ) ) } \
  /* SafeMultiShiftTrsm */ \
  ElError ElSafeMultiShiftTrsm_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG shifts, ElMatrix_ ## SIG B, \
    ElMatrix_ ## SIG scales ) \
  { EL_TRY( \
      SafeMultiShiftTrsm( \
	CReflect(side), CReflect(uplo), CReflect(orientation), \
	CReflect(alpha), *CReflect(A), \
	*CReflect(shifts), *CReflect(B), \
	*CReflect(scales) ) ) } \
  ElError ElSafeMultiShiftTrsmDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG shifts, ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG scales ) \
  { EL_TRY( \
      SafeMultiShiftTrsm( \
	CReflect(side), CReflect(uplo), CReflect(orientation), \
	CReflect(alpha), *CReflect(A), \
	*CReflect(shifts), *CReflect(B), \
	*CReflect(scales) ) ) } \
  /* QuasiTrsm */	      \
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
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
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
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
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
            CReflect(alpha), *CReflect(A), \
            CReflect(beta), *CReflect(C) ) ) } \
  ElError ElSyrkSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstSparseMatrix_ ## SIG A, \
    CREFLECT(F) beta,  ElSparseMatrix_ ## SIG C ) \
  { EL_TRY( \
      Syrk( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), *CReflect(A), CReflect(beta), *CReflect(C) ) ) } \
  ElError ElSyrkDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(F) alpha, ElConstDistSparseMatrix_ ## SIG A, \
    CREFLECT(F) beta,  ElDistSparseMatrix_ ## SIG C ) \
  { EL_TRY( \
      Syrk( CReflect(uplo), CReflect(orientation), \
            CReflect(alpha), *CReflect(A), CReflect(beta), *CReflect(C) ) ) } \
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
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
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
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
  /* Trrk */ \
  ElError ElTrrk_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    CREFLECT(F) beta,  ElMatrix_ ## SIG C ) \
  { EL_TRY( Trrk( CReflect(uplo), CReflect(orientA), CReflect(orientB), \
      CReflect(alpha), *CReflect(A), *CReflect(B), \
      CReflect(beta), *CReflect(C) ) ) } \
  ElError ElTrrkDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(F) beta,  ElDistMatrix_ ## SIG C ) \
  { EL_TRY( Trrk( CReflect(uplo), CReflect(orientA), CReflect(orientB), \
      CReflect(alpha), *CReflect(A), *CReflect(B), \
      CReflect(beta), *CReflect(C) ) ) } \
  /* Trr2k */ \
  /*
  ElError ElTrr2k_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElOrientation orientA, ElOrientation orientB, \
    ElOrientation orientC, ElOrientation orientD, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    CREFLECT(F) beta,  ElConstMatrix_ ## SIG C, ElConstMatrix_ ## SIG D, \
    CREFLECT(F) gamma, ElMatrix_ ## SIG E ) \
  { EL_TRY( Trr2k( CReflect(uplo), \
      CReflect(orientA), CReflect(orientB), \
      CReflect(orientC), CReflect(orientD), \
      CReflect(alpha), *CReflect(A), *CReflect(B), \
      CReflect(beta),  *CReflect(C), *CReflect(D), \
      CReflect(gamma), *CReflect(E) ) ) } \
  */ \
  ElError ElTrr2kDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElOrientation orientA, ElOrientation orientB, \
    ElOrientation orientC, ElOrientation orientD, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    CREFLECT(F) beta,  ElConstDistMatrix_ ## SIG C, \
                       ElConstDistMatrix_ ## SIG D, \
    CREFLECT(F) gamma, ElDistMatrix_ ## SIG E ) \
  { EL_TRY( Trr2k( CReflect(uplo), \
      CReflect(orientA), CReflect(orientB), \
      CReflect(orientC), CReflect(orientD), \
      CReflect(alpha), *CReflect(A), *CReflect(B), \
      CReflect(beta),  *CReflect(C), *CReflect(D), \
      CReflect(gamma), *CReflect(E) ) ) } \
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
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
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
        CReflect(alpha), *CReflect(A), *CReflect(B) ) ) } \
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
        *CReflect(A), *CReflect(B) ) ) } \
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
        *CReflect(A), *CReflect(B) ) ) }

#define C_PROTO_INT(SIG,T) C_PROTO_BASE(SIG,SIG,T)

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_BASE(SIG,SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* Trdtrmm */ \
  ElError ElTrdtrmm_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( Trdtrmm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElTrdtrmmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Trdtrmm( CReflect(uplo), *CReflect(A) ) ) } \
  /* TrdtrmmQuasi */ \
  ElError ElTrdtrmmQuasi_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG dOff ) \
  { EL_TRY( \
      Trdtrmm( CReflect(uplo), *CReflect(A), *CReflect(dOff) ) ) } \
  ElError ElTrdtrmmQuasiDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG dOff ) \
  { EL_TRY( \
      Trdtrmm( CReflect(uplo), *CReflect(A), *CReflect(dOff) ) ) } \
  /* Trtrmm */ \
  ElError ElTrtrmm_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( Trtrmm( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElTrtrmmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Trtrmm( CReflect(uplo), *CReflect(A) ) ) }

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
            CReflect(alpha), *CReflect(A), *CReflect(B), \
            CReflect(beta), *CReflect(C) ) ) } \
  /* Herk */ \
  ElError ElHerk_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    Base<T> alpha, ElConstMatrix_ ## SIG A, \
    Base<T> beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Herk( CReflect(uplo), CReflect(orientation), \
            alpha, *CReflect(A), beta, *CReflect(C) ) ) } \
  ElError ElHerkDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    Base<T> alpha, ElConstDistMatrix_ ## SIG A, \
    Base<T> beta, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Herk( CReflect(uplo), CReflect(orientation), \
            alpha, *CReflect(A), beta, *CReflect(C) ) ) } \
  ElError ElHerkSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    Base<T> alpha, ElConstSparseMatrix_ ## SIG A, \
    Base<T> beta,  ElSparseMatrix_ ## SIG C ) \
  { EL_TRY( \
      Herk( CReflect(uplo), CReflect(orientation), \
            alpha, *CReflect(A), beta, *CReflect(C) ) ) } \
  ElError ElHerkDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    Base<T> alpha, ElConstDistSparseMatrix_ ## SIG A, \
    Base<T> beta,  ElDistSparseMatrix_ ## SIG C ) \
  { EL_TRY( \
      Herk( CReflect(uplo), CReflect(orientation), \
            alpha, *CReflect(A), beta, *CReflect(C) ) ) } \
  /* Her2k */ \
  ElError ElHer2k_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    Base<T> beta, ElMatrix_ ## SIG C ) \
  { EL_TRY( \
      Her2k( CReflect(uplo), CReflect(orientation), \
             CReflect(alpha), *CReflect(A), *CReflect(B), \
             beta, *CReflect(C) ) ) } \
  ElError ElHer2kDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG B, \
    Base<T> beta,  ElDistMatrix_ ## SIG C ) \
  { EL_TRY( \
      Her2k( CReflect(uplo), CReflect(orientation), \
             CReflect(alpha), *CReflect(A), *CReflect(B), \
             beta, *CReflect(C) ) ) } \
  /* Trdtrmm */ \
  ElError ElTrdtrmm_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trdtrmm( CReflect(uplo), *CReflect(A), conjugate ) ) } \
  ElError ElTrdtrmmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trdtrmm( CReflect(uplo), *CReflect(A), conjugate ) ) } \
  /* TrdtrmmQuasi */ \
  ElError ElTrdtrmmQuasi_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG dOff, bool conjugate ) \
  { EL_TRY( \
      Trdtrmm( CReflect(uplo), *CReflect(A), *CReflect(dOff), conjugate ) ) } \
  ElError ElTrdtrmmQuasiDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG dOff, bool conjugate ) \
  { EL_TRY( \
      Trdtrmm( CReflect(uplo), *CReflect(A), *CReflect(dOff), conjugate ) ) } \
  /* Trtrmm */ \
  ElError ElTrtrmm_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trtrmm( CReflect(uplo), *CReflect(A), conjugate ) ) } \
  ElError ElTrtrmmDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( Trtrmm( CReflect(uplo), *CReflect(A), conjugate ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
