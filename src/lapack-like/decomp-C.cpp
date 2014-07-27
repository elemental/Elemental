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

/* HermitianSdcCtrl */
ElError ElHermitianSdcCtrlDefault_s( ElHermitianSdcCtrl_s* ctrl )
{
    ctrl->cutoff = 256;
    ctrl->maxInnerIts = 2;
    ctrl->maxOuterIts = 10;
    ctrl->tol = 0;
    ctrl->spreadFactor = 1e-6f;
    ctrl->progress = false;
    return EL_SUCCESS;
}
ElError ElHermitianSdcCtrlDefault_d( ElHermitianSdcCtrl_d* ctrl )
{
    ctrl->cutoff = 256;
    ctrl->maxInnerIts = 2;
    ctrl->maxOuterIts = 10;
    ctrl->tol = 0;
    ctrl->spreadFactor = 1e-6;
    ctrl->progress = false;
    return EL_SUCCESS;
}

/* HermitianEigSubset */
ElError ElHermitianEigSubsetDefault_s( ElHermitianEigSubset_s* subset )
{
    subset->indexSubset = false;
    subset->lowerIndex = 0;
    subset->upperIndex = 0;
    subset->rangeSubset = false;
    subset->lowerBound = 0;
    subset->upperBound = 0;
    return EL_SUCCESS;
}
ElError ElHermitianEigSubsetDefault_d( ElHermitianEigSubset_d* subset )
{
    subset->indexSubset = false;
    subset->lowerIndex = 0;
    subset->upperIndex = 0;
    subset->rangeSubset = false;
    subset->lowerBound = 0;
    subset->upperBound = 0;
    return EL_SUCCESS;
}

/* HermitianEigCtrl */
ElError ElHermitianEigCtrlDefault_s( ElHermitianEigCtrl_s* ctrl )
{
    ElHermitianTridiagCtrlDefault( &ctrl->tridiagCtrl );
    ElHermitianSdcCtrlDefault_s( &ctrl->sdcCtrl );
    ctrl->useSdc = false;
    return EL_SUCCESS;
}
ElError ElHermitianEigCtrlDefault_d( ElHermitianEigCtrl_d* ctrl )
{
    ElHermitianTridiagCtrlDefault( &ctrl->tridiagCtrl );
    ElHermitianSdcCtrlDefault_d( &ctrl->sdcCtrl );
    ctrl->useSdc = false;
    return EL_SUCCESS;
}

/* PolarCtrl */
ElError ElPolarCtrlDefault( ElPolarCtrl* ctrl )
{
    ctrl->qdwh = false;
    ctrl->colPiv = false;
    ctrl->maxIts = 20;
    ctrl->numIts = 0;
    return EL_SUCCESS;
}

/* HessQrCtrl */
ElError ElHessQrCtrlDefault( ElHessQrCtrl* ctrl )
{
    ctrl->aed = false;
    ctrl->blockHeight = DefaultBlockHeight();
    ctrl->blockWidth = DefaultBlockWidth();
    return EL_SUCCESS;
}

/* SdcCtrl */
ElError ElSdcCtrlDefault_s( ElSdcCtrl_s* ctrl )
{
    ctrl->cutoff = 256;
    ctrl->maxInnerIts = 2;
    ctrl->maxOuterIts = 10;
    ctrl->tol = 0;
    ctrl->spreadFactor = 1e-6f;
    ctrl->random = true;
    ctrl->progress = false;
    ElSignCtrlDefault_s( &ctrl->signCtrl );
    return EL_SUCCESS;
}
ElError ElSdcCtrlDefault_d( ElSdcCtrl_d* ctrl )
{
    ctrl->cutoff = 256;
    ctrl->maxInnerIts = 2;
    ctrl->maxOuterIts = 10;
    ctrl->tol = 0;
    ctrl->spreadFactor = 1e-6;
    ctrl->random = true;
    ctrl->progress = false;
    ElSignCtrlDefault_d( &ctrl->signCtrl );
    return EL_SUCCESS;
}

/* SchurCtrl */
ElError ElSchurCtrlDefault_s( ElSchurCtrl_s* ctrl )
{
    ctrl->useSdc = false;
    ElHessQrCtrlDefault( &ctrl->qrCtrl );
    ElSdcCtrlDefault_s( &ctrl->sdcCtrl );
    return EL_SUCCESS; 
}
ElError ElSchurCtrlDefault_d( ElSchurCtrl_d* ctrl )
{
    ctrl->useSdc = false;
    ElHessQrCtrlDefault( &ctrl->qrCtrl );
    ElSdcCtrlDefault_d( &ctrl->sdcCtrl );
    return EL_SUCCESS; 
}

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* HermitianEig
     ============ */ \
  /* Return all eigenvalues
     ---------------------- */ \
  ElError ElHermitianEig_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElSortType sort ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                          CReflect(sort) ) ) } \
  ElError ElHermitianEigDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElSortType sort ) \
  { EL_TRY( HermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      CReflect(sort) ) ) } \
  /* Return all eigenpairs
     --------------------- */ \
  ElError ElHermitianEigPair_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                          *CReflect(Z), CReflect(sort) ) ) } \
  ElError ElHermitianEigPairDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElDistMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( HermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      DM_CAST(F,Z), CReflect(sort) ) ) } \
  /* Return a subset of eigenvalues 
     ------------------------------ */ \
  ElError ElHermitianEigPartial_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                          CReflect(sort), CReflect(subset) ) ) } \
  ElError ElHermitianEigPartialDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      CReflect(sort), CReflect(subset) ) ) } \
  /* Return a subset of eigenpairs
     ----------------------------- */ \
  ElError ElHermitianEigPairPartial_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                          *CReflect(Z), CReflect(sort), CReflect(subset) ) ) } \
  ElError ElHermitianEigPairPartialDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElDistMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      DM_CAST(F,Z), CReflect(sort), CReflect(subset) ) ) } \
  /* SkewHermitianEig
     ================ */ \
  /* Return all eigenvalues
     ---------------------- */ \
  ElError ElSkewHermitianEig_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElSortType sort ) \
  { EL_TRY( SkewHermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                              CReflect(sort) ) ) } \
  ElError ElSkewHermitianEigDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElSortType sort ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      CReflect(sort) ) ) } \
  /* Return a subset of eigenvalues 
     ------------------------------ */ \
  ElError ElSkewHermitianEigPartial_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( SkewHermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                          CReflect(sort), CReflect(subset) ) ) } \
  ElError ElSkewHermitianEigPartialDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      CReflect(sort), CReflect(subset) ) ) } \
  /* HermitianGenDefEig
     ================== */ \
  /* Return all eigenvalues
     ---------------------- */ \
  ElError ElHermitianGenDefEig_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, ElMatrix_ ## SIGBASE w, \
    ElSortType sort ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w), CReflect(sort) ) ) } \
  ElError ElHermitianGenDefEigDist_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIGBASE w, ElSortType sort ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), DM_CAST(F,A), DM_CAST(F,B), \
      DM_VR_STAR_CAST(Base<F>,w), CReflect(sort) ) ) } \
  /* Return all eigenpairs
     --------------------- */ \
  ElError ElHermitianGenDefEigPair_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w), *CReflect(Z), CReflect(sort) ) ) } \
  ElError ElHermitianGenDefEigPairDist_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), DM_CAST(F,A), DM_CAST(F,B), \
      DM_VR_STAR_CAST(Base<F>,w), DM_CAST(F,Z), CReflect(sort) ) ) } \
  /* Return a subset of eigenvalues 
     ------------------------------ */ \
  ElError ElHermitianGenDefEigPartial_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, ElMatrix_ ## SIGBASE w, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w), CReflect(sort), CReflect(subset) ) ) } \
  ElError ElHermitianGenDefEigPartialDist_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIGBASE w, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), DM_CAST(F,A), DM_CAST(F,B), \
      DM_VR_STAR_CAST(Base<F>,w), CReflect(sort), CReflect(subset) ) ) } \
  /* Return a subset of eigenpairs
     ----------------------------- */ \
  ElError ElHermitianGenDefEigPairPartial_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, \
    ElMatrix_ ## SIGBASE w, ElMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w), *CReflect(Z), CReflect(sort), CReflect(subset) ) ) } \
  ElError ElHermitianGenDefEigPairPartialDist_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), DM_CAST(F,A), DM_CAST(F,B), \
      DM_VR_STAR_CAST(Base<F>,w), DM_CAST(F,Z), CReflect(sort), \
      CReflect(subset) ) ) } \
  /* HermitianTridiagEig
     =================== */ \
  /* Compute all eigenvalues
     ----------------------- */ \
  ElError ElHermitianTridiagEig_ ## SIG \
  ( ElMatrix_ ## SIGBASE d, ElMatrix_ ## SIG e, \
    ElMatrix_ ## SIGBASE w, ElSortType sort ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(e), *CReflect(w), \
      CReflect(sort) ) ) } \
  /* Compute all eigenpairs
     ---------------------- */ \
  ElError ElHermitianTridiagEigPair_ ## SIG \
  ( ElMatrix_ ## SIGBASE d, ElMatrix_ ## SIG e, \
    ElMatrix_ ## SIGBASE w, ElMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(e), *CReflect(w), *CReflect(Z), \
      CReflect(sort) ) ) } \
  /* Compute a subset of eigenvalues
     ------------------------------- */ \
  ElError ElHermitianTridiagEigPartial_ ## SIG \
  ( ElMatrix_ ## SIGBASE d, ElMatrix_ ## SIG e, \
    ElMatrix_ ## SIGBASE w, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(e), *CReflect(w), \
      CReflect(sort), CReflect(subset) ) ) } \
  /* Compute a subset of eigenpairs
     ------------------------------ */ \
  ElError ElHermitianTridiagEigPairPartial_ ## SIG \
  ( ElMatrix_ ## SIGBASE d, ElMatrix_ ## SIG e, \
    ElMatrix_ ## SIGBASE w, ElMatrix_ ## SIG Z, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(e), *CReflect(w), *CReflect(Z), \
      CReflect(sort), CReflect(subset) ) ) } \
  /* Polar decomposition
     =================== */ \
  /* Compute just the polar factor
     ----------------------------- */ \
  ElError ElPolar_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Polar( *CReflect(A) ) ) } \
  ElError ElPolarDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Polar( DM_CAST(F,A) ) ) } \
   ElError ElHermitianPolar_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianPolarDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), DM_CAST(F,A) ) ) } \
  /* Compute the entire polar decomposition
     -------------------------------------- */ \
  ElError ElPolarDecomp_ ## SIG ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG P ) \
  { EL_TRY( Polar( *CReflect(A), *CReflect(P) ) ) } \
  ElError ElPolarDecompDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG P ) \
  { EL_TRY( Polar( DM_CAST(F,A), DM_CAST(F,P) ) ) } \
  ElError ElHermitianPolarDecomp_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIG P ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), *CReflect(A), *CReflect(P) ) ) } \
  ElError ElHermitianPolarDecompDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG P ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), DM_CAST(F,A), DM_CAST(F,P) ) ) } \
  /* Singular Value Decomposition
     ============================ */ \
  /* Compute the singular values
     --------------------------- */ \
  ElError ElSingularValues_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE s ) \
  { EL_TRY( SVD( *CReflect(A), *CReflect(s) ) ) } \
  ElError ElSingularValuesDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE s ) \
  { EL_TRY( SVD( DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,s) ) ) } \
  /* Compute the full SVD
     -------------------- */ \
  ElError ElSVD_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE s, ElMatrix_ ## SIG V ) \
  { EL_TRY( SVD( *CReflect(A), *CReflect(s), *CReflect(V) ) ) } \
  ElError ElSVDDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE s, \
    ElDistMatrix_ ## SIG V ) \
  { EL_TRY( SVD( DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,s), DM_CAST(F,V) ) ) }

#define C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F) \
  /* Schur decomposition
     =================== */ \
  /* Compute the eigenvalues (and possibly Schur factor) 
     --------------------------------------------------- */ \
  ElError ElSchur_ ## SIGBASE \
  ( ElMatrix_ ## SIGBASE A, ElMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( *CReflect(A), *CReflect(w), fullTriangle ) ) } \
  ElError ElSchurDist_ ## SIGBASE \
  ( ElDistMatrix_ ## SIGBASE A, ElDistMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( \
      DM_CAST(Base<F>,A), DM_VR_STAR_CAST(F,w), fullTriangle ) ) } \
  ElError ElSchur_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( *CReflect(A), *CReflect(w), fullTriangle ) ) } \
  ElError ElSchurDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( DM_CAST(F,A), DM_VR_STAR_CAST(F,w), fullTriangle ) ) } \
  /* Compute the eigvalues and Schur vectors (and possibly Schur factor)
     ------------------------------------------------------------------- */ \
  ElError ElSchurDecomp_ ## SIGBASE \
  ( ElMatrix_ ## SIGBASE A, ElMatrix_ ## SIG w, ElMatrix_ ## SIGBASE Q, \
    bool fullTriangle ) \
  { EL_TRY( Schur( \
      *CReflect(A), *CReflect(w), *CReflect(Q), fullTriangle ) ) } \
  ElError ElSchurDecompDist_ ## SIGBASE \
  ( ElDistMatrix_ ## SIGBASE A, ElDistMatrix_ ## SIG w, \
    ElDistMatrix_ ## SIGBASE Q, bool fullTriangle ) \
  { EL_TRY( Schur( \
      DM_CAST(Base<F>,A), DM_VR_STAR_CAST(F,w), DM_CAST(Base<F>,Q), \
      fullTriangle ) ) } \
  ElError ElSchurDecomp_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG w, \
    ElMatrix_ ## SIG Q, bool fullTriangle ) \
  { EL_TRY( Schur( \
      *CReflect(A), *CReflect(w), *CReflect(Q), fullTriangle ) ) } \
  ElError ElSchurDecompDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG w, \
    ElDistMatrix_ ## SIG Q, bool fullTriangle ) \
  { EL_TRY( Schur( \
      DM_CAST(F,A), DM_VR_STAR_CAST(F,w), DM_CAST(F,Q), fullTriangle ) ) }

#define C_PROTO_DOUBLE_ONLY(SIG,SIGBASE,F) \
  /* HermitianTridiagEig
     =================== */ \
  /* Compute all eigenvalues
     ----------------------- */ \
  ElError ElHermitianTridiagEigDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG e, \
    ElDistMatrix_ ## SIGBASE w, ElSortType sort ) \
  { EL_TRY( HermitianTridiagEig( \
      DM_VR_STAR_CAST_CONST(Base<F>,d), DM_VR_STAR_CAST_CONST(F,e), \
      DM_VR_STAR_CAST(Base<F>,w), CReflect(sort) ) ) } \
  /* Compute all eigenpairs
     ---------------------- */ \
  ElError ElHermitianTridiagEigPairDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG e, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( HermitianTridiagEig( \
      DM_VR_STAR_CAST_CONST(Base<F>,d), DM_VR_STAR_CAST_CONST(F,e), \
      DM_VR_STAR_CAST(Base<F>,w), DM_STAR_VR_CAST(F,Z), \
      CReflect(sort) ) ) } \
  /* Compute a subset of eigenvalues
     ------------------------------- */ \
  ElError ElHermitianTridiagEigPartialDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG e, \
    ElDistMatrix_ ## SIGBASE w, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianTridiagEig( \
      DM_VR_STAR_CAST_CONST(Base<F>,d), DM_VR_STAR_CAST_CONST(F,e), \
      DM_VR_STAR_CAST(Base<F>,w), CReflect(sort), CReflect(subset) ) ) } \
  /* Compute a subset of eigenpairs
     ------------------------------ */ \
  ElError ElHermitianTridiagEigPairPartialDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG e, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Z, \
    ElSortType sort, ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( HermitianTridiagEig( \
      DM_VR_STAR_CAST_CONST(Base<F>,d), DM_VR_STAR_CAST_CONST(F,e), \
      DM_VR_STAR_CAST(Base<F>,w), DM_STAR_VR_CAST(F,Z), \
      CReflect(sort), CReflect(subset) ) ) }

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F)

#define C_PROTO_DOUBLE \
  C_PROTO_FIELD(d,d,double) \
  C_PROTO_DOUBLE_ONLY(d,d,double)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F) \
  /* SkewHermitianEig
     ================ */ \
  /* Return all eigenpairs
     --------------------- */ \
  ElError ElSkewHermitianEigPair_ ## SIGBASE \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIGBASE A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), *CReflect(Z), \
      CReflect(sort) ) ) } \
  ElError ElSkewHermitianEigPair_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), *CReflect(Z), \
      CReflect(sort) ) ) } \
  ElError ElSkewHermitianEigPairDist_ ## SIGBASE \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIGBASE A, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), DM_CAST(Base<F>,A), DM_VR_STAR_CAST(Base<F>,w), \
      DM_CAST(F,Z), CReflect(sort) ) ) } \
  ElError ElSkewHermitianEigPairDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElDistMatrix_ ## SIG Z, ElSortType sort ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      DM_CAST(F,Z), CReflect(sort) ) ) } \
  /* Return a subset of eigenpairs
     ----------------------------- */ \
  ElError ElSkewHermitianEigPairPartial_ ## SIGBASE \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIGBASE A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), \
      *CReflect(Z), CReflect(sort), CReflect(subset) ) ) } \
  ElError ElSkewHermitianEigPairPartial_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), \
      *CReflect(Z), CReflect(sort), CReflect(subset) ) ) } \
  ElError ElSkewHermitianEigPairPartialDist_ ## SIGBASE \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIGBASE A, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), DM_CAST(Base<F>,A), DM_VR_STAR_CAST(Base<F>,w), \
      DM_CAST(F,Z), CReflect(sort), CReflect(subset) ) ) } \
  ElError ElSkewHermitianEigPairPartialDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElDistMatrix_ ## SIG Z, ElSortType sort, \
    ElHermitianEigSubset_ ## SIGBASE subset ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), DM_CAST(F,A), DM_VR_STAR_CAST(Base<F>,w), \
      DM_CAST(F,Z), CReflect(sort), CReflect(subset) ) ) }

#define C_PROTO_COMPLEX_DOUBLE \
  C_PROTO_FIELD(z,d,Complex<double>) \
  C_PROTO_DOUBLE_ONLY(z,d,Complex<double>) \
  C_PROTO_COMPLEX_ONLY(z,d,Complex<double>)

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
