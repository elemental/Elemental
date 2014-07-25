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
      DM_CAST(F,Z), CReflect(sort), CReflect(subset) ) ) }


#define C_PROTO_REAL(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* TODO */

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* TODO */

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
