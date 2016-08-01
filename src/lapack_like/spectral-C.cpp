/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <El.h>
using namespace El;

extern "C" {

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

/* HermitianTridiagEigQRCtrl */
ElError ElHermitianTridiagEigQRCtrlDefault( ElHermitianTridiagEigQRCtrl* ctrl )
{
    ctrl->maxIterPerEig = 30;
    ctrl->demandConverged = true;
    ctrl->fullAccuracyTwoByTwo = true;
    return EL_SUCCESS;
}

/* HermitianTridiagEigCtrl */
ElError ElHermitianTridiagEigCtrlDefault_s( ElHermitianTridiagEigCtrl_s* ctrl )
{
    ctrl->wantEigVecs = false;
    ctrl->accumulateEigVecs = false;
    ctrl->sort = EL_ASCENDING;
    ElHermitianEigSubsetDefault_s( &ctrl->subset );
    ctrl->progress = false;
    ctrl->useQR = false;
    ElHermitianTridiagEigQRCtrlDefault( &ctrl->qrCtrl );
    return EL_SUCCESS;
}

ElError ElHermitianTridiagEigCtrlDefault_d( ElHermitianTridiagEigCtrl_d* ctrl )
{
    ctrl->wantEigVecs = false;
    ctrl->accumulateEigVecs = false;
    ctrl->sort = EL_ASCENDING;
    ElHermitianEigSubsetDefault_d( &ctrl->subset );
    ctrl->progress = false;
    ctrl->useQR = false;
    ElHermitianTridiagEigQRCtrlDefault( &ctrl->qrCtrl );
    return EL_SUCCESS;
}

/* HermitianSDCCtrl */
ElError ElHermitianSDCCtrlDefault_s( ElHermitianSDCCtrl_s* ctrl )
{
    ctrl->cutoff = 256;
    ctrl->maxInnerIts = 2;
    ctrl->maxOuterIts = 10;
    ctrl->tol = 0;
    ctrl->spreadFactor = 1e-6f;
    ctrl->progress = false;
    return EL_SUCCESS;
}
ElError ElHermitianSDCCtrlDefault_d( ElHermitianSDCCtrl_d* ctrl )
{
    ctrl->cutoff = 256;
    ctrl->maxInnerIts = 2;
    ctrl->maxOuterIts = 10;
    ctrl->tol = 0;
    ctrl->spreadFactor = 1e-6;
    ctrl->progress = false;
    return EL_SUCCESS;
}

/* HermitianEigCtrl */
ElError ElHermitianEigCtrlDefault_s( ElHermitianEigCtrl_s* ctrl )
{
    ElHermitianTridiagCtrlDefault_s( &ctrl->tridiagCtrl );
    ElHermitianTridiagEigCtrlDefault_s( &ctrl->tridiagEigCtrl );
    ElHermitianSDCCtrlDefault_s( &ctrl->sdcCtrl );
    ctrl->useSDC = false;
    return EL_SUCCESS;
}
ElError ElHermitianEigCtrlDefault_d( ElHermitianEigCtrl_d* ctrl )
{
    ElHermitianTridiagCtrlDefault_d( &ctrl->tridiagCtrl );
    ElHermitianTridiagEigCtrlDefault_d( &ctrl->tridiagEigCtrl );
    ElHermitianSDCCtrlDefault_d( &ctrl->sdcCtrl );
    ctrl->useSDC = false;
    return EL_SUCCESS;
}
ElError ElHermitianEigCtrlDefault_c( ElHermitianEigCtrl_c* ctrl )
{
    ElHermitianTridiagCtrlDefault_c( &ctrl->tridiagCtrl );
    ElHermitianTridiagEigCtrlDefault_s( &ctrl->tridiagEigCtrl );
    ElHermitianSDCCtrlDefault_s( &ctrl->sdcCtrl );
    ctrl->useSDC = false;
    return EL_SUCCESS;
}
ElError ElHermitianEigCtrlDefault_z( ElHermitianEigCtrl_z* ctrl )
{
    ElHermitianTridiagCtrlDefault_z( &ctrl->tridiagCtrl );
    ElHermitianTridiagEigCtrlDefault_d( &ctrl->tridiagEigCtrl );
    ElHermitianSDCCtrlDefault_d( &ctrl->sdcCtrl );
    ctrl->useSDC = false;
    return EL_SUCCESS;
}

/* QDWHCtrl */
ElError ElQDWHCtrlDefault( ElQDWHCtrl* ctrl )
{
    ctrl->colPiv = false;
    ctrl->maxIts = 20;
    return EL_SUCCESS;
}

/* PolarCtrl */
ElError ElPolarCtrlDefault( ElPolarCtrl* ctrl )
{
    ctrl->qdwh = false;
    ElQDWHCtrlDefault( &ctrl->qdwhCtrl );
    return EL_SUCCESS;
}

/* BidiagSVDCtrl */
ElError ElBidiagSVDQRCtrlDefault( ElBidiagSVDQRCtrl* ctrl )
{
    ctrl->maxIterPerVal = 6;
    ctrl->demandConverged = true;
    ctrl->looseMinSingValEst = true;
    ctrl->useFLAME = false;
    ctrl->useLAPACK = false;
    return EL_SUCCESS;
}

ElError ElBidiagSVDCtrlDefault_s( ElBidiagSVDCtrl_s* ctrl )
{
    ctrl->wantU = true;
    ctrl->wantV = true;
    ctrl->accumulateU = false;
    ctrl->accumulateV = false;
    ctrl->approach = EL_THIN_SVD;
    ctrl->tolType = EL_RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl->tol = 0;
    ctrl->progress = false;
    ElBidiagSVDQRCtrlDefault( &ctrl->qrCtrl );
    return EL_SUCCESS;
}
ElError ElBidiagSVDCtrlDefault_d( ElBidiagSVDCtrl_d* ctrl )
{
    ctrl->wantU = true;
    ctrl->wantV = true;
    ctrl->accumulateU = false;
    ctrl->accumulateV = false;
    ctrl->approach = EL_THIN_SVD;
    ctrl->tolType = EL_RELATIVE_TO_MAX_SING_VAL_TOL;
    ctrl->tol = 0;
    ctrl->progress = false;
    ElBidiagSVDQRCtrlDefault( &ctrl->qrCtrl );
    return EL_SUCCESS;
}

/* SVDCtrl */
ElError ElSVDCtrlDefault_s( ElSVDCtrl_s* ctrl )
{
    ctrl->overwrite = false;
    ctrl->time = false;

    ctrl->useLAPACK = false;
    ctrl->useScaLAPACK = false;

    ctrl->valChanRatio = 1.2;
    ctrl->fullChanRatio = 1.5;

    ElBidiagSVDCtrlDefault_s( &ctrl->bidiagSVDCtrl );

    return EL_SUCCESS;
}
ElError ElSVDCtrlDefault_d( ElSVDCtrl_d* ctrl )
{
    ctrl->overwrite = false;
    ctrl->time = false;

    ctrl->useLAPACK = false;
    ctrl->useScaLAPACK = false;

    ctrl->valChanRatio = 1.2;
    ctrl->fullChanRatio = 1.5;

    ElBidiagSVDCtrlDefault_d( &ctrl->bidiagSVDCtrl );

    return EL_SUCCESS;
}

/* HessQRCtrl */
ElError ElHessQRCtrlDefault( ElHessQRCtrl* ctrl )
{
    ctrl->distAED = false;
    ctrl->blockHeight = DefaultBlockHeight();
    ctrl->blockWidth = DefaultBlockWidth();
    return EL_SUCCESS;
}

/* SDCCtrl */
ElError ElSDCCtrlDefault_s( ElSDCCtrl_s* ctrl )
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
ElError ElSDCCtrlDefault_d( ElSDCCtrl_d* ctrl )
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
    ctrl->useSDC = false;
    ElHessQRCtrlDefault( &ctrl->qrCtrl );
    ElSDCCtrlDefault_s( &ctrl->sdcCtrl );
    ctrl->time = false;
    return EL_SUCCESS; 
}
ElError ElSchurCtrlDefault_d( ElSchurCtrl_d* ctrl )
{
    ctrl->useSDC = false;
    ElHessQRCtrlDefault( &ctrl->qrCtrl );
    ElSDCCtrlDefault_d( &ctrl->sdcCtrl );
    ctrl->time = false;
    return EL_SUCCESS; 
}

ElError ElSnapshotCtrlDefault( ElSnapshotCtrl* ctrl )
{
    ctrl->realSize = 0;
    ctrl->imagSize = 0;
    ctrl->imgSaveFreq = -1;
    ctrl->numSaveFreq = -1;
    ctrl->imgDispFreq = -1;
    ctrl->imgSaveCount = 0;
    ctrl->numSaveCount = 0;
    ctrl->imgDispCount = 0;
    ctrl->imgBase = "ps";
    ctrl->numBase = "ps";
    ctrl->imgFormat = EL_PNG;
    ctrl->numFormat = EL_ASCII_MATLAB;
    ctrl->itCounts = true;
    return EL_SUCCESS;
}
ElError ElSnapshotCtrlDestroy( const ElSnapshotCtrl* ctrl )
{
    delete ctrl->imgBase;
    delete ctrl->numBase;
    delete ctrl;
    return EL_SUCCESS;
}

ElError ElPseudospecCtrlDefault_s( ElPseudospecCtrl_s* ctrl )
{
    ctrl->norm = EL_PS_TWO_NORM;
    ctrl->blockWidth = 10;
    ctrl->schur = true;
    ctrl->forceComplexSchur = false;
    ctrl->forceComplexPs = false;
    ElSchurCtrlDefault_s( &ctrl->schurCtrl );
    ctrl->maxIts = 50;
    ctrl->tol = 1e-6f;
    ctrl->deflate = true;
    ctrl->arnoldi = true;
    ctrl->basisSize = 10;
    ctrl->reorthog = true;
    ctrl->progress = false;
    ElSnapshotCtrlDefault( &ctrl->snapCtrl );
    return EL_SUCCESS;
}
ElError ElPseudospecCtrlDestroy_s( const ElPseudospecCtrl_s* ctrl )
{
    ElSnapshotCtrlDestroy( &ctrl->snapCtrl );
    delete ctrl;
    return EL_SUCCESS;
}

ElError ElPseudospecCtrlDefault_d( ElPseudospecCtrl_d* ctrl )
{
    ctrl->norm = EL_PS_TWO_NORM;
    ctrl->blockWidth = 10;
    ctrl->schur = true;
    ctrl->forceComplexSchur = false;
    ctrl->forceComplexPs = false;
    ElSchurCtrlDefault_d( &ctrl->schurCtrl );
    ctrl->maxIts = 50;
    ctrl->tol = 1e-6;
    ctrl->deflate = true;
    ctrl->arnoldi = true;
    ctrl->basisSize = 10;
    ctrl->reorthog = true;
    ctrl->progress = false;
    ElSnapshotCtrlDefault( &ctrl->snapCtrl );
    return EL_SUCCESS;
}
ElError ElPseudospecCtrlDestroy_d( const ElPseudospecCtrl_d* ctrl )
{
    ElSnapshotCtrlDestroy( &ctrl->snapCtrl );
    delete ctrl;
    return EL_SUCCESS;
}

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* HermitianEig
     ============ */ \
  /* Return all eigenvalues
     ---------------------- */ \
  ElError ElHermitianEig_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w) ) ) } \
  ElError ElHermitianEigDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w) ) ) } \
  /* Return all eigenpairs
     --------------------- */ \
  ElError ElHermitianEigPair_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Q ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                          *CReflect(Q) ) ) } \
  ElError ElHermitianEigPairDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElDistMatrix_ ## SIG Q ) \
  { EL_TRY( HermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w), \
                          *CReflect(Q) ) ) } \
  /* TODO: Expert versions */ \
  /* SkewHermitianEig
     ================ */ \
  /* Return all eigenvalues
     ---------------------- */ \
  ElError ElSkewHermitianEig_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w ) \
  { EL_TRY( SkewHermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w) ) ) } \
  ElError ElSkewHermitianEigDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w ) \
  { EL_TRY( SkewHermitianEig( CReflect(uplo), *CReflect(A), *CReflect(w) ) ) } \
  /* TODO: Expert versions */ \
  /* HermitianGenDefEig
     ================== */ \
  /* Return all eigenvalues
     ---------------------- */ \
  ElError ElHermitianGenDefEig_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, ElMatrix_ ## SIGBASE w ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w) ) ) } \
  ElError ElHermitianGenDefEigDist_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIGBASE w ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w) ) ) } \
  /* Return all eigenpairs
     --------------------- */ \
  ElError ElHermitianGenDefEigPair_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Q ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w), *CReflect(Q) ) ) } \
  ElError ElHermitianGenDefEigPairDist_ ## SIG \
  ( ElPencil pencil, ElUpperOrLower uplo, \
    ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Q ) \
  { EL_TRY( HermitianGenDefEig( \
      CReflect(pencil), CReflect(uplo), *CReflect(A), *CReflect(B), \
      *CReflect(w), *CReflect(Q) ) ) } \
  /* TODO: Expert versions */ \
  /* HermitianTridiagEig
     =================== */ \
  /* Compute all eigenvalues
     ----------------------- */ \
  ElError ElHermitianTridiagEig_ ## SIG \
  ( ElMatrix_ ## SIGBASE d, ElMatrix_ ## SIG dSub, \
    ElMatrix_ ## SIGBASE w ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(dSub), *CReflect(w) ) ) } \
  ElError ElHermitianTridiagEigDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG dSub, \
    ElDistMatrix_ ## SIGBASE w ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(dSub), *CReflect(w) ) ) } \
  /* Compute all eigenpairs
     ---------------------- */ \
  ElError ElHermitianTridiagEigPair_ ## SIG \
  ( ElMatrix_ ## SIGBASE d, ElMatrix_ ## SIG dSub, \
    ElMatrix_ ## SIGBASE w, ElMatrix_ ## SIG Q ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(dSub), *CReflect(w), *CReflect(Q) ) ) } \
  ElError ElHermitianTridiagEigPairDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIGBASE d, ElConstDistMatrix_ ## SIG dSub, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Q ) \
  { EL_TRY( HermitianTridiagEig( \
      *CReflect(d), *CReflect(dSub), *CReflect(w), *CReflect(Q) ) ) } \
  /* TODO: Expert versions */ \
  /* Polar decomposition
     =================== */ \
  /* Compute just the polar factor
     ----------------------------- */ \
  ElError ElPolar_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Polar( *CReflect(A) ) ) } \
  ElError ElPolarDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Polar( *CReflect(A) ) ) } \
   ElError ElHermitianPolar_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElHermitianPolarDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), *CReflect(A) ) ) } \
  /* Compute the entire polar decomposition
     -------------------------------------- */ \
  ElError ElPolarDecomp_ ## SIG ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG P ) \
  { EL_TRY( Polar( *CReflect(A), *CReflect(P) ) ) } \
  ElError ElPolarDecompDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG P ) \
  { EL_TRY( Polar( *CReflect(A), *CReflect(P) ) ) } \
  ElError ElHermitianPolarDecomp_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIG P ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), *CReflect(A), *CReflect(P) ) ) } \
  ElError ElHermitianPolarDecompDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG P ) \
  { EL_TRY( HermitianPolar( CReflect(uplo), *CReflect(A), *CReflect(P) ) ) } \
  /* Singular Value Decomposition
     ============================ */ \
  /* Compute the singular values
     --------------------------- */ \
  ElError ElSingularValues_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIGBASE s ) \
  { EL_TRY( SVD( *CReflect(A), *CReflect(s) ) ) } \
  ElError ElSingularValuesDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIGBASE s ) \
  { EL_TRY( SVD( *CReflect(A), *CReflect(s) ) ) } \
  /* Expert versions
     ^^^^^^^^^^^^^^^ */ \
  ElError ElSingularValuesXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIGBASE s, \
    ElSVDCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SVD( *CReflect(A), *CReflect(s), CReflect(ctrl) ) ) } \
  ElError ElTSQRSingularValues_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIGBASE s ) \
  { EL_TRY( svd::TSQR( *CReflect(A), *CReflect(s) ) ) } \
  /* Compute the full SVD
     -------------------- */ \
  ElError ElSVD_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG U, \
    ElMatrix_ ## SIGBASE s, \
    ElMatrix_ ## SIG V ) \
  { EL_TRY( SVD( *CReflect(A), *CReflect(U), *CReflect(s), *CReflect(V) ) ) } \
  ElError ElSVDDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG U, \
    ElDistMatrix_ ## SIGBASE s, \
    ElDistMatrix_ ## SIG V ) \
  { EL_TRY( SVD( *CReflect(A), *CReflect(U), *CReflect(s), *CReflect(V) ) ) } \
  /* Expert versions
     ^^^^^^^^^^^^^^^ */ \
  ElError ElSVDX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG U, \
    ElMatrix_ ## SIGBASE s, \
    ElMatrix_ ## SIG V, \
    ElSVDCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY(SVD( *CReflect(A), *CReflect(U), *CReflect(s), *CReflect(V), \
      CReflect(ctrl) )) } \
  ElError ElSVDXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG U, \
    ElDistMatrix_ ## SIGBASE s, \
    ElDistMatrix_ ## SIG V, \
    ElSVDCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY(SVD( *CReflect(A), *CReflect(U), *CReflect(s), *CReflect(V), \
      CReflect(ctrl) )) } \
  ElError ElTSQRSVD_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG U, \
    ElDistMatrix_ ## SIGBASE s, \
    ElDistMatrix_ ## SIG V ) \
  { EL_TRY( svd::TSQR( *CReflect(A), \
      *CReflect(U), *CReflect(s), *CReflect(V) ) ) } \
  /* Hermitian Singular Value Decomposition
     ====================================== */ \
  /* Compute the singular values
     --------------------------- */ \
  ElError ElHermitianSingularValues_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIGBASE s ) \
  { EL_TRY( HermitianSVD( CReflect(uplo), *CReflect(A), *CReflect(s) ) ) } \
  ElError ElHermitianSingularValuesDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIGBASE s ) \
  { EL_TRY( HermitianSVD( CReflect(uplo), *CReflect(A), *CReflect(s) ) ) } \
  /* Compute the full SVD
     -------------------- */ \
  ElError ElHermitianSVD_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG U, \
    ElMatrix_ ## SIGBASE s, \
    ElMatrix_ ## SIG V ) \
  { EL_TRY( HermitianSVD( CReflect(uplo), *CReflect(A), \
      *CReflect(U), *CReflect(s), *CReflect(V) ) ) } \
  ElError ElHermitianSVDDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG U, \
    ElDistMatrix_ ## SIGBASE s, \
    ElDistMatrix_ ## SIG V ) \
  { EL_TRY( HermitianSVD( CReflect(uplo), *CReflect(A), \
      *CReflect(U), *CReflect(s), *CReflect(V) ) ) } \
  /* Product Lanczos
     =============== */ \
  ElError ElProductLanczosSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIGBASE T, \
    ElInt basisSize ) \
  { EL_TRY( ProductLanczos( *CReflect(A), *CReflect(T), basisSize ) ) } \
  ElError ElProductLanczosDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE T, \
    ElInt basisSize ) \
  { EL_TRY( ProductLanczos( *CReflect(A), *CReflect(T), basisSize ) ) } \
  ElError ElProductLanczosDecompSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG V, \
    ElMatrix_ ## SIGBASE T,        ElMatrix_ ## SIG v, \
    Base<F>* beta, ElInt basisSize ) \
  { EL_TRY( *beta = ProductLanczosDecomp \
                    ( *CReflect(A), *CReflect(V), \
                      *CReflect(T), *CReflect(v), basisSize ) ) } \
  ElError ElProductLanczosDecompDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG V, \
    ElDistMatrix_ ## SIGBASE T,        ElDistMultiVec_ ## SIG v, \
    Base<F>* beta, ElInt basisSize ) \
  { EL_TRY( *beta = ProductLanczosDecomp \
                    ( *CReflect(A), *CReflect(V), \
                      *CReflect(T), *CReflect(v), basisSize ) ) } \
  /* Extremal singular value estimation
     ================================== */ \
  ElError ElExtremalSingValEstSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElInt basisSize, \
    Base<F>* sigMin, Base<F>* sigMax ) \
  { EL_TRY( auto extremal = ExtremalSingValEst( *CReflect(A), basisSize ); \
            *sigMin = extremal.first; *sigMax = extremal.second ) } \
  ElError ElExtremalSingValEstDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElInt basisSize, \
    Base<F>* sigMin, Base<F>* sigMax ) \
  { EL_TRY( auto extremal = ExtremalSingValEst( *CReflect(A), basisSize ); \
            *sigMin = extremal.first; *sigMax = extremal.second ) } \
  /* Triangular eigenvectors
     ======================= */ \
  ElError ElTriangEig_ ## SIG \
  ( ElMatrix_ ## SIG U, ElMatrix_ ## SIG X ) \
  { EL_TRY( TriangEig( *CReflect(U), *CReflect(X) ) ) } \
  ElError ElTriangEigDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG U, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( TriangEig( *CReflect(U), *CReflect(X) ) ) }

#define C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F) \
  /* Eigenvalue decomposition
     ======================== */ \
  ElError ElEig_ ## SIGBASE \
  ( ElMatrix_ ## SIGBASE A, \
    ElMatrix_ ## SIG w, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( Eig( *CReflect(A), *CReflect(w), *CReflect(X) ) ) } \
  ElError ElEig_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG w, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( Eig( *CReflect(A), *CReflect(w), *CReflect(X) ) ) } \
  ElError ElEigDist_ ## SIGBASE \
  ( ElDistMatrix_ ## SIGBASE A, \
    ElDistMatrix_ ## SIG w, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( Eig( *CReflect(A), *CReflect(w), *CReflect(X) ) ) } \
  ElError ElEigDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG w, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( Eig( *CReflect(A), *CReflect(w), *CReflect(X) ) ) } \
  /* Schur decomposition
     =================== */ \
  /* Compute the eigenvalues (and possibly Schur factor) 
     --------------------------------------------------- */ \
  ElError ElSchur_ ## SIGBASE \
  ( ElMatrix_ ## SIGBASE A, ElMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( *CReflect(A), *CReflect(w), fullTriangle ) ) } \
  ElError ElSchurDist_ ## SIGBASE \
  ( ElDistMatrix_ ## SIGBASE A, ElDistMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( *CReflect(A), *CReflect(w), fullTriangle ) ) } \
  ElError ElSchur_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( *CReflect(A), *CReflect(w), fullTriangle ) ) } \
  ElError ElSchurDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG w, bool fullTriangle ) \
  { EL_TRY( Schur( *CReflect(A), *CReflect(w), fullTriangle ) ) } \
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
      *CReflect(A), *CReflect(w), *CReflect(Q), fullTriangle ) ) } \
  ElError ElSchurDecomp_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG w, \
    ElMatrix_ ## SIG Q, bool fullTriangle ) \
  { EL_TRY( Schur( \
      *CReflect(A), *CReflect(w), *CReflect(Q), fullTriangle ) ) } \
  ElError ElSchurDecompDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG w, \
    ElDistMatrix_ ## SIG Q, bool fullTriangle ) \
  { EL_TRY( Schur( *CReflect(A), *CReflect(w), *CReflect(Q), fullTriangle ) ) }\
  /* SkewHermitianEig
     ================ */ \
  /* Return all eigenpairs
     --------------------- */ \
  ElError ElSkewHermitianEigPair_ ## SIGBASE \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIGBASE A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Q ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), *CReflect(Q) ) ) } \
  ElError ElSkewHermitianEigPair_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_ ## SIGBASE w, \
    ElMatrix_ ## SIG Q ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), *CReflect(Q) ) ) } \
  ElError ElSkewHermitianEigPairDist_ ## SIGBASE \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIGBASE A, \
    ElDistMatrix_ ## SIGBASE w, ElDistMatrix_ ## SIG Q ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), *CReflect(Q) ) ) } \
  ElError ElSkewHermitianEigPairDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE w, \
    ElDistMatrix_ ## SIG Q ) \
  { EL_TRY( SkewHermitianEig( \
      CReflect(uplo), *CReflect(A), *CReflect(w), *CReflect(Q) ) ) } \
  /* Pseudospectra
     ============= */ \
  /* (Pseudo-)Spectral portrait
     -------------------------- */ \
  ElError ElSpectralPortrait_ ## SIGBASE \
  ( ElConstMatrix_ ## SIGBASE A, ElMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( *CReflect(A), *CReflect(invNormMap), \
        realSize, imagSize, box ); \
      *boxC = CReflect(box) ) } \
  ElError ElSpectralPortrait_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( *CReflect(A), *CReflect(invNormMap), \
        realSize, imagSize, box ); \
      *boxC = CReflect(box) ) } \
  ElError ElSpectralPortraitDist_ ## SIGBASE \
  ( ElConstDistMatrix_ ## SIGBASE A, ElDistMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( *CReflect(A), *CReflect(invNormMap), \
        realSize, imagSize, box ); \
      *boxC = CReflect(box) ) } \
  ElError ElSpectralPortraitDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( *CReflect(A), *CReflect(invNormMap), \
        realSize, imagSize, box ); \
      *boxC = CReflect(box) ) } \
  /* Expert version */ \
  ElError ElSpectralPortraitX_ ## SIGBASE \
  ( ElConstMatrix_ ## SIGBASE A, ElMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC, \
    ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( \
        *CReflect(A), *CReflect(invNormMap), realSize, imagSize, box, \
        CReflect(ctrl) ); \
      *boxC = CReflect(box) ) } \
  ElError ElSpectralPortraitX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC, \
    ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( \
        *CReflect(A), *CReflect(invNormMap), realSize, imagSize, box, \
        CReflect(ctrl) ); \
      *boxC = CReflect(box) ) } \
  ElError ElSpectralPortraitXDist_ ## SIGBASE \
  ( ElConstDistMatrix_ ## SIGBASE A, ElDistMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC, \
    ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( \
        *CReflect(A), *CReflect(invNormMap), realSize, imagSize, box, \
        CReflect(ctrl) ); \
      *boxC = CReflect(box) ) } \
  ElError ElSpectralPortraitXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE invNormMap, \
    ElInt realSize, ElInt imagSize, ElSpectralBox_ ## SIGBASE* boxC, \
    ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( \
      SpectralBox<Base<F>> box; \
      SpectralPortrait( \
        *CReflect(A), *CReflect(invNormMap), realSize, imagSize, box, \
        CReflect(ctrl) ); \
      *boxC = CReflect(box) ) } \
  /* (Pseudo-)Spectral window
     ------------------------ */ \
  ElError ElSpectralWindow_ ## SIGBASE \
  ( ElConstMatrix_ ## SIGBASE A, ElMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, realSize, imagSize ) ) } \
  ElError ElSpectralWindow_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, realSize, imagSize ) ) } \
  ElError ElSpectralWindowDist_ ## SIGBASE \
  ( ElConstDistMatrix_ ## SIGBASE A, ElDistMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, realSize, imagSize ) ) } \
  ElError ElSpectralWindowDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, realSize, imagSize ) ) } \
  /* Expert version */ \
  ElError ElSpectralWindowX_ ## SIGBASE \
  ( ElConstMatrix_ ## SIGBASE A, ElMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, \
      realSize, imagSize, CReflect(ctrl) ) ) } \
  ElError ElSpectralWindowX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, \
      realSize, imagSize, CReflect(ctrl) ) ) } \
  ElError ElSpectralWindowXDist_ ## SIGBASE \
  ( ElConstDistMatrix_ ## SIGBASE A, ElDistMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, \
      realSize, imagSize, CReflect(ctrl) ) ) } \
  ElError ElSpectralWindowXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE invNormMap, \
    CREFLECT(F) center, Base<F> realWidth, Base<F> imagWidth, \
    ElInt realSize, ElInt imagSize, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralWindow( \
      *CReflect(A), *CReflect(invNormMap), \
      CReflect(center), realWidth, imagWidth, \
      realSize, imagSize, CReflect(ctrl) ) ) } \
  /* (Pseudo-)Spectral Cloud
     ----------------------- */ \
  ElError ElSpectralCloud_ ## SIGBASE \
  ( ElConstMatrix_ ## SIGBASE A, ElConstMatrix_ ## SIG shifts, \
    ElMatrix_ ## SIGBASE invNorms ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms) ) ) } \
  ElError ElSpectralCloud_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG shifts, \
    ElMatrix_ ## SIGBASE invNorms ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms) ) ) } \
  ElError ElSpectralCloudDist_ ## SIGBASE \
  ( ElConstDistMatrix_ ## SIGBASE A, ElConstDistMatrix_ ## SIG shifts, \
    ElDistMatrix_ ## SIGBASE invNorms ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms) ) ) } \
  ElError ElSpectralCloudDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG shifts, \
    ElDistMatrix_ ## SIGBASE invNorms ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms) ) ) } \
  /* Expert version */ \
  ElError ElSpectralCloudX_ ## SIGBASE \
  ( ElConstMatrix_ ## SIGBASE A, ElConstMatrix_ ## SIG shifts, \
    ElMatrix_ ## SIGBASE invNorms, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms), \
      CReflect(ctrl) ) ) } \
  ElError ElSpectralCloudX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG shifts, \
    ElMatrix_ ## SIGBASE invNorms, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms), \
      CReflect(ctrl) ) ) } \
  ElError ElSpectralCloudXDist_ ## SIGBASE \
  ( ElConstDistMatrix_ ## SIGBASE A, ElConstDistMatrix_ ## SIG shifts, \
    ElDistMatrix_ ## SIGBASE invNorms, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms), \
      CReflect(ctrl) ) ) } \
  ElError ElSpectralCloudXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG shifts, \
    ElDistMatrix_ ## SIGBASE invNorms, ElPseudospecCtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( SpectralCloud( \
      *CReflect(A), *CReflect(shifts), *CReflect(invNorms), \
      CReflect(ctrl) ) ) }

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F) \

#define EL_NO_INT_PROTO
#include <El/macros/CInstantiate.h>

} // extern "C"
