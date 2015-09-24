/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

/* Basis Pursuit
   ============= */
ElError ElBPADMMCtrlDefault_s( ElBPADMMCtrl_s* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->usePinv = false;
    ctrl->pinvTol = 0;
    ctrl->progress = true;
    return EL_SUCCESS;
}

ElError ElBPADMMCtrlDefault_d( ElBPADMMCtrl_d* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->usePinv = false;
    ctrl->pinvTol = 0;
    ctrl->progress = true;
    return EL_SUCCESS;
}

ElError ElBPCtrlDefault_s( ElBPCtrl_s* ctrl, bool isSparse )
{
    ctrl->useIPM = true;
    ctrl->useSOCP = false;
    ElBPADMMCtrlDefault_s( &ctrl->admmCtrl );
    ElLPDirectCtrlDefault_s( &ctrl->lpIPMCtrl, isSparse );
    ElSOCPDirectCtrlDefault_s( &ctrl->socpIPMCtrl );
    return EL_SUCCESS;
}

ElError ElBPCtrlDefault_d( ElBPCtrl_d* ctrl, bool isSparse )
{
    ctrl->useIPM = true;
    ctrl->useSOCP = false;
    ElBPADMMCtrlDefault_d( &ctrl->admmCtrl );
    ElLPDirectCtrlDefault_d( &ctrl->lpIPMCtrl, isSparse );
    ElSOCPDirectCtrlDefault_d( &ctrl->socpIPMCtrl );
    return EL_SUCCESS;
}

ElError ElBPCtrlDefault_c( ElBPCtrl_c* ctrl )
{
    ElSOCPDirectCtrlDefault_s( &ctrl->ipmCtrl );
    return EL_SUCCESS;
}

ElError ElBPCtrlDefault_z( ElBPCtrl_z* ctrl )
{
    ElSOCPDirectCtrlDefault_d( &ctrl->ipmCtrl );
    return EL_SUCCESS;
}

/* Basis Pursuit Denoising / LASSO
   =============================== */
ElError ElBPDNADMMCtrlDefault_s( ElBPDNADMMCtrl_s* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->inv = true;
    ctrl->progress = true;
    return EL_SUCCESS;
}

ElError ElBPDNADMMCtrlDefault_d( ElBPDNADMMCtrl_d* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->inv = true;
    ctrl->progress = true;
    return EL_SUCCESS;
}

ElError ElBPDNCtrlDefault_s( ElBPDNCtrl_s* ctrl )
{
    ctrl->useIPM = true;
    ElBPDNADMMCtrlDefault_s( &ctrl->admmCtrl );
    ElQPAffineCtrlDefault_s( &ctrl->ipmCtrl );
    return EL_SUCCESS;
}

ElError ElBPDNCtrlDefault_d( ElBPDNCtrl_d* ctrl )
{
    ctrl->useIPM = true;
    ElBPDNADMMCtrlDefault_d( &ctrl->admmCtrl );
    ElQPAffineCtrlDefault_d( &ctrl->ipmCtrl );
    return EL_SUCCESS;
}

/* Non-negative Least Squares
   ========================== */
ElError ElNNLSCtrlDefault_s( ElNNLSCtrl_s* ctrl )
{
    ctrl->approach = EL_NNLS_SOCP;
    ElADMMCtrlDefault_s( &ctrl->admmCtrl );
    ElQPDirectCtrlDefault_s( &ctrl->qpCtrl );
    ElSOCPAffineCtrlDefault_s( &ctrl->socpCtrl );
    return EL_SUCCESS;
}

ElError ElNNLSCtrlDefault_d( ElNNLSCtrl_d* ctrl )
{
    ctrl->approach = EL_NNLS_SOCP;
    ElADMMCtrlDefault_d( &ctrl->admmCtrl );
    ElQPDirectCtrlDefault_d( &ctrl->qpCtrl );
    ElSOCPAffineCtrlDefault_d( &ctrl->socpCtrl );
    return EL_SUCCESS;
}

/* Non-negative Matrix factorization
   ================================= */
ElError ElNMFCtrlDefault_s( ElNMFCtrl_s* ctrl )
{
    ElNNLSCtrlDefault_s( &ctrl->nnlsCtrl );
    ctrl->maxIter = 20;
    return EL_SUCCESS;
}

ElError ElNMFCtrlDefault_d( ElNMFCtrl_d* ctrl )
{
    ElNNLSCtrlDefault_d( &ctrl->nnlsCtrl );
    ctrl->maxIter = 20;
    return EL_SUCCESS;
}

/* Robust Principal Component Analysis
   =================================== */
ElError ElRPCACtrlDefault_s( ElRPCACtrl_s* ctrl )
{
    ctrl->useALM = true;
    ctrl->usePivQR = false;
    ctrl->progress = true;
    ctrl->numPivSteps = 7;
    ctrl->maxIts = 1000;
    ctrl->tau = 0;
    ctrl->beta = 1;
    ctrl->rho = 6;
    ctrl->tol = 1e-5;
    return EL_SUCCESS;
}

ElError ElRPCACtrlDefault_d( ElRPCACtrl_d* ctrl )
{
    ctrl->useALM = true;
    ctrl->usePivQR = false;
    ctrl->progress = true;
    ctrl->numPivSteps = 7;
    ctrl->maxIts = 1000;
    ctrl->tau = 0;
    ctrl->beta = 1;
    ctrl->rho = 6;
    ctrl->tol = 1e-5;
    return EL_SUCCESS;
}

/* Sparse Inverse Covariance Selection
   =================================== */
ElError ElSparseInvCovCtrlDefault_s( ElSparseInvCovCtrl_s* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->progress = true;
    return EL_SUCCESS;
}

ElError ElSparseInvCovCtrlDefault_d( ElSparseInvCovCtrl_d* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->progress = true;
    return EL_SUCCESS;
}

/* ModelFit
   ======== */
ElError ElModelFitCtrlDefault_s( ElModelFitCtrl_s* ctrl )
{
    ctrl->rho = 1;
    ctrl->maxIter = 500;
    ctrl->inv = true;
    ctrl->progress = true;
    return EL_SUCCESS;
}

ElError ElModelFitCtrlDefault_d( ElModelFitCtrl_d* ctrl )
{
    ctrl->rho = 1;
    ctrl->maxIter = 500;
    ctrl->inv = true;
    ctrl->progress = true;
    return EL_SUCCESS;
}

/* Support Vector Machine
   ====================== */
ElError ElSVMCtrlDefault_s( ElSVMCtrl_s* ctrl )
{
    ctrl->useIPM = true;
    ElModelFitCtrlDefault_s( &ctrl->modelFitCtrl );
    ElQPAffineCtrlDefault_s( &ctrl->ipmCtrl );
    return EL_SUCCESS;
}

ElError ElSVMCtrlDefault_d( ElSVMCtrl_d* ctrl )
{
    ctrl->useIPM = true;
    ElModelFitCtrlDefault_d( &ctrl->modelFitCtrl );
    ElQPAffineCtrlDefault_d( &ctrl->ipmCtrl );
    return EL_SUCCESS;
}

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Basis Pursuit
     ============= */ \
  ElError ElBP_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElBPDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElBPSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElBPDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( BP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElBPX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElBPCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x, ElBPCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElBPCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x, ElBPCtrl_ ## SIG ctrl ) \
  { EL_TRY( BP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  /* Robust Principal Component Analysis
     =================================== */ \
  ElError ElRPCA_ ## SIG \
  ( ElConstMatrix_ ## SIG M, ElMatrix_ ## SIG L, ElMatrix_ ## SIG S ) \
  { EL_TRY( RPCA( *CReflect(M), *CReflect(L), *CReflect(S) ) ) } \
  ElError ElRPCADist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG M, \
    ElDistMatrix_ ## SIG L, ElDistMatrix_ ## SIG S ) \
  { EL_TRY( RPCA( *CReflect(M), *CReflect(L), *CReflect(S) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElRPCAX_ ## SIG \
  ( ElConstMatrix_ ## SIG M, ElMatrix_ ## SIG L, ElMatrix_ ## SIG S, \
    ElRPCACtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( \
      RPCA( *CReflect(M), *CReflect(L), *CReflect(S), CReflect(ctrl) ) ) } \
  ElError ElRPCAXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG M, \
    ElDistMatrix_ ## SIG L, ElDistMatrix_ ## SIG S, \
    ElRPCACtrl_ ## SIGBASE ctrl ) \
  { EL_TRY( \
      RPCA( *CReflect(M), *CReflect(L), *CReflect(S), CReflect(ctrl) ) ) } \
  /* Sparse inverse covariance selection
     =================================== */ \
  ElError ElSparseInvCov_ ## SIG \
  ( ElConstMatrix_ ## SIG D, Base<F> lambda, ElMatrix_ ## SIG Z, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = SparseInvCov( *CReflect(D), lambda, *CReflect(Z) ) ) } \
  ElError ElSparseInvCovDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG D, Base<F> lambda, ElDistMatrix_ ## SIG Z, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = SparseInvCov( *CReflect(D), lambda, *CReflect(Z) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElSparseInvCovX_ ## SIG \
  ( ElConstMatrix_ ## SIG D, Base<F> lambda, ElMatrix_ ## SIG Z, \
    ElSparseInvCovCtrl_ ## SIGBASE ctrl, ElInt* numIts ) \
  { EL_TRY( *numIts = \
      SparseInvCov( *CReflect(D), lambda, *CReflect(Z), CReflect(ctrl) ) ) } \
  ElError ElSparseInvCovXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG D, Base<F> lambda, ElDistMatrix_ ## SIG Z, \
    ElSparseInvCovCtrl_ ## SIGBASE ctrl, ElInt* numIts ) \
  { EL_TRY( *numIts = \
      SparseInvCov( *CReflect(D), lambda, *CReflect(Z), CReflect(ctrl) ) ) }

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* Chebyshev point
     =============== */ \
  ElError ElCP_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElCPDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElCPSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElCPDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( CP( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElCPX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElCPXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElCPXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElCPXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( CP \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  /* Dantzig Selector
     ================ */ \
  ElError ElDS_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElDSDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElDSSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElDSDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( DS( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElDSX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElDSXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElDSXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElDSXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( DS \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Least Absolute Value regression
     =============================== */ \
  ElError ElLAV_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElLAVDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElLAVSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  ElError ElLAVDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( LAV( *CReflect(A), *CReflect(b), *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElLAVX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElLAVXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElLAVXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElLAVXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElDistMultiVec_ ## SIG x, ElLPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LAV \
      ( *CReflect(A), *CReflect(b), *CReflect(x), CReflect(ctrl) ) ) } \
  /* Basis Pursuit Denoising
     ======================= */ \
  ElError ElBPDN_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElBPDNDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElBPDNSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElBPDNDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( BPDN( *CReflect(A), *CReflect(b), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElBPDNX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElBPDNCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPDNXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda, ElDistMatrix_ ## SIG x, ElBPDNCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPDNXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda, ElMatrix_ ## SIG x, ElBPDNCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElBPDNXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda, ElDistMultiVec_ ## SIG x, ElBPDNCtrl_ ## SIG ctrl ) \
  { EL_TRY( BPDN \
      ( *CReflect(A), *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Elastic net
     =========== */ \
  ElError ElEN_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  ElError ElENDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  ElError ElENSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  ElError ElENDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( EN( *CReflect(A), *CReflect(b), \
      lambda1, lambda2, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElENX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElENXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElENXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real lambda1, Real lambda2, ElMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElENXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real lambda1, Real lambda2, ElDistMultiVec_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( EN \
      ( *CReflect(A), *CReflect(b), lambda1, lambda2, \
        *CReflect(x), CReflect(ctrl) ) ) } \
  /* Support Vector Machine (soft-margin) 
     ==================================== */ \
  ElError ElSVM_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  ElError ElSVMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG d, \
    Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  ElError ElSVMSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  ElError ElSVMDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG d, \
    Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( SVM( *CReflect(A), *CReflect(d), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElSVMX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x, ElSVMCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElSVMXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG d, \
    Real lambda, ElDistMatrix_ ## SIG x, ElSVMCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElSVMXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG d, \
    Real lambda, ElMatrix_ ## SIG x, ElSVMCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElSVMXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG d, \
    Real lambda, ElDistMultiVec_ ## SIG x, ElSVMCtrl_ ## SIG ctrl ) \
  { EL_TRY( SVM \
      ( *CReflect(A), *CReflect(d), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Total variation denoising 
     ========================= */ \
  ElError ElTV_ ## SIG \
  ( ElConstMatrix_ ## SIG b, Real lambda, ElMatrix_ ## SIG x ) \
  { EL_TRY( TV( *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElTVDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG b, Real lambda, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( TV( *CReflect(b), lambda, *CReflect(x) ) ) } \
  ElError ElTVDistSparse_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG b, Real lambda, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( TV( *CReflect(b), lambda, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElTVX_ ## SIG \
  ( ElConstMatrix_ ## SIG b, Real lambda, ElMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( TV \
      ( *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElTVXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG b, Real lambda, ElDistMatrix_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( TV \
      ( *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElTVXDistSparse_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG b, Real lambda, ElDistMultiVec_ ## SIG x, \
    ElQPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( TV \
      ( *CReflect(b), lambda, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Logistic regression
     =================== */ \
  ElError ElLogisticRegression_ ## SIG \
  ( ElConstMatrix_ ## SIG G, ElConstMatrix_ ## SIG q, \
    ElMatrix_ ## SIG z, Real gamma, ElRegularization penalty, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = LogisticRegression( *CReflect(G), *CReflect(q), \
      *CReflect(z), gamma, CReflect(penalty) ) ) } \
  ElError ElLogisticRegressionDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG G, ElConstDistMatrix_ ## SIG q, \
    ElDistMatrix_ ## SIG z, Real gamma, ElRegularization penalty, \
    ElInt* numIts ) \
  { EL_TRY( *numIts = LogisticRegression( *CReflect(G), *CReflect(q), \
      *CReflect(z), gamma, CReflect(penalty) ) ) } \
  /* ModelFit 
     ======== */ \
  ElError ElModelFit_ ## SIG \
  ( void (*lossProx)(ElMatrix_ ## SIG,Real), \
    void (*regProx)(ElMatrix_ ## SIG,Real), \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG w, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( Matrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( Matrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( function<void(Matrix<Real>&,Real)>(lossLambda), \
          function<void(Matrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElModelFitDist_ ## SIG \
  ( void (*lossProx)(ElDistMatrix_ ## SIG,Real), \
    void (*regProx)(ElDistMatrix_ ## SIG,Real), \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG w, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( function<void(DistMatrix<Real>&,Real)>(lossLambda), \
          function<void(DistMatrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Expert versions
     --------------- */ \
  ElError ElModelFitX_ ## SIG \
  ( void (*lossProx)(ElMatrix_ ## SIG,Real), \
    void (*regProx)(ElMatrix_ ## SIG,Real), \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG w, ElModelFitCtrl_ ## SIG ctrl, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( Matrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( Matrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( function<void(Matrix<Real>&,Real)>(lossLambda), \
          function<void(Matrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w), CReflect(ctrl) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElModelFitXDist_ ## SIG \
  ( void (*lossProx)(ElDistMatrix_ ## SIG,Real), \
    void (*regProx)(ElDistMatrix_ ## SIG,Real), \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG w, ElModelFitCtrl_ ## SIG ctrl, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( function<void(DistMatrix<Real>&,Real)>(lossLambda), \
          function<void(DistMatrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w), CReflect(ctrl) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Non-negative matrix factorization
     ================================= */ \
  ElError ElNMF_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElNMFDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y) ) ) } \
  /* Expert versions */ \
  ElError ElNMFX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG X, ElMatrix_ ## SIG Y, \
    ElNMFCtrl_ ## SIG ctrl ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y), \
      CReflect(ctrl) ) ) } \
  ElError ElNMFXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y, \
    ElNMFCtrl_ ## SIG ctrl ) \
  { EL_TRY( NMF( *CReflect(A), *CReflect(X), *CReflect(Y), \
      CReflect(ctrl) ) ) } \
  /* Robust least squares
     ==================== */ \
  ElError ElRLS_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x ) \
  { EL_TRY( RLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  ElError ElRLSDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real rho, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( RLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  ElError ElRLSSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x ) \
  { EL_TRY( RLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  ElError ElRLSDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real rho, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( RLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElRLSX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElRLSXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real rho, ElDistMatrix_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElRLSXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElRLSXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real rho, ElDistMultiVec_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Robust non-negative least squares
     ================================= */ \
  ElError ElRNNLS_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x ) \
  { EL_TRY( RNNLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  ElError ElRNNLSDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real rho, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( RNNLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  ElError ElRNNLSSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x ) \
  { EL_TRY( RNNLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  ElError ElRNNLSDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real rho, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( RNNLS( *CReflect(A), *CReflect(b), rho, *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElRNNLSX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RNNLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElRNNLSXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    Real rho, ElDistMatrix_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RNNLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElRNNLSXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    Real rho, ElMatrix_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RNNLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElRNNLSXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    Real rho, ElDistMultiVec_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( RNNLS( \
      *CReflect(A), *CReflect(b), rho, *CReflect(x), CReflect(ctrl) ) ) } \
  /* Non-negative least squares
     ========================== */ \
  ElError ElNNLS_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElNNLSDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElNNLSSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  ElError ElNNLSDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElDistMultiVec_ ## SIG X ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X) ) ) } \
  /* Expert version 
     -------------- */ \
  ElError ElNNLSX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X, ElNNLSCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  ElError ElNNLSXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG X, ElNNLSCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  ElError ElNNLSXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG X, ElNNLSCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  ElError ElNNLSXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElDistMultiVec_ ## SIG X, ElNNLSCtrl_ ## SIG ctrl ) \
  { EL_TRY( NNLS( *CReflect(A), *CReflect(B), *CReflect(X), \
      CReflect(ctrl) ) ) } \
  /* Long-only portfolio
     =================== */ \
  ElError ElLongOnlyPortfolioSparse_ ## SIG \
  ( ElConstMatrix_ ## SIG d, ElConstSparseMatrix_ ## SIG F, \
    ElConstMatrix_ ## SIG c, Real gamma, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( LongOnlyPortfolio( \
      *CReflect(d), *CReflect(F), \
      *CReflect(c), CReflect(gamma), \
      *CReflect(x) ) ) } \
  ElError ElLongOnlyPortfolioDistSparse_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG d, ElConstDistSparseMatrix_ ## SIG F, \
    ElConstDistMultiVec_ ## SIG c, Real gamma, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( LongOnlyPortfolio( \
      *CReflect(d), *CReflect(F), \
      *CReflect(c), CReflect(gamma), \
      *CReflect(x) ) ) } \
  /* Expert versions
     --------------- */ \
  ElError ElLongOnlyPortfolioXSparse_ ## SIG \
  ( ElConstMatrix_ ## SIG d, ElConstSparseMatrix_ ## SIG F, \
    ElConstMatrix_ ## SIG c, Real gamma, \
    ElMatrix_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LongOnlyPortfolio( \
      *CReflect(d), *CReflect(F), \
      *CReflect(c), CReflect(gamma), \
      *CReflect(x), CReflect(ctrl) ) ) } \
  ElError ElLongOnlyPortfolioXDistSparse_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG d, ElConstDistSparseMatrix_ ## SIG F, \
    ElConstDistMultiVec_ ## SIG c, Real gamma, \
    ElDistMultiVec_ ## SIG x, ElSOCPAffineCtrl_ ## SIG ctrl ) \
  { EL_TRY( LongOnlyPortfolio( \
      *CReflect(d), *CReflect(F), \
      *CReflect(c), CReflect(gamma), \
      *CReflect(x), CReflect(ctrl) ) ) }

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F)

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
