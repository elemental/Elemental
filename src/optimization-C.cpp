/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

ElError ElLPPrimalADMMCtrlDefault_s( ElLPPrimalADMMCtrl_s* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->inv = true;
    ctrl->print = true;
    return EL_SUCCESS;
}

ElError ElLPPrimalADMMCtrlDefault_d( ElLPPrimalADMMCtrl_d* ctrl )
{
    ctrl->rho = 1;
    ctrl->alpha = 1.2;
    ctrl->maxIter = 500;
    ctrl->absTol = 1e-6;
    ctrl->relTol = 1e-4;
    ctrl->inv = true;
    ctrl->print = true;
    return EL_SUCCESS;
}

ElError ElLPPrimalIPFLineSearchCtrlDefault_s
( ElLPPrimalIPFLineSearchCtrl_s* ctrl )
{
    ctrl->gamma = 1e-3;
    ctrl->beta = 2;
    ctrl->psi = 100;
    ctrl->print = false;
}

ElError ElLPPrimalIPFLineSearchCtrlDefault_d
( ElLPPrimalIPFLineSearchCtrl_d* ctrl )
{
    ctrl->gamma = 1e-3;
    ctrl->beta = 2;
    ctrl->psi = 100;
    ctrl->print = false;
}

ElError ElLPPrimalIPFCtrlDefault_s
( ElLPPrimalIPFCtrl_s* ctrl, bool isSparse )
{
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT
                              : EL_LP_PRIMAL_NORMAL_KKT );
    ElLPPrimalIPFLineSearchCtrlDefault_s( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPPrimalIPFCtrlDefault_d
( ElLPPrimalIPFCtrl_d* ctrl, bool isSparse )
{
    ctrl->tol = 1e-8;
    ctrl->maxIts = 1000;
    ctrl->centering = 0.9;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT
                              : EL_LP_PRIMAL_NORMAL_KKT );
    ElLPPrimalIPFLineSearchCtrlDefault_d( &ctrl->lineSearchCtrl );
    ctrl->print = false;
    return EL_SUCCESS;
}

ElError ElLPPrimalMehrotraCtrlDefault_s
( ElLPPrimalMehrotraCtrl_s* ctrl, bool isSparse )
{
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT 
                              : EL_LP_PRIMAL_NORMAL_KKT );;
    return EL_SUCCESS;
}

ElError ElLPPrimalMehrotraCtrlDefault_d
( ElLPPrimalMehrotraCtrl_d* ctrl, bool isSparse )
{
    ctrl->tol = 1e-8;
    ctrl->maxIts = 100;
    ctrl->maxStepRatio = 0.99;
    ctrl->print = false;
    ctrl->system = ( isSparse ? EL_LP_PRIMAL_AUGMENTED_KKT 
                              : EL_LP_PRIMAL_NORMAL_KKT );;
    return EL_SUCCESS;
}

ElError ElLPPrimalCtrlDefault_s( ElLPPrimalCtrl_s* ctrl, bool isSparse )
{
    ctrl->approach = EL_LP_MEHROTRA;
    ElLPPrimalADMMCtrlDefault_s( &ctrl->admmCtrl );
    ElLPPrimalIPFCtrlDefault_s( &ctrl->ipfCtrl, isSparse );
    ElLPPrimalMehrotraCtrlDefault_s( &ctrl->mehrotraCtrl, isSparse );
    ctrl->initialized = false;
    return EL_SUCCESS;
}

ElError ElLPPrimalCtrlDefault_d
( ElLPPrimalCtrl_d* ctrl, bool isSparse )
{
    ctrl->approach = EL_LP_MEHROTRA;
    ElLPPrimalADMMCtrlDefault_d( &ctrl->admmCtrl );
    ElLPPrimalIPFCtrlDefault_d( &ctrl->ipfCtrl, isSparse );
    ElLPPrimalMehrotraCtrlDefault_d( &ctrl->mehrotraCtrl, isSparse );
    ctrl->initialized = false;
    return EL_SUCCESS;
}

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Basis pursuit
     ============= */ \
  ElError ElBasisPursuit_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = BasisPursuit \
      ( *CReflect(A), *CReflect(b), *CReflect(z) ) ) } \
  ElError ElBasisPursuitDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = BasisPursuit \
      ( *CReflect(A), *CReflect(b), *CReflect(z) ) ) } \
  /* Least Absolute Shrinkage and Selection Operator (LASSO)
     ======================================================= */ \
  ElError ElLasso_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, Base<F> lambda, \
    ElMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = Lasso( *CReflect(A), *CReflect(b), lambda, \
      *CReflect(z) ) ) } \
  ElError ElLassoDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, Base<F> lambda, \
    ElDistMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = Lasso( *CReflect(A), *CReflect(b), lambda, \
      *CReflect(z) ) ) } \
  /* Robust Principal Component Analysis
     =================================== */ \
  ElError ElRPCA_ ## SIG \
  ( ElConstMatrix_ ## SIG M, ElMatrix_ ## SIG L, ElMatrix_ ## SIG S ) \
  { EL_TRY( RPCA( *CReflect(M), *CReflect(L), *CReflect(S) ) ) } \
  ElError ElRPCADist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG M, \
    ElDistMatrix_ ## SIG L, ElDistMatrix_ ## SIG S ) \
  { EL_TRY( RPCA( *CReflect(M), *CReflect(L), *CReflect(S) ) ) } \
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
  /* Utilities
     ========= */ \
  /* Coherence
     --------- */ \
  ElError ElCoherence_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* coherence ) \
  { EL_TRY( *coherence = Coherence(*CReflect(A)) ) } \
  ElError ElCoherenceDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* coherence ) \
  { EL_TRY( *coherence = Coherence(*CReflect(A)) ) } \
  /* Covariance
     ---------- */ \
  ElError ElCovariance_ ## SIG \
  ( ElConstMatrix_ ## SIG D, ElMatrix_ ## SIG S ) \
  { EL_TRY( Covariance( *CReflect(D), *CReflect(S) ) ) } \
  ElError ElCovarianceDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG D, ElDistMatrix_ ## SIG S ) \
  { EL_TRY( Covariance( *CReflect(D), *CReflect(S) ) ) } \
  /* Frobenius-norm prox
     ------------------- */ \
  ElError ElFrobeniusProx_ ## SIG ( ElMatrix_ ## SIG A, Base<F> rho ) \
  { EL_TRY( FrobeniusProx( *CReflect(A), rho ) ) } \
  ElError ElFrobeniusProxDist_ ## SIG ( ElDistMatrix_ ## SIG A, Base<F> rho ) \
  { EL_TRY( FrobeniusProx( *CReflect(A), rho ) ) } \
  /* Log barrier
     ----------- */ \
  ElError ElLogBarrier_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, Base<F>* barrier ) \
  { EL_TRY( *barrier = LogBarrier( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElLogBarrierDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, Base<F>* barrier ) \
  { EL_TRY( *barrier = LogBarrier( CReflect(uplo), *CReflect(A) ) ) } \
  /* Log-det divergence
     ------------------ */ \
  ElError ElLogDetDiv_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG B, Base<F>* div ) \
  { EL_TRY( *div = LogDetDiv( CReflect(uplo), *CReflect(A), *CReflect(B) ) ) } \
  ElError ElLogDetDivDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG B, Base<F>* div ) \
  { EL_TRY( *div = LogDetDiv( CReflect(uplo), *CReflect(A), *CReflect(B) ) ) } \
  /* Singular-value soft-thresholding
     -------------------------------- */ \
  ElError ElSVT_ ## SIG \
  ( ElMatrix_ ## SIG A, Base<F> rho, bool relative ) \
  { EL_TRY( SVT( *CReflect(A), rho, relative ) ) } \
  ElError ElSVTDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, Base<F> rho, bool relative ) \
  { EL_TRY( SVT( *CReflect(A), rho, relative ) ) } \
  /* Soft-thresholding
     ----------------- */ \
  ElError ElSoftThreshold_ ## SIG \
  ( ElMatrix_ ## SIG A, Base<F> rho, bool relative ) \
  { EL_TRY( SoftThreshold( *CReflect(A), rho, relative ) ) } \
  ElError ElSoftThresholdDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, Base<F> rho, bool relative ) \
  { EL_TRY( SoftThreshold( *CReflect(A), rho, relative ) ) } \

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* Linear program
     ============== */ \
  /* Primal conic form
     ----------------- */ \
  ElError ElLPPrimal_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElLPPrimalDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElLPPrimalSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElLPPrimalDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  /* Expert version
     ^^^^^^^^^^^^^^ */ \
  ElError ElLPPrimalX_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,      ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z, \
    ElLPPrimalCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), CReflect(ctrl) ) ) } \
  ElError ElLPPrimalXDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x,      ElDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG z, \
    ElLPPrimalCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), CReflect(ctrl) ) ) } \
  ElError ElLPPrimalXSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x,            ElMatrix_ ## SIG y, \
    ElMatrix_ ## SIG z, \
    ElLPPrimalCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), CReflect(ctrl) ) ) } \
  ElError ElLPPrimalXDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x,          ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z, \
    ElLPPrimalCtrl_ ## SIG ctrl ) \
  { EL_TRY( LP( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), CReflect(ctrl) ) ) } \
  /* Dual conic form
     --------------- */ \
  /* TODO */ \
  /* Self-dual conic form
     -------------------- */ \
  /* TODO */ \
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
    ElMatrix_ ## SIG w, Real rho, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( Matrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( Matrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( std::function<void(Matrix<Real>&,Real)>(lossLambda), \
          std::function<void(Matrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w), rho ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElModelFitDist_ ## SIG \
  ( void (*lossProx)(ElDistMatrix_ ## SIG,Real), \
    void (*regProx)(ElDistMatrix_ ## SIG,Real), \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElDistMatrix_ ## SIG w, Real rho, ElInt* numIts ) \
  { try { \
      auto lossLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { lossProx(CReflect(&B),tau); }; \
      auto regLambda = \
        [&]( DistMatrix<Real>& B, Real tau ) { regProx(CReflect(&B),tau); }; \
      *numIts = ModelFit \
        ( std::function<void(DistMatrix<Real>&,Real)>(lossLambda), \
          std::function<void(DistMatrix<Real>&,Real)>(regLambda), \
          *CReflect(A), *CReflect(b), *CReflect(w), rho ); \
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
  /* Non-negative least squares
     ========================== */ \
  ElError ElNonNegativeLeastSquares_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG Y, \
    ElMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = NonNegativeLeastSquares( *CReflect(A), *CReflect(Y), \
      *CReflect(Z) ) ) } \
  ElError ElNonNegativeLeastSquaresDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG Y, \
    ElDistMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = NonNegativeLeastSquares( *CReflect(A), *CReflect(Y), \
      *CReflect(Z) ) ) } \
  /* Quadratic program
     ================= */ \
  /* Primal conic form
     ----------------- */ \
  ElError ElQPPrimal_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( QP( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  ElError ElQPPrimalDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( QP( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  ElError ElQPPrimalSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( QP( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  ElError ElQPPrimalDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( QP( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  /* Infeasible Path-Following
     ^^^^^^^^^^^^^^^^^^^^^^^^^ */ \
  ElError ElQPPrimalIPF_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x, ElMatrix_ ## SIG y, ElMatrix_ ## SIG z ) \
  { EL_TRY( qp::primal::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c),  \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElQPPrimalIPFDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG y, ElDistMatrix_ ## SIG z ) \
  { EL_TRY( \
      qp::primal::IPFCtrl<Real> ctrl; \
      ctrl.print = true; \
      qp::primal::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), ctrl ) ) } \
  ElError ElQPPrimalIPFSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x, ElMatrix_ ## SIG y, ElMatrix_ ## SIG z ) \
  { EL_TRY( qp::primal::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElQPPrimalIPFDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x, ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z ) \
  { EL_TRY( \
      qp::primal::IPFCtrl<Real> ctrl; \
      ctrl.print = true; \
      qp::primal::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), ctrl ) ) } \
  /* Mehrotra Predictor-Corrector
     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */ \
  ElError ElQPPrimalMehrotra_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x, ElMatrix_ ## SIG y, ElMatrix_ ## SIG z ) \
  { EL_TRY( qp::primal::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c),  \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElQPPrimalMehrotraDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG y, ElDistMatrix_ ## SIG z ) \
  { EL_TRY( \
      qp::primal::MehrotraCtrl<Real> ctrl; \
      ctrl.print = true; \
      qp::primal::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), ctrl ) ) } \
  ElError ElQPPrimalMehrotraSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x, ElMatrix_ ## SIG y, ElMatrix_ ## SIG z ) \
  { EL_TRY( qp::primal::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z) ) ) } \
  ElError ElQPPrimalMehrotraDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x, ElDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG z ) \
  { EL_TRY( \
      qp::primal::MehrotraCtrl<Real> ctrl; \
      ctrl.print = true; \
      qp::primal::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(y), *CReflect(z), ctrl ) ) } \
  /* Box form (no linear equalities)
     ------------------------------- */ \
  ElError ElQPBoxADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG C, \
    Real lb, Real ub, ElMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = qp::box::ADMM( *CReflect(Q), *CReflect(C), lb, ub, \
      *CReflect(Z) ) ) } \
  ElError ElQPBoxADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG C, \
    Real lb, Real ub, ElDistMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = qp::box::ADMM( *CReflect(Q), *CReflect(C), lb, ub, \
      *CReflect(Z) ) ) } \
  /* Dual conic form
     --------------- */ \
  /* TODO */ \
  /* Self-dual conic form
     -------------------- */ \
  /* TODO */ \
  /* Support Vector Machine
     ====================== */ \
  ElError ElSVM_ ## SIG \
  ( ElConstMatrix_ ## SIG G, ElConstMatrix_ ## SIG q, \
    ElMatrix_ ## SIG z, Real gamma, ElInt* numIts ) \
  { EL_TRY( *numIts = \
      SVM( *CReflect(G), *CReflect(q), *CReflect(z), gamma ) ) } \
  ElError ElSVMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG G, ElConstDistMatrix_ ## SIG q, \
    ElDistMatrix_ ## SIG z, Real gamma, ElInt* numIts ) \
  { EL_TRY( *numIts = \
      SVM( *CReflect(G), *CReflect(q), *CReflect(z), gamma ) ) } \
  /* Utilities
     ========= */ \
  /* Clip
     ---- */ \
  ElError ElLowerClip_ ## SIG ( ElMatrix_ ## SIG X, Real lowerBound ) \
  { EL_TRY( LowerClip( *CReflect(X), lowerBound ) ) } \
  ElError ElLowerClipDist_ ## SIG ( ElDistMatrix_ ## SIG X, Real lowerBound ) \
  { EL_TRY( LowerClip( *CReflect(X), lowerBound ) ) } \
  ElError ElUpperClip_ ## SIG ( ElMatrix_ ## SIG X, Real upperBound ) \
  { EL_TRY( UpperClip( *CReflect(X), upperBound ) ) } \
  ElError ElUpperClipDist_ ## SIG ( ElDistMatrix_ ## SIG X, Real upperBound ) \
  { EL_TRY( UpperClip( *CReflect(X), upperBound ) ) } \
  ElError ElClip_ ## SIG \
  ( ElMatrix_ ## SIG X, Real lowerBound, Real upperBound ) \
  { EL_TRY( Clip( *CReflect(X), lowerBound, upperBound ) ) } \
  ElError ElClipDist_ ## SIG \
  ( ElDistMatrix_ ## SIG X, Real lowerBound, Real upperBound ) \
  { EL_TRY( Clip( *CReflect(X), lowerBound, upperBound ) ) } \
  /* Hinge-loss prox 
     --------------- */ \
  ElError ElHingeLossProx_ ## SIG ( ElMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( HingeLossProx( *CReflect(A), rho ) ) } \
  ElError ElHingeLossProxDist_ ## SIG ( ElDistMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( HingeLossProx( *CReflect(A), rho ) ) } \
  /* Logistic prox
     ------------- */ \
  ElError ElLogisticProx_ ## SIG ( ElMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( LogisticProx( *CReflect(A), rho ) ) } \
  ElError ElLogisticProxDist_ ## SIG ( ElDistMatrix_ ## SIG A, Real rho ) \
  { EL_TRY( LogisticProx( *CReflect(A), rho ) ) }

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
