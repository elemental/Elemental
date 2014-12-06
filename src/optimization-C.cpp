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
  ElError ElLinearProgram_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, ElMatrix_ ## SIG x ) \
  { EL_TRY( LinearProgram( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x) ) ) } \
  ElError ElLinearProgramDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElConstDistMatrix_ ## SIG c, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( LinearProgram( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x) ) ) } \
  ElError ElLinearProgramSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, ElMatrix_ ## SIG x ) \
  { EL_TRY( LinearProgram( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x) ) ) } \
  ElError ElLinearProgramDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistMultiVec_ ## SIG b, \
    ElConstDistMultiVec_ ## SIG c, ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( LinearProgram( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x) ) ) } \
  /* IPF
     --- */ \
  ElError ElLinearProgramIPF_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::IPF( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElLinearProgramIPFDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG s, ElDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::IPFCtrl<Real> ctrl(false); \
      ctrl.print = true; \
      ctrl.lineSearchCtrl.print = true; \
      lin_prog::IPF( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  ElError ElLinearProgramIPFSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::IPF( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElLinearProgramIPFDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG s, ElDistMultiVec_ ## SIG x, \
    ElDistMultiVec_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::IPFCtrl<Real> ctrl(true); \
      ctrl.print = true; \
      lin_prog::IPF( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  /* Mehrotra
     -------- */ \
  ElError ElLinearProgramMehrotra_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::Mehrotra( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElLinearProgramMehrotraDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG s, ElDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::MehrotraCtrl<Real> ctrl(false); \
      ctrl.print = true; \
      lin_prog::Mehrotra( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  ElError ElLinearProgramMehrotraSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::Mehrotra( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElLinearProgramMehrotraDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG s, ElDistMultiVec_ ## SIG x, \
    ElDistMultiVec_ ## SIG l ) \
  { EL_TRY( \
      lin_prog::MehrotraCtrl<Real> ctrl(true); \
      ctrl.print = true; \
      lin_prog::Mehrotra( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  /* ADMM
     ---- */ \
  ElError ElLinearProgramADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, ElMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = lin_prog::ADMM( \
      *CReflect(A), *CReflect(b), *CReflect(c), *CReflect(z) ) ) } \
  ElError ElLinearProgramADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElConstDistMatrix_ ## SIG c, ElDistMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = lin_prog::ADMM( \
      *CReflect(A), *CReflect(b), *CReflect(c), *CReflect(z) ) ) } \
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
  ElError ElQuadraticProgram_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( QuadraticProgram( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  ElError ElQuadraticProgramDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG x ) \
  { EL_TRY( QuadraticProgram( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  ElError ElQuadraticProgramSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( QuadraticProgram( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  ElError ElQuadraticProgramDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x ) \
  { EL_TRY( QuadraticProgram( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), *CReflect(x) ) ) } \
  /* IPF 
     --- */ \
  ElError ElQuadraticProgramIPF_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( quad_prog::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c),  \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElQuadraticProgramIPFDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG s, ElDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG l ) \
  { EL_TRY( \
      quad_prog::IPFCtrl<Real> ctrl; \
      ctrl.print = true; \
      quad_prog::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  ElError ElQuadraticProgramIPFSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( quad_prog::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElQuadraticProgramIPFDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG s, ElDistMultiVec_ ## SIG x, \
    ElDistMultiVec_ ## SIG l ) \
  { EL_TRY( \
      quad_prog::IPFCtrl<Real> ctrl; \
      ctrl.print = true; \
      quad_prog::IPF( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  /* Mehrotra
     -------- */ \
  ElError ElQuadraticProgramMehrotra_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( quad_prog::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c),  \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElQuadraticProgramMehrotraDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG b, ElConstDistMatrix_ ## SIG c, \
    ElDistMatrix_ ## SIG s, ElDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG l ) \
  { EL_TRY( \
      quad_prog::MehrotraCtrl<Real> ctrl; \
      ctrl.print = true; \
      quad_prog::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  ElError ElQuadraticProgramMehrotraSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG Q, ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG s, ElMatrix_ ## SIG x, ElMatrix_ ## SIG l ) \
  { EL_TRY( quad_prog::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l) ) ) } \
  ElError ElQuadraticProgramMehrotraDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG Q, ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG s, ElDistMultiVec_ ## SIG x, \
    ElDistMultiVec_ ## SIG l ) \
  { EL_TRY( \
      quad_prog::MehrotraCtrl<Real> ctrl; \
      ctrl.print = true; \
      quad_prog::Mehrotra( \
      *CReflect(Q), *CReflect(A), \
      *CReflect(b), *CReflect(c), \
      *CReflect(s), *CReflect(x), *CReflect(l), ctrl ) ) } \
  /* ADMM (non-conic form)
     --------------------- */ \
  ElError ElQuadraticProgramADMM_ ## SIG \
  ( ElConstMatrix_ ## SIG Q, ElConstMatrix_ ## SIG C, \
    Real lb, Real ub, ElMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = quad_prog::ADMM( *CReflect(Q), *CReflect(C), lb, ub, \
      *CReflect(Z) ) ) } \
  ElError ElQuadraticProgramADMMDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG Q, ElConstDistMatrix_ ## SIG C, \
    Real lb, Real ub, ElDistMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = quad_prog::ADMM( *CReflect(Q), *CReflect(C), lb, ub, \
      *CReflect(Z) ) ) } \
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
