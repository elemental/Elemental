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
  ElError ElLinearProgramFormAugmentedSystem_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG l, ElConstMatrix_ ## SIG s, \
    Real tau, ElSparseMatrix_ ## SIG J, ElMatrix_ ## SIG y ) \
  { EL_TRY( lin_prog::FormAugmentedSystem( \
      *CReflect(A),*CReflect(b),*CReflect(c), \
      *CReflect(x),*CReflect(l),*CReflect(s), \
      tau,*CReflect(J),*CReflect(y)) ) } \
  ElError ElLinearProgramFormAugmentedSystemDist_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG x, ElConstDistMultiVec_ ## SIG l, \
    ElConstDistMultiVec_ ## SIG s, \
    Real tau, ElDistSparseMatrix_ ## SIG J, ElDistMultiVec_ ## SIG y ) \
  { EL_TRY( lin_prog::FormAugmentedSystem( \
      *CReflect(A),*CReflect(b),*CReflect(c), \
      *CReflect(x),*CReflect(l),*CReflect(s), \
      tau,*CReflect(J),*CReflect(y)) ) } \
  ElError ElLinearProgramFormNormalSystem_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG l, ElConstMatrix_ ## SIG s, \
    Real tau, ElSparseMatrix_ ## SIG J, ElMatrix_ ## SIG y ) \
  { EL_TRY( lin_prog::FormNormalSystem( \
      *CReflect(A),*CReflect(b),*CReflect(c), \
      *CReflect(x),*CReflect(l),*CReflect(s), \
      tau,*CReflect(J),*CReflect(y)) ) } \
  ElError ElLinearProgramFormNormalSystemDist_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG x, ElConstDistMultiVec_ ## SIG l, \
    ElConstDistMultiVec_ ## SIG s, \
    Real tau, ElDistSparseMatrix_ ## SIG J, ElDistMultiVec_ ## SIG y ) \
  { EL_TRY( lin_prog::FormNormalSystem( \
      *CReflect(A),*CReflect(b),*CReflect(c), \
      *CReflect(x),*CReflect(l),*CReflect(s), \
      tau,*CReflect(J),*CReflect(y)) ) } \
  ElError ElLinearProgramSolveNormalSystem_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG l, ElConstMatrix_ ## SIG s, \
    Real tau, ElConstSparseMatrix_ ## SIG J, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG dx, ElMatrix_ ## SIG dl, ElMatrix_ ## SIG ds ) \
  { EL_TRY( lin_prog::SolveNormalSystem( \
      *CReflect(A),*CReflect(b),*CReflect(c), \
      *CReflect(x),*CReflect(l),*CReflect(s), \
      tau,*CReflect(J),*CReflect(y), \
      *CReflect(dx), *CReflect(dl), *CReflect(ds)) ) } \
  ElError ElLinearProgramSolveNormalSystemDist_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG x, ElConstDistMultiVec_ ## SIG l, \
    ElConstDistMultiVec_ ## SIG s, \
    Real tau, ElConstDistSparseMatrix_ ## SIG J, \
    ElConstDistMultiVec_ ## SIG y, \
    ElDistMultiVec_ ## SIG dx, ElDistMultiVec_ ## SIG dl, \
    ElDistMultiVec_ ## SIG ds ) \
  { EL_TRY( lin_prog::SolveNormalSystem( \
      *CReflect(A),*CReflect(b),*CReflect(c), \
      *CReflect(x),*CReflect(l),*CReflect(s), \
      tau,*CReflect(J),*CReflect(y), \
      *CReflect(dx), *CReflect(dl), *CReflect(ds)) ) } \
  ElError ElLinearProgramIPFLineSearch_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG l, \
    ElConstMatrix_ ## SIG s, \
    ElConstMatrix_ ## SIG dx, ElConstMatrix_ ## SIG dl, \
    ElConstMatrix_ ## SIG ds, \
    Real gamma, Real beta, Real psi, bool print, Real* alpha ) \
  { EL_TRY( *alpha = lin_prog::IPFLineSearch( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(l), *CReflect(s), \
      *CReflect(dx), *CReflect(dl), *CReflect(ds), \
      gamma, beta, psi, print ) ) } \
  ElError ElLinearProgramIPFLineSearchDist_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElConstDistMultiVec_ ## SIG x, ElConstDistMultiVec_ ## SIG l, \
    ElConstDistMultiVec_ ## SIG s, \
    ElConstDistMultiVec_ ## SIG dx, ElConstDistMultiVec_ ## SIG dl, \
    ElConstDistMultiVec_ ## SIG ds, \
    Real gamma, Real beta, Real psi, bool print, Real* alpha ) \
  { EL_TRY( *alpha = lin_prog::IPFLineSearch( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(l), *CReflect(s), \
      *CReflect(dx), *CReflect(dl), *CReflect(ds), \
      gamma, beta, psi, print ) ) } \
  ElError ElLinearProgramIPF_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG b, ElConstMatrix_ ## SIG c, \
    ElMatrix_ ## SIG x, ElMatrix_ ## SIG l, ElMatrix_ ## SIG s, \
    Real muTol, Real rbTol, Real rcTol, ElInt maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ) \
  { EL_TRY( lin_prog::IPF( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(l), *CReflect(s), \
      muTol, rbTol, rcTol, maxIts, sigma, gamma, beta, psi, print ) ) } \
  ElError ElLinearProgramIPFDist_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, \
    ElConstDistMultiVec_ ## SIG b, ElConstDistMultiVec_ ## SIG c, \
    ElDistMultiVec_ ## SIG x, ElDistMultiVec_ ## SIG l, \
    ElDistMultiVec_ ## SIG s, \
    Real muTol, Real rbTol, Real rcTol, ElInt maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ) \
  { EL_TRY( lin_prog::IPF( \
      *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(x), *CReflect(l), *CReflect(s), \
      muTol, rbTol, rcTol, maxIts, sigma, gamma, beta, psi, print ) ) } \
  ElError ElLinearProgram_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG b, \
    ElConstMatrix_ ## SIG c, ElMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = LinearProgram( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(z) ) ) } \
  ElError ElLinearProgramDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG b, \
    ElConstDistMatrix_ ## SIG c, ElDistMatrix_ ## SIG z, ElInt* numIts ) \
  { EL_TRY( *numIts = LinearProgram( *CReflect(A), *CReflect(b), *CReflect(c), \
      *CReflect(z) ) ) } \
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
  ( ElConstMatrix_ ## SIG P, ElConstMatrix_ ## SIG S, \
    Real lb, Real ub, ElMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = QuadraticProgram( *CReflect(P), *CReflect(S), lb, ub, \
      *CReflect(Z) ) ) } \
  ElError ElQuadraticProgramDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG P, ElConstDistMatrix_ ## SIG S, \
    Real lb, Real ub, ElDistMatrix_ ## SIG Z, ElInt* numIts ) \
  { EL_TRY( *numIts = QuadraticProgram( *CReflect(P), *CReflect(S), lb, ub, \
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
