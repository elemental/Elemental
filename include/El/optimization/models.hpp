/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_MODELS_HPP
#define EL_OPTIMIZATION_MODELS_HPP

namespace El {

namespace RegularizationNS {
enum Regularization {
  NO_PENALTY,
  L1_PENALTY,
  L2_PENALTY
};
}
using namespace RegularizationNS;

// TODO: Modify the following routines to use control structures instead

// Basis pursuit: min || x ||_1 such that A x = b
// ==============================================

namespace bp {

// Put this into BP as an optional backend
template<typename F>
Int ADMM
( const Matrix<F>& A, const Matrix<F>& b,
  Matrix<F>& z,
  Base<F> rho=1., Base<F> alpha=1.2, Int maxIter=500, Base<F> absTol=1e-6, 
  Base<F> relTol=1e-4, bool usePinv=false, Base<F> pinvTol=0, 
  bool progress=true );
template<typename F>
Int ADMM
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& b,
        AbstractDistMatrix<F>& z,
  Base<F> rho=1., Base<F> alpha=1.2, Int maxIter=500, Base<F> absTol=1e-6, 
  Base<F> relTol=1e-4, bool usePinv=false, Base<F> pinvTol=0,
  bool progress=true );

} // namespace bp

template<typename Real>
void BP
( const Matrix<Real>& A, const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void BP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
        AbstractDistMatrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void BP
( const SparseMatrix<Real>& A, const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );
template<typename Real>
void BP
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b,
        DistMultiVec<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );

// Chebyshev point: min || A x - b||_oo
// ====================================
template<typename Real>
void CP
( const Matrix<Real>& A, const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void CP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void CP
( const SparseMatrix<Real>& A, const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void CP
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

// Least Absolute Value: min || A x - b||_1
// ========================================
template<typename Real>
void LAV
( const Matrix<Real>& A, const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LAV
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LAV
( const SparseMatrix<Real>& A, const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LAV
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

// Dantzig selector: min || x ||_1 s.t. || A^T (b - A x) ||_oo <= lambda
// =====================================================================
// TODO: 
//  Add the ability to switch between the (DS1) and (DS2) affine LP
//  formulations described in 
//
//  Michael Friedlander and Michael Saunders,
//  "Discussion: The Dantzig Selector: Statistical estimation when p is much
//   larger than n",
//  The Annals of Statistics, Vol. 35, No. 6, pp. 2385--2391, 2007.
//
// Elemental currently defaults to (DS1) for dense matrices and (DS2) for
// sparse matrices.
//

template<typename Real>
void DS
( const Matrix<Real>& A, const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void DS
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void DS
( const SparseMatrix<Real>& A, const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void DS
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b,
        Real lambda,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

// Logistic Regression
// ===================
template<typename Real>
Int LogisticRegression
( const Matrix<Real>& G, const Matrix<Real>& q, Matrix<Real>& z,
  Real gamma, Regularization penalty=L1_PENALTY,
  Real rho=1, Int maxIter=500, bool inv=true, bool progress=true );
template<typename Real>
Int LogisticRegression
( const AbstractDistMatrix<Real>& G, const AbstractDistMatrix<Real>& q, 
        AbstractDistMatrix<Real>& z,
  Real gamma, Regularization penalty=L1_PENALTY,
  Real rho=1, Int maxIter=500, bool inv=true, bool progress=true );

// Fit a model with using a loss function plus regularization
// ==========================================================
// TODO: Implement these functions
template<typename Real>
Int ModelFit
( function<void(Matrix<Real>&,Real)> lossProx,
  function<void(Matrix<Real>&,Real)> regProx,
  const Matrix<Real>& A, const Matrix<Real>& b, Matrix<Real>& w,
  Real rho, Int maxIter=1000, bool inv=true, bool progress=true );
template<typename Real>
Int ModelFit
( function<void(DistMatrix<Real>&,Real)> lossProx,
  function<void(DistMatrix<Real>&,Real)> regProx,
  const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& w,
  Real rho, Int maxIter=1000, bool inv=true, bool progress=true );

// Non-negative matrix factorization
// =================================
// TODO: Generalize to complex
template<typename Real>
void NMF
( const Matrix<Real>& A, 
        Matrix<Real>& X, Matrix<Real>& Y,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void NMF
( const AbstractDistMatrix<Real>& A, 
        AbstractDistMatrix<Real>& X, AbstractDistMatrix<Real>& Y,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
// TODO: Sparse versions

// Non-negative least squares
// ==========================
// TODO: Generalize to complex
// NOTE: The following can solve a *sequence* of NNLS problems

namespace nnls {

template<typename Real>
Int ADMM
( const Matrix<Real>& A, const Matrix<Real>& B, 
        Matrix<Real>& X,
  const qp::box::ADMMCtrl<Real>& ctrl=qp::box::ADMMCtrl<Real>() );
template<typename Real>
Int ADMM
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& B, 
        AbstractDistMatrix<Real>& X, 
  const qp::box::ADMMCtrl<Real>& ctrl=qp::box::ADMMCtrl<Real>() );

} // namespace nnls

template<typename Real>
void NNLS
( const Matrix<Real>& A, const Matrix<Real>& B, 
        Matrix<Real>& X,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void NNLS
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& B, 
        AbstractDistMatrix<Real>& X, 
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void NNLS
( const SparseMatrix<Real>& A, const Matrix<Real>& B, 
        Matrix<Real>& X,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void NNLS
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& B, 
        DistMultiVec<Real>& X,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );

// Basis pursuit denoising (BPDN), a.k.a.,
// Least absolute selection and shrinkage operator (Lasso):
//   min (1/2) || b - A x ||_2^2 + lambda || x ||_1
// ================================================

namespace bpdn {

template<typename F>
Int ADMM
( const Matrix<F>& A, const Matrix<F>& b, Base<F> lambda,
  Matrix<F>& z,
  Base<F> rho=1, Base<F> alpha=1.2, Int maxIter=500, Base<F> absTol=1e-6, 
  Base<F> relTol=1e-4, bool inv=true, bool progress=true );
template<typename F>
Int ADMM
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& b, 
  Base<F> lambda, AbstractDistMatrix<F>& z,
  Base<F> rho=1, Base<F> alpha=1.2, Int maxIter=500, Base<F> absTol=1e-6, 
  Base<F> relTol=1e-4, bool inv=true, bool progress=true );

} // namespace bpdn

template<typename Real>
void BPDN
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void BPDN
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void BPDN
( const SparseMatrix<Real>& A, const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void BPDN
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b,
        Real lambda,
        DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

// Elastic net (EN): 
//   min || b - A x ||_2^2 + lambda_1 || x ||_1 + lambda_2 || x ||_2^2
// ===================================================================

template<typename Real>
void EN
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Real lambda1,          Real lambda2,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void EN
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
        Real lambda1,                      Real lambda2,
        AbstractDistMatrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void EN
( const SparseMatrix<Real>& A, const Matrix<Real>& b,
        Real lambda1,                Real lambda2,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void EN
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b,
        Real lambda1,                    Real lambda2,
        DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

// Robust Principal Component Analysis (RPCA)
// ==========================================

template<typename Real>
struct RPCACtrl
{
    bool useALM;
    bool usePivQR;
    bool progress;

    Int numPivSteps;
    Int maxIts;

    Real tau;
    Real beta;
    Real rho;
    Real tol;

    RPCACtrl() 
    : useALM(true), usePivQR(false), progress(true), 
      numPivSteps(75), maxIts(1000),
      tau(0), beta(1), rho(6), tol(1e-5)
    { }
};

template<typename F>
void RPCA
( const Matrix<F>& M, Matrix<F>& L, Matrix<F>& S,
  const RPCACtrl<Base<F>>& ctrl=RPCACtrl<Base<F>>() );

template<typename F>
void RPCA
( const AbstractDistMatrix<F>& M, AbstractDistMatrix<F>& L, 
        AbstractDistMatrix<F>& S,
  const RPCACtrl<Base<F>>& ctrl=RPCACtrl<Base<F>>() );

// Sparse inverse covariance selection
// ===================================
template<typename F>
Int SparseInvCov
( const Matrix<F>& D, Base<F> lambda, Matrix<F>& Z,
  Base<F> rho=1., Base<F> alpha=1.2, Int maxIter=500,
  Base<F> absTol=1e-6, Base<F> relTol=1e-4, bool progress=true );
template<typename F>
Int SparseInvCov
( const AbstractDistMatrix<F>& D, Base<F> lambda, AbstractDistMatrix<F>& Z,
  Base<F> rho=1., Base<F> alpha=1.2, Int maxIter=500,
  Base<F> absTol=1e-6, Base<F> relTol=1e-4, bool progress=true );

// Support Vector Machine (soft-margin)
// ====================================
// TODO: Use the formulation described in 
//       <http://abel.ee.ucla.edu/cvxopt/applications/svm/>
//
// min_{w,beta,z} (1/2) || w ||_2^2 + lambda 1^T z
//
// s.t. | -diag(d) A, -d, -I | | w    | <= | -1 |
//      |       0,     0, -I | | beta |    |  0 |
//                             | z    |
//
// The output, x, is set to the concatenation of w and beta, x := [w; beta].
//

namespace svm {

// TODO: Put this into SVM as an optional backend
template<typename Real>
Int ADMM
( const Matrix<Real>& A, const Matrix<Real>& d, 
        Real lambda,           Matrix<Real>& x,
  Real rho=1, Int maxIter=500, bool inv=true, bool progress=true );
template<typename Real>
Int ADMM
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& d, 
        Real lambda,                       AbstractDistMatrix<Real>& x,
  Real rho=1, Int maxIter=500, bool inv=true, bool progress=true );

} // namespace svm

// TODO: Switch to explicitly returning w, beta, and z, as it is difficult
//       for users to unpack a DistMultiVec
template<typename Real>
void SVM
( const Matrix<Real>& A, const Matrix<Real>& d, 
        Real lambda,           Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void SVM
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& d, 
        Real lambda,                       AbstractDistMatrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void SVM
( const SparseMatrix<Real>& A, const Matrix<Real>& d, 
        Real lambda,                 Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void SVM
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& d, 
        Real lambda,                     DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

// 1D total variation denoising (TV):
//
//   min (1/2) || b - x ||_2^2 + lambda || D x ||_1,
//
// where D is the 1D finite-difference operator.
// =================================================
// We follow the formulation used within CVXOPT:
//
//   min (1/2) || b - x ||_2^2 + lambda 1^T y
//   s.t. -y <= D x <= y,
//
// where x is in R^n and y is in R^(n-1).
//

template<typename Real>
void TV
( const AbstractDistMatrix<Real>& b,
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void TV
( const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void TV
( const DistMultiVec<Real>& b,
        Real lambda,
        DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_MODELS_HPP
