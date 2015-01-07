/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_HPP
#define EL_OPTIMIZATION_HPP

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

// Coherence
// =========
template<typename F>
Base<F> Coherence( const Matrix<F>& A );
template<typename F>
Base<F> Coherence( const AbstractDistMatrix<F>& A );

// Least Absolute Shrinkage and Selection Operator (LASSO)
// =======================================================
template<typename F>
Int Lasso
( const Matrix<F>& A, const Matrix<F>& b, Base<F> lambda,
  Matrix<F>& z,
  Base<F> rho=1, Base<F> alpha=1.2, Int maxIter=500, Base<F> absTol=1e-6, 
  Base<F> relTol=1e-4, bool inv=true, bool progress=true );
template<typename F>
Int Lasso
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& b, 
  Base<F> lambda, AbstractDistMatrix<F>& z,
  Base<F> rho=1, Base<F> alpha=1.2, Int maxIter=500, Base<F> absTol=1e-6, 
  Base<F> relTol=1e-4, bool inv=true, bool progress=true );

// Linear program
// ==============

namespace LPApproachNS {
enum LPApproach {
    LP_ADMM,
    LP_IPF,
    LP_IPF_SELFDUAL,     // NOTE: Not yet supported
    LP_MEHROTRA,
    LP_MEHROTRA_SELFDUAL // NOTE: Not yet supported
};
} // namespace LPApproachNS
using namespace LPApproachNS;

namespace lp {

// Infeasible Path-Following Interior Point Method (IPF)
// -----------------------------------------------------
template<typename Real>
struct IPFLineSearchCtrl {
    Real gamma;
    Real beta;
    Real psi;
    Real stepRatio;
    bool print;

    IPFLineSearchCtrl()
    : gamma(1e-3), beta(2), psi(100), stepRatio(1.5), print(false)
    { }
};

namespace direct {

// Attempt to solve a pair of Linear Programs in "direct" conic form:
//
//   min c^T x, 
//   s.t. A x = b, x >= 0
//
//   max -b^T y
//   s.t. A^T y -z + c = 0, z >= 0
//

namespace KKTSystemNS {
enum KKTSystem {
  FULL_KKT,
  AUGMENTED_KKT,
  NORMAL_KKT
};
}
using namespace KKTSystemNS;

template<typename Real>
struct IPFCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real centering; 
    KKTSystem system;

    lp::IPFLineSearchCtrl<Real> lineSearchCtrl;

    bool print;

    IPFCtrl( bool isSparse ) 
    : primalInitialized(false), dualInitialized(false),
      tol(1e-8), maxIts(1000), centering(0.9), print(false)
    { system = ( isSparse ? AUGMENTED_KKT : NORMAL_KKT ); }
};

// Mehrotra's Predictor-Corrector Infeasible Interior Point Method
// ---------------------------------------------------------------
template<typename Real>
struct MehrotraCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real maxStepRatio;
    KKTSystem system;
    bool print;

    // TODO: Add a user-definable (muAff,mu) -> sigma function to replace
    //       the default, (muAff/mu)^3 

    MehrotraCtrl( bool isSparse )
    : primalInitialized(false), dualInitialized(false), 
      tol(1e-8), maxIts(1000), maxStepRatio(0.99), print(false)
    { system = ( isSparse ? AUGMENTED_KKT : NORMAL_KKT ); }
};

// Alternating Direction Method of Multipliers (ADMM)
// --------------------------------------------------
template<typename Real>
struct ADMMCtrl
{
    Real rho;
    Real alpha;
    Int maxIter;
    Real absTol;
    Real relTol;
    bool inv;
    bool print;

    ADMMCtrl()
    : rho(1), alpha(1.2), maxIter(500), absTol(1e-6), relTol(1e-4), inv(true),
      print(true)
    { }
};

// Control structure for the high-level "direct" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach;
    ADMMCtrl<Real> admmCtrl;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl( bool isSparse ) 
    : approach(LP_MEHROTRA), ipfCtrl(isSparse), mehrotraCtrl(isSparse)
    { }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Linear Programs in "affine" conic form:
//
//   min c^T x, 
//   s.t. A x = b, G x + s = h, s >= 0
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z >= 0
//

// Infeasible Path-Following Interior Point Method (IPF)
// -----------------------------------------------------
template<typename Real>
struct IPFCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real centering; 

    lp::IPFLineSearchCtrl<Real> lineSearchCtrl;

    bool print;

    IPFCtrl() 
    : primalInitialized(false), dualInitialized(false),
      tol(1e-8), maxIts(1000), centering(0.9), print(false)
    { }
};

// Mehrotra's Predictor-Corrector Infeasible Interior Point Method
// ---------------------------------------------------------------
template<typename Real>
struct MehrotraCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real maxStepRatio;
    bool print;

    // TODO: Add a user-definable (muAff,mu) -> sigma function to replace
    //       the default, (muAff/mu)^3 

    MehrotraCtrl()
    : primalInitialized(false), dualInitialized(false),
      tol(1e-8), maxIts(1000), maxStepRatio(0.99), print(false)
    { }
};

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl() : approach(LP_MEHROTRA) { }
};

} // namespace affine

} // namespace lp

// Direct conic form
// -----------------
template<typename Real>
void LP
( const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void LP
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void LP
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,       const Matrix<Real>& c,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );
template<typename Real>
void LP
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );

// Affine conic form
// -----------------
template<typename Real>
void LP
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,       Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,             Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,           DistMultiVec<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

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
// For now, we will default to the (DS2) algorithm since it does not require
// explicitly forming A^T A.

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
( std::function<void(Matrix<Real>&,Real)> lossProx,
  std::function<void(Matrix<Real>&,Real)> regProx,
  const Matrix<Real>& A, const Matrix<Real>& b, Matrix<Real>& w,
  Real rho, Int maxIter=1000, bool inv=true, bool progress=true );
template<typename Real>
Int ModelFit
( std::function<void(DistMatrix<Real>&,Real)> lossProx,
  std::function<void(DistMatrix<Real>&,Real)> regProx,
  const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& w,
  Real rho, Int maxIter=1000, bool inv=true, bool progress=true );

// Non-negative matrix factorization
// =================================
// TODO: Generalize to complex
template<typename Real>
void NMF( const Matrix<Real>& A, Matrix<Real>& X, Matrix<Real>& Y );
template<typename Real>
void NMF
( const AbstractDistMatrix<Real>& A, AbstractDistMatrix<Real>& X, 
        AbstractDistMatrix<Real>& Y );

// Non-negative least squares
// ==========================
// TODO: Generalize to complex
template<typename Real>
Int NonNegativeLeastSquares
( const Matrix<Real>& A, const Matrix<Real>& Y, Matrix<Real>& Z,
  Real rho=1., Real alpha=1.2, Int maxIter=500, Real absTol=1e-6,
  Real relTol=1e-4, bool inv=true, bool progress=true );
template<typename Real>
Int NonNegativeLeastSquares
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& Y, 
        AbstractDistMatrix<Real>& Z, 
  Real rho=1., Real alpha=1.2, Int maxIter=500, Real absTol=1e-6,
  Real relTol=1e-4, bool inv=true, bool progress=true );

// Quadratic program
// =================

namespace QPApproachNS {
enum QPApproach {
    QP_ADMM,
    QP_IPF,
    QP_IPF_SELFDUAL,     // NOTE: Not yet supported
    QP_MEHROTRA,
    QP_MEHROTRA_SELFDUAL // NOTE: Not yet supported
};
} // namespace QPApproachNS
using namespace QPApproachNS;

namespace qp {

// Infeasible Path-Following Interior Point Method (IPF)
// -----------------------------------------------------
template<typename Real>
struct IPFLineSearchCtrl {
    Real gamma;
    Real beta;
    Real psi;
    Real stepRatio;
    bool print;

    IPFLineSearchCtrl()
    : gamma(1e-3), beta(2), psi(100), stepRatio(1.5), print(false)
    { }
};

namespace direct {

// Attempt to solve a pair of Quadratic Programs in "direct" conic form:
//
//   min (1/2) x^T Q x + c^T x, 
//   s.t. A x = b, x >= 0
//
//   max (1/2) (A^T y - z + c)^T pinv(Q) (A^T y -z + c) - b^T y
//   s.t. A^T y - z + c in range(Q), z >= 0
//

namespace KKTSystemNS {
enum KKTSystem {
  FULL_KKT,
  AUGMENTED_KKT
};
}
using namespace KKTSystemNS;

template<typename Real>
struct IPFCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real centering; 
    KKTSystem system;

    qp::IPFLineSearchCtrl<Real> lineSearchCtrl;

    bool print;

    IPFCtrl() 
    : primalInitialized(false), dualInitialized(false),
      tol(1e-8), maxIts(1000), centering(0.9), system(AUGMENTED_KKT), 
      print(false)
    { }
};

// Mehrotra's Predictor-Corrector Infeasible Interior Point Method
// ---------------------------------------------------------------
template<typename Real>
struct MehrotraCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real maxStepRatio;
    KKTSystem system;
    bool print;

    // TODO: Add a user-definable (muAff,mu) -> sigma function to replace
    //       the default, (muAff/mu)^3 

    MehrotraCtrl()
    : primalInitialized(false), dualInitialized(false), 
      tol(1e-8), maxIts(1000), maxStepRatio(0.99), system(AUGMENTED_KKT),
      print(false)
    { }
};

// Control structure for the high-level "direct" conic-form QP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    QPApproach approach;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl() : approach(QP_MEHROTRA) { }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Quadratic Programs in "affine" conic form:
//
//   min (1/2) x^T Q x + c^T x, 
//   s.t. A x = b, G x + s = h, s >= 0
//
//   max (1/2) (A^T y + G^T z + c)^T pinv(Q) (A^T y + G^T z + c)  -b^T y - h^T z
//   s.t. A^T y + G^T z + c in range(Q), z >= 0
//

// Infeasible Path-Following Interior Point Method (IPF)
// -----------------------------------------------------
template<typename Real>
struct IPFCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real centering; 

    qp::IPFLineSearchCtrl<Real> lineSearchCtrl;

    bool print;

    IPFCtrl() 
    : primalInitialized(false), dualInitialized(false),
      tol(1e-8), maxIts(1000), centering(0.9), print(false)
    { }
};

// Mehrotra's Predictor-Corrector Infeasible Interior Point Method
// ---------------------------------------------------------------
template<typename Real>
struct MehrotraCtrl {
    bool primalInitialized;
    bool dualInitialized;
    Real tol;
    Int maxIts;
    Real maxStepRatio;
    bool print;

    // TODO: Add a user-definable (muAff,mu) -> sigma function to replace
    //       the default, (muAff/mu)^3 

    MehrotraCtrl()
    : primalInitialized(false), dualInitialized(false),
      tol(1e-8), maxIts(1000), maxStepRatio(0.99), print(false)
    { }
};

// Control structure for the high-level "affine" conic-form QP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    QPApproach approach;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl() : approach(QP_MEHROTRA) { }
};

} // namespace affine

namespace box {

// Solve (a set of) quadratic programs of the form 
//   min (1/2) x' Q x + c' x, subject to l_b <= x <= u_b
//    x

// Alternating Direction Method of Multipliers
// -------------------------------------------
template<typename Real>
struct ADMMCtrl
{
    Real rho;
    Real alpha;
    Int maxIter;
    Real absTol;
    Real relTol;
    bool inv;
    bool print;

    ADMMCtrl()
    : rho(1), alpha(1.2), maxIter(500), absTol(1e-6), relTol(1e-4), inv(true),
      print(true)
    { }
};

template<typename Real>
Int ADMM
( const Matrix<Real>& Q, const Matrix<Real>& C, 
  Real lb, Real ub, Matrix<Real>& Z, 
  const ADMMCtrl<Real>& ctrl=ADMMCtrl<Real>() );
template<typename Real>
Int ADMM
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& C, 
  Real lb, Real ub, AbstractDistMatrix<Real>& Z,
  const ADMMCtrl<Real>& ctrl=ADMMCtrl<Real>() );

} // namespace box

} // namespace qp

// Direct conic form
// -----------------
template<typename Real>
void QP
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,       const Matrix<Real>& c,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );

// Affine conic form
// -----------------
template<typename Real>
void QP
( const Matrix<Real>& Q,
  const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,       Matrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,             Matrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,           DistMultiVec<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

// Basis pursuit denoising (BPDN): 
//   min (1/2) || b - A x ||_2^2 + lambda || x ||_1
// ================================================

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

// Support Vector Machine
// ======================
// NOTE: Currently only soft-margin SVM is supported.

namespace svm {

// TODO: Put this into SVM as an optional backend
template<typename Real>
Int ADMM
( const Matrix<Real>& G, const Matrix<Real>& q, 
        Real gamma,            Matrix<Real>& z,
  Real rho=1, Int maxIter=500, bool inv=true, bool progress=true );
template<typename Real>
Int ADMM
( const AbstractDistMatrix<Real>& G, const AbstractDistMatrix<Real>& q, 
        Real gamma,                        AbstractDistMatrix<Real>& z,
  Real rho=1, Int maxIter=500, bool inv=true, bool progress=true );

} // namespace svm

// NOTE: These routines are currently stubs
template<typename Real>
void SVM
( const Matrix<Real>& G, const Matrix<Real>& q, 
        Real gamma,            Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void SVM
( const AbstractDistMatrix<Real>& G, const AbstractDistMatrix<Real>& q, 
        Real gamma,                        AbstractDistMatrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void SVM
( const SparseMatrix<Real>& G, const Matrix<Real>& q, 
        Real gamma,                  Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void SVM
( const DistSparseMatrix<Real>& G, const DistMultiVec<Real>& q, 
        Real gamma,                      DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

} // namespace El

#include "./optimization/prox.hpp"
#include "./optimization/util.hpp"

#endif // ifndef EL_OPTIMIZATION_HPP
