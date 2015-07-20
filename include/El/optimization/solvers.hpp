/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_SOLVERS_HPP
#define EL_OPTIMIZATION_SOLVERS_HPP

namespace El {

namespace KKTSystemNS {
enum KKTSystem {
  FULL_KKT,
  AUGMENTED_KKT,
  NORMAL_KKT
};
}
using namespace KKTSystemNS;

// Infeasible Path-Following Interior Point Method (IPF)
// =====================================================
template<typename Real>
struct IPFLineSearchCtrl 
{
    Real gamma=1e-3;
    Real beta=2;
    Real psi=100;
    Real stepRatio=1.5;
    bool print=false;
};

template<typename Real>
struct IPFCtrl 
{
    bool primalInit=false, dualInit=false;
    Real minTol=Pow(Epsilon<Real>(),Real(0.3));
    Real targetTol=Pow(Epsilon<Real>(),Real(0.5));
    Int maxIts=1000;
    Real centering=0.9; 
    KKTSystem system=FULL_KKT;
    IPFLineSearchCtrl<Real> lineSearchCtrl;

    RegQSDCtrl<Real> qsdCtrl;
    bool outerEquil=true, innerEquil=true;
    Int basisSize = 6;
    bool print=false;
    bool time=false;
};

// Mehrotra's Predictor-Corrector Infeasible Interior Point Method
// ===============================================================
template<typename Real>
inline Real StepLengthCentrality
( Real mu, Real muAff, Real alphaAffPri, Real alphaAffDual )
{ return Pow(1-Min(alphaAffPri,alphaAffDual),Real(3)); }

template<typename Real>
inline Real MehrotraCentrality
( Real mu, Real muAff, Real alphaAffPri, Real alphaAffDual )
{ return Min(Pow(muAff/mu,Real(3)),Real(1)); }

template<typename Real>
struct MehrotraCtrl 
{
    bool primalInit=false, dualInit=false;
    Real minTol=Pow(Epsilon<Real>(),Real(0.3));
    Real targetTol=Pow(Epsilon<Real>(),Real(0.5));
    Int maxIts=1000;
    Real maxStepRatio=0.99;
    KKTSystem system=FULL_KKT;

    RegQSDCtrl<Real> qsdCtrl;
    bool outerEquil=true, innerEquil=true;
    Int basisSize = 6;
    bool print=false;
    bool time=false;

    // TODO: Add a user-definable (muAff,mu) -> sigma function to replace
    //       the default, (muAff/mu)^3 
};

// Alternating Direction Method of Multipliers
// ===========================================
template<typename Real>
struct ADMMCtrl
{
    Real rho=1;
    Real alpha=1.2;
    Int maxIter=500;
    Real absTol=1e-6;
    Real relTol=1e-4;
    bool inv=true;
    bool print=true;
};

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

namespace direct {

// Attempt to solve a pair of Linear Programs in "direct" conic form:
//
//   min c^T x, 
//   s.t. A x = b, x >= 0
//
//   max -b^T y
//   s.t. A^T y -z + c = 0, z >= 0
//

// Control structure for the high-level "direct" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach=LP_MEHROTRA;
    ADMMCtrl<Real> admmCtrl;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl( bool isSparse ) 
    { 
        ipfCtrl.system = ( isSparse ? AUGMENTED_KKT : NORMAL_KKT );
        mehrotraCtrl.system = ( isSparse ? AUGMENTED_KKT : NORMAL_KKT );
    }
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

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach=LP_MEHROTRA;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;
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

namespace direct {

// Attempt to solve a pair of Quadratic Programs in "direct" conic form:
//
//   min (1/2) x^T Q x + c^T x, 
//   s.t. A x = b, x >= 0
//
//   max (1/2) (A^T y - z + c)^T pinv(Q) (A^T y -z + c) - b^T y
//   s.t. A^T y - z + c in range(Q), z >= 0
//

// Control structure for the high-level "direct" conic-form QP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    QPApproach approach=QP_MEHROTRA;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl()
    {
        ipfCtrl.system = AUGMENTED_KKT;
        mehrotraCtrl.system = AUGMENTED_KKT;
    }
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

// Control structure for the high-level "affine" conic-form QP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    QPApproach approach=QP_MEHROTRA;
    IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;
};

} // namespace affine

namespace box {

// Solve (a set of) quadratic programs of the form 
//   min (1/2) x' Q x + c' x, subject to l_b <= x <= u_b
//    x

// Alternating Direction Method of Multipliers
// -------------------------------------------
template<typename Real>
Int ADMM
( const Matrix<Real>& Q, const Matrix<Real>& C, 
  Real lb, Real ub, Matrix<Real>& X, 
  const ADMMCtrl<Real>& ctrl=ADMMCtrl<Real>() );
template<typename Real>
Int ADMM
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& C, 
  Real lb, Real ub, AbstractDistMatrix<Real>& X,
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

// Second-order Cone Program
// =========================
namespace SOCPApproachNS {
enum SOCPApproach {
  SOCP_ADMM,             // NOTE: Not yet supported
  SOCP_IPF,              // NOTE: Not yet supported
  SOCP_IPF_SELFDUAL,     // NOTE: Not yet supported
  SOCP_MEHROTRA,
  SOCP_MEHROTRA_SELFDUAL // NOTE: Not yet supported
};
} // namespace SOCPApproachNS
using namespace SOCPApproachNS;

namespace socp {
namespace direct {

// Attempt to solve a pair of Second-Order Cone Programs in "direct" conic form:
//
//   min c^T x, 
//   s.t. A x = b, x in K,
//
//   max -b^T y
//   s.t. A^T y - z + c = 0, z in K,
//
// where the cone K is a product of second-order cones.
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    SOCPApproach approach=SOCP_MEHROTRA;
    //IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl()
    {
        mehrotraCtrl.system = AUGMENTED_KKT;
        mehrotraCtrl.minTol = Pow(Epsilon<Real>(),Real(0.25));
        mehrotraCtrl.targetTol = Pow(Epsilon<Real>(),Real(0.5));
    }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Second-Order Cone Programs in "affine" conic form:
//
//   min c^T x, 
//   s.t. A x = b, G x + s = h, s in K,
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z in K,
//
// where the cone K is a product of second-order cones.
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    SOCPApproach approach=SOCP_MEHROTRA;
    //IPFCtrl<Real> ipfCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl()
    {
        mehrotraCtrl.minTol = Pow(Epsilon<Real>(),Real(0.25));
        mehrotraCtrl.targetTol = Pow(Epsilon<Real>(),Real(0.5));
    }
};

} // namespace affine
} // namespace socp

// Direct conic form
// -----------------
template<typename Real>
void SOCP
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,  
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, 
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x,       
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,       
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z, 
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );

// Affine conic form
// -----------------
template<typename Real>
void SOCP
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,  
        Matrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, 
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x,       
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,       
        AbstractDistMatrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const SparseMatrix<Real>& A, 
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z, 
        DistMultiVec<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_SOLVERS_HPP
