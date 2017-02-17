/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_SOLVERS_QP_HPP
#define EL_OPTIMIZATION_SOLVERS_QP_HPP

#include <El/optimization/solvers/util.hpp>

namespace El {

namespace QPApproachNS {
enum QPApproach {
  QP_ADMM,
  QP_MEHROTRA
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
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl() { mehrotraCtrl.system = AUGMENTED_KKT; }
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
    MehrotraCtrl<Real> mehrotraCtrl;
};

} // namespace affine

namespace box {

// Solve (a set of) quadratic programs of the form
//   min (1/2) x' Q x + c' x, subject to l_b <= x <= u_b
//    x

// Alternating Direction Method of Multipliers
// -------------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Int ADMM
( const Matrix<Real>& Q,
  const Matrix<Real>& C,
        Real lb,
        Real ub,
        Matrix<Real>& X,
  const ADMMCtrl<Real>& ctrl=ADMMCtrl<Real>() );
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Int ADMM
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& C,
        Real lb,
        Real ub,
        AbstractDistMatrix<Real>& X,
  const ADMMCtrl<Real>& ctrl=ADMMCtrl<Real>() );

} // namespace box

} // namespace qp

// Direct conic form
// -----------------
template<typename Real>
void QP
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );

// Affine conic form
// -----------------
template<typename Real>
void QP
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_SOLVERS_QP_HPP
