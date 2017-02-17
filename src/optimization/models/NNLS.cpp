/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./NNLS/SOCP.hpp"
#include "./NNLS/QP.hpp"
#include "./NNLS/ADMM.hpp"

namespace El {

// SOCP formulation
// ----------------
//
// Solve each problem
//
//   min || A x - b ||_2
//   s.t. x >= 0
//
// via the SOCP
//
//   min t
//   s.t. || A x - b ||_2 <= t, x >= 0.
//
// QP formulation
// --------------
//
// Solve each problem
//
//   min || A x - b ||_2
//   s.t. x >= 0
//
// by transforming it into the explicit QP
//
//   min (1/2) x^T (A^T A) x + (-A^T b)^T x
//   s.t. x >= 0.
//
// Note that the matrix A^T A is cached amongst all instances
// (and this caching is the reason NNLS supports X and B as matrices).
//

template<typename Real>
void NNLS
( const Matrix<Real>& A,
  const Matrix<Real>& B,
        Matrix<Real>& X,
  const NNLSCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == NNLS_SOCP )
        nnls::SOCP( A, B, X, ctrl.socpCtrl );
    else if( ctrl.approach == NNLS_QP )
        nnls::QP( A, B, X, ctrl.qpCtrl );
    else
        nnls::ADMM( A, B, X, ctrl.admmCtrl );
}

template<typename Real>
void NNLS
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& B,
        AbstractDistMatrix<Real>& X,
  const NNLSCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == NNLS_SOCP )
        nnls::SOCP( A, B, X, ctrl.socpCtrl );
    else if( ctrl.approach == NNLS_QP )
        nnls::QP( A, B, X, ctrl.qpCtrl );
    else
        nnls::ADMM( A, B, X, ctrl.admmCtrl );
}

template<typename Real>
void NNLS
( const SparseMatrix<Real>& A,
  const Matrix<Real>& B,
        Matrix<Real>& X,
  const NNLSCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == NNLS_SOCP )
        nnls::SOCP( A, B, X, ctrl.socpCtrl );
    else if( ctrl.approach == NNLS_QP )
        nnls::QP( A, B, X, ctrl.qpCtrl );
    else
        LogicError("ADMM NNLS not yet supported for sparse matrices");
}

template<typename Real>
void NNLS
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& B,
        DistMultiVec<Real>& X,
  const NNLSCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.approach == NNLS_SOCP )
        nnls::SOCP( A, B, X, ctrl.socpCtrl );
    else if( ctrl.approach == NNLS_QP )
        nnls::QP( A, B, X, ctrl.qpCtrl );
    else
        LogicError("ADMM NNLS not yet supported for sparse matrices");
}

#define PROTO(Real) \
  template void NNLS \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& B, \
          Matrix<Real>& X, \
    const NNLSCtrl<Real>& ctrl ); \
  template void NNLS \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& B, \
          AbstractDistMatrix<Real>& X, \
    const NNLSCtrl<Real>& ctrl ); \
  template void NNLS \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& B, \
          Matrix<Real>& X, \
    const NNLSCtrl<Real>& ctrl ); \
  template void NNLS \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& B, \
          DistMultiVec<Real>& X, \
    const NNLSCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
