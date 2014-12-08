/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void QP
( const Matrix<Real>& Q, const Matrix<Real>& A,
  const Matrix<Real>& b, const Matrix<Real>& c, 
        Matrix<Real>& x, const qp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QP"))
    if( ctrl.approach == QP_MEHROTRA )
    {
        // TODO: Use better initializations
        Matrix<Real> y, z; 
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
    else if( ctrl.approach == QP_IPF )
    {
        // TODO: Use better initializations
        Matrix<Real> y, z; 
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::IPF( Q, A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

template<typename Real>
void QP
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, 
        AbstractDistMatrix<Real>& x, const qp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QP"))
    if( ctrl.approach == QP_MEHROTRA )
    {
        // TODO: Use better initializations
        DistMatrix<Real> y(A.Grid()), z(A.Grid());
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
    else if( ctrl.approach == QP_IPF )
    {
        // TODO: Use better initializations
        DistMatrix<Real> y(A.Grid()), z(A.Grid());
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::IPF( Q, A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

template<typename Real>
void QP
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A,
  const Matrix<Real>& b, const Matrix<Real>& c, 
        Matrix<Real>& x, const qp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QP"))
    if( ctrl.approach == QP_MEHROTRA )
    {
        // TODO: Use better initializations
        Matrix<Real> y, z; 
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
    else if( ctrl.approach == QP_IPF )
    {
        // TODO: Use better initializations
        Matrix<Real> y, z; 
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::IPF( Q, A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, 
        DistMultiVec<Real>& x, const qp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QP"))
    if( ctrl.approach == QP_MEHROTRA )
    {
        // TODO: Use better initializations
        DistMultiVec<Real> y(A.Comm()), z(A.Comm()); 
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
    else if( ctrl.approach == QP_IPF )
    {
        // TODO: Use better initializations
        DistMultiVec<Real> y(A.Comm()), z(A.Comm()); 
        Zeros( y, A.Height(), 1 );
        Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
        qp::primal::IPF( Q, A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

#define PROTO(Real) \
  template void QP \
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x, const qp::primal::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x, const qp::primal::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x, const qp::primal::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, const qp::primal::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
