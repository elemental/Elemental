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
void LP
( const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c, 
        Matrix<Real>& x, const lp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_ADMM )
    {
       lp::primal::ADMM( A, b, c, x, ctrl.admmCtrl );
    }
    else if( ctrl.approach == LP_IPF )
    {
       // TODO: Use better initializations
       Matrix<Real> y, z;
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else if( ctrl.approach == LP_MEHROTRA )
    {
       // TODO: Use better initializations
       Matrix<Real> y, z;
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
}

template<typename Real>
void LP
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x, const lp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_ADMM )
    {
       lp::primal::ADMM( A, b, c, x, ctrl.admmCtrl );
    }
    else if( ctrl.approach == LP_IPF )
    {
       // TODO: Use better initializations
       DistMatrix<Real> y(A.Grid()), z(A.Grid());
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else if( ctrl.approach == LP_MEHROTRA )
    {
       // TODO: Use better initializations
       DistMatrix<Real> y(A.Grid()), z(A.Grid());
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
}

template<typename Real>
void LP
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c, 
        Matrix<Real>& x, const lp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
    {
       // TODO: Use better initializations
       Matrix<Real> y, z;
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else // ctrl.approach == LP_MEHROTRA
    {
       // TODO: Use better initializations
       Matrix<Real> y, z;
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
}

template<typename Real>
void LP
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, 
        DistMultiVec<Real>& x, const lp::primal::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
    {
       // TODO: Use better initializations
       DistMultiVec<Real> y(A.Comm()), z(A.Comm());
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    }
    else // ctrl.approach == LP_MEHROTRA
    {
       // TODO: Use better initializations
       DistMultiVec<Real> y(A.Comm()), z(A.Comm());
       Zeros( y, A.Height(), 1 );
       Uniform( z, A.Width(), 1, Real(0.5), Real(0.49) );
       lp::primal::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    }
}

#define PROTO(Real) \
  template void LP \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x, const lp::primal::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x, const lp::primal::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x, const lp::primal::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, const lp::primal::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
