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
void QuadraticProgram
( const Matrix<Real>& Q, const Matrix<Real>& A,
  const Matrix<Real>& b, const Matrix<Real>& c, 
        Matrix<Real>& x,
  const QuadProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( ctrl.alg == QUAD_PROG_MEHROTRA )
    {
        // TODO: Use better initializations
        Matrix<Real> s, l; 
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::Mehrotra( Q, A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
    else if( ctrl.alg == QUAD_PROG_IPF )
    {
        // TODO: Use better initializations
        Matrix<Real> s, l; 
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::IPF( Q, A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

template<typename Real>
void QuadraticProgram
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, 
        AbstractDistMatrix<Real>& x,
  const QuadProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( ctrl.alg == QUAD_PROG_MEHROTRA )
    {
        // TODO: Use better initializations
        DistMatrix<Real> s(A.Grid()), l(A.Grid());
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::Mehrotra( Q, A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
    else if( ctrl.alg == QUAD_PROG_IPF )
    {
        // TODO: Use better initializations
        DistMatrix<Real> s(A.Grid()), l(A.Grid());
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::IPF( Q, A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

template<typename Real>
void QuadraticProgram
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A,
  const Matrix<Real>& b, const Matrix<Real>& c, 
        Matrix<Real>& x,
  const QuadProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( ctrl.alg == QUAD_PROG_MEHROTRA )
    {
        // TODO: Use better initializations
        Matrix<Real> s, l; 
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::Mehrotra( Q, A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
    else if( ctrl.alg == QUAD_PROG_IPF )
    {
        // TODO: Use better initializations
        Matrix<Real> s, l; 
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::IPF( Q, A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

template<typename Real>
void QuadraticProgram
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, 
        DistMultiVec<Real>& x,
  const QuadProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuadraticProgram"))
    if( ctrl.alg == QUAD_PROG_MEHROTRA )
    {
        // TODO: Use better initializations
        DistMultiVec<Real> s(A.Comm()), l(A.Comm()); 
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::Mehrotra( Q, A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
    else if( ctrl.alg == QUAD_PROG_IPF )
    {
        // TODO: Use better initializations
        DistMultiVec<Real> s(A.Comm()), l(A.Comm()); 
        Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
        Zeros( l, A.Height(), 1 );
        quad_prog::IPF( Q, A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else
        LogicError("ADMM is not yet supported for conic form QPs");
}

#define PROTO(Real) \
  template void QuadraticProgram \
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x, \
    const QuadProgCtrl<Real>& ctrl ); \
  template void QuadraticProgram \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x, \
    const QuadProgCtrl<Real>& ctrl ); \
  template void QuadraticProgram \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x, \
    const QuadProgCtrl<Real>& ctrl ); \
  template void QuadraticProgram \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, \
    const QuadProgCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
