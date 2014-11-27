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
void LinearProgram
( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, 
  Matrix<Real>& x, const LinProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))
    if( ctrl.alg == LIN_PROG_ADMM )
    {
       lin_prog::ADMM( A, b, c, x, ctrl.admmCtrl );
    }
    else if( ctrl.alg == LIN_PROG_IPF )
    {
       Matrix<Real> s, l;
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::IPF( A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else if( ctrl.alg == LIN_PROG_MEHROTRA )
    {
       Matrix<Real> s, l;
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::Mehrotra( A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
}

template<typename Real>
void LinearProgram
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,       AbstractDistMatrix<Real>& x,
  const LinProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))
    if( ctrl.alg == LIN_PROG_ADMM )
    {
       lin_prog::ADMM( A, b, c, x, ctrl.admmCtrl );
    }
    else if( ctrl.alg == LIN_PROG_IPF )
    {
       DistMatrix<Real> s(A.Grid()), l(A.Grid());
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::IPF( A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else if( ctrl.alg == LIN_PROG_MEHROTRA )
    {
       DistMatrix<Real> s(A.Grid()), l(A.Grid());
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::Mehrotra( A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
}

template<typename Real>
void LinearProgram
( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, 
  Matrix<Real>& x, const LinProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))
    if( ctrl.alg == LIN_PROG_IPF )
    {
       Matrix<Real> s, l;
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::IPF( A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else // ctrl.alg == LIN_PROG_MEHROTRA
    {
       Matrix<Real> s, l;
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::Mehrotra( A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
}

template<typename Real>
void LinearProgram
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, 
  DistMultiVec<Real>& x, const LinProgCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))
    if( ctrl.alg == LIN_PROG_IPF )
    {
       DistMultiVec<Real> s(A.Comm()), l(A.Comm());
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::IPF( A, b, c, s, x, l, ctrl.ipfCtrl );
    }
    else // ctrl.alg == LIN_PROG_MEHROTRA
    {
       DistMultiVec<Real> s(A.Comm()), l(A.Comm());
       Uniform( s, A.Width(), 1, Real(0.5), Real(0.49) );
       Zeros( l, A.Height(), 1 );
       lin_prog::Mehrotra( A, b, c, s, x, l, ctrl.mehrotraCtrl );
    }
}

#define PROTO(Real) \
  template void LinearProgram \
  ( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x, const LinProgCtrl<Real>& ctrl ); \
  template void LinearProgram \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c,       AbstractDistMatrix<Real>& x, \
    const LinProgCtrl<Real>& ctrl ); \
  template void LinearProgram \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x, const LinProgCtrl<Real>& ctrl ); \
  template void LinearProgram \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& x, const LinProgCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
