/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Solve the problem 
//
//   min || A x - b ||_2 
//   s.t. x >= 0
//
// by transforming it into the explicit QP
//
//   min (1/2) x^T (A^T A) x + (-A^T b)^T x
//   s.t. x >= 0
//

template<typename Real>
void NNLS
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x, 
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    const Int n = A.Width();
    Matrix<Real> Q, AHat, bHat, c;

    Herk( LOWER, ADJOINT, Real(1), A, Q );

    Zeros( c, n, 1 );
    Gemv( ADJOINT, Real(-1), A, b, Real(0), c );

    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );

    Matrix<Real> y, z;
    QP( Q, AHat, bHat, c, x, y, z, ctrl );
}

template<typename Real>
void NNLS
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    auto APtr = ReadProxy<Real,MC,MR>( &APre );
    auto& A = *APtr;

    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<Real> Q(g), AHat(g), bHat(g), c(g);

    Herk( LOWER, ADJOINT, Real(1), A, Q );

    Zeros( c, n, 1 );
    Gemv( ADJOINT, Real(-1), A, b, Real(0), c );

    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );
    
    DistMatrix<Real> y(g), z(g);
    QP( Q, AHat, bHat, c, x, y, z, ctrl );
}

template<typename Real>
void NNLS
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x, 
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    const Int n = A.Width();
    SparseMatrix<Real> Q, AHat;
    Matrix<Real> bHat, c;

    Herk( LOWER, ADJOINT, Real(1), A, Q );
    MakeHermitian( LOWER, Q );

    Zeros( c, n, 1 );
    Multiply( ADJOINT, Real(-1), A, b, Real(0), c );

    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );

    Matrix<Real> y, z;
    QP( Q, AHat, bHat, c, x, y, z, ctrl );
}

template<typename Real>
void NNLS
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x, 
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> Q(comm), AHat(comm);
    DistMultiVec<Real> bHat(comm), c(comm);

    Herk( LOWER, ADJOINT, Real(1), A, Q );
    MakeHermitian( LOWER, Q );

    Zeros( c, n, 1 );
    Multiply( ADJOINT, Real(-1), A, b, Real(0), c );

    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );

    DistMultiVec<Real> y(comm), z(comm);
    QP( Q, AHat, bHat, c, x, y, z, ctrl );
}

#define PROTO(Real) \
  template void NNLS \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void NNLS \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void NNLS \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void NNLS \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          DistMultiVec<Real>& x, \
    const qp::direct::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
