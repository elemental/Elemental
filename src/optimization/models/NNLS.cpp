/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

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
( const Matrix<Real>& A, const Matrix<Real>& B, 
        Matrix<Real>& X, 
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    Matrix<Real> Q, AHat, bHat, c;

    Herk( LOWER, ADJOINT, Real(1), A, Q );
    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );
    Zeros( X,    n, k );

    Matrix<Real> y, z;
    for( Int j=0; j<k; ++j )
    {
        auto x = X( IR(0,n), IR(j,j+1) );
        auto b = B( IR(0,m), IR(j,j+1) );

        Zeros( c, n, 1 );
        Gemv( ADJOINT, Real(-1), A, b, Real(0), c );

        QP( Q, AHat, bHat, c, x, y, z, ctrl );
    }
}

template<typename Real>
void NNLS
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& BPre, 
        AbstractDistMatrix<Real>& XPre,
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    auto APtr = ReadProxy<Real,MC,MR>(&APre);      auto& A = *APtr;
    auto BPtr = ReadProxy<Real,MC,MR>(&BPre);      auto& B = *BPtr;
    auto XPtr = ReadWriteProxy<Real,MC,MR>(&XPre); auto& X = *XPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    const Grid& g = A.Grid();
    DistMatrix<Real> Q(g), AHat(g), bHat(g), c(g);

    Herk( LOWER, ADJOINT, Real(1), A, Q );
    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );
    Zeros( X,    n, k );

    DistMatrix<Real> y(g), z(g);
    for( Int j=0; j<k; ++j )
    {
        auto x = X( IR(0,n), IR(j,j+1) );
        auto b = B( IR(0,m), IR(j,j+1) );

        Zeros( c, n, 1 );
        Gemv( ADJOINT, Real(-1), A, b, Real(0), c );

        QP( Q, AHat, bHat, c, x, y, z, ctrl );
    }
}

template<typename Real>
void NNLS
( const SparseMatrix<Real>& A, const Matrix<Real>& B, 
        Matrix<Real>& X, 
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    SparseMatrix<Real> Q, AHat;
    Matrix<Real> bHat, c;

    Herk( LOWER, ADJOINT, Real(1), A, Q );
    MakeHermitian( LOWER, Q );
    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );
    Zeros( X,    n, k );

    Matrix<Real> y, z;
    for( Int j=0; j<k; ++j )
    {
        auto x = X( IR(0,n), IR(j,j+1) );
        auto b = B( IR(0,m), IR(j,j+1) );

        Zeros( c, n, 1 );
        Multiply( ADJOINT, Real(-1), A, b, Real(0), c );

        QP( Q, AHat, bHat, c, x, y, z, ctrl );
    }
}

template<typename Real>
void NNLS
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& B, 
        DistMultiVec<Real>& X, 
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("NNLS"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> Q(comm), AHat(comm);
    DistMultiVec<Real> bHat(comm), c(comm);

    Herk( LOWER, ADJOINT, Real(1), A, Q );
    MakeHermitian( LOWER, Q );
    Zeros( AHat, 0, n );
    Zeros( bHat, 0, 1 );
    Zeros( X,    n, k );

    DistMultiVec<Real> q(comm), y(comm), z(comm);
    auto& qLoc = q.Matrix();
    auto& XLoc = X.Matrix();
    auto& BLoc = B.LockedMatrix();
    for( Int j=0; j<k; ++j )
    {
        auto xLoc = XLoc( IR(0,XLoc.Height()), IR(j,j+1) );
        auto bLoc = BLoc( IR(0,BLoc.Height()), IR(j,j+1) );

        Zeros( c, n, 1 );
        Zeros( q, m, 1 );
        qLoc = bLoc;
        Multiply( ADJOINT, Real(-1), A, q, Real(0), c );

        Zeros( q, n, 1 );
        qLoc = xLoc;
        QP( Q, AHat, bHat, c, q, y, z, ctrl );
        xLoc = qLoc;
    }
}

#define PROTO(Real) \
  template void NNLS \
  ( const Matrix<Real>& A, const Matrix<Real>& B, \
          Matrix<Real>& X, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void NNLS \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& B, \
          AbstractDistMatrix<Real>& X, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void NNLS \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& B, \
          Matrix<Real>& X, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void NNLS \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& B, \
          DistMultiVec<Real>& X, \
    const qp::direct::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
