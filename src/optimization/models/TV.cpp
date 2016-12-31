/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// 1D total variation denoising (TV):
//
//   min (1/2) || b - x ||_2^2 + lambda || D x ||_1,
//
// where D is the 1D finite-difference operator. We follow the formulation used
// within CVXOPT:
//
//   min (1/2) || b - x ||_2^2 + lambda 1^T t
//   s.t. -t <= D x <= t,
//
// where x is in R^n and t is in R^(n-1).
//

namespace El {

template<typename Real>
void TV
( const AbstractDistMatrix<Real>& b,
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    DistMultiVec<Real> bDMV(b.Grid()), xDMV(b.Grid());
    bDMV = b;
    xDMV = x;
    TV( bDMV, lambda, xDMV, ctrl );
    Copy( xDMV, x );
}

template<typename Real>
void TV
( const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = b.Height();
    const Range<Int> xInd(0,n), tInd(n,2*n-1);

    SparseMatrix<Real> Q, A, G;
    Matrix<Real> c, bHat, h;

    // Q := | I 0 |
    //      | 0 0 |
    // ============
    Zeros( Q, 2*n-1, 2*n-1 );
    Q.Reserve( n );
    for( Int e=0; e<n; ++e )
        Q.QueueUpdate( e, e, Real(1) );
    Q.ProcessQueues();

    // c := [-b;lambda]
    // =================
    Zeros( c, 2*n-1, 1 );
    auto cx = c( xInd, ALL );
    cx = b; cx *= -1;
    auto ct = c( tInd, ALL );
    Fill( ct, lambda );

    // A := []
    // =======
    Zeros( A, 0, 2*n-1 );

    // bHat := []
    // ==========
    Zeros( bHat, 0, 1 );

    // G := |  D -I |
    //      | -D -I |
    // ==============
    Zeros( G, 2*(n-1), 2*n-1 );
    G.Reserve( 6*(n-1) );
    for( Int e=0; e<n-1; ++e )
    {
        // Queue D
        G.QueueUpdate( e, e,   Real(1) );
        G.QueueUpdate( e, e+1, Real(-1) );
        // Queue -D
        G.QueueUpdate( e+n-1, e,   Real(-1) );
        G.QueueUpdate( e+n-1, e+1, Real(1) );
        // Queue the -I's
        G.QueueUpdate( e,     e+n, Real(-1) );
        G.QueueUpdate( e+n-1, e+n, Real(-1) );
    }
    G.ProcessQueues();

    // h := 0
    // ======
    Zeros( h, 2*(n-1), 1 );

    // Solve the affine QP
    // ===================
    Matrix<Real> xHat, y, z, s;
    QP( Q, A, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    x = xHat( xInd, ALL );
}

template<typename Real>
void TV
( const DistMultiVec<Real>& b,
        Real lambda,
        DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = b.Height();
    const Grid& grid = b.Grid();

    DistSparseMatrix<Real> Q(grid), A(grid), G(grid);
    DistMultiVec<Real> c(grid), bHat(grid), h(grid);

    auto& bLoc = b.LockedMatrix();
    auto& cLoc = c.Matrix();

    // Q := | I 0 |
    //      | 0 0 |
    // ============
    Zeros( Q, 2*n-1, 2*n-1 );
    Int numIdentUpdates = 0;
    for( Int iLoc=0; iLoc<Q.LocalHeight(); ++iLoc )
        if( Q.GlobalRow(iLoc) < n )
            ++numIdentUpdates;
        else
            break;
    Q.Reserve( numIdentUpdates );
    for( Int iLoc=0; iLoc<Q.LocalHeight(); ++iLoc )
    {
        const Int i = Q.GlobalRow(iLoc);
        if( i < n )
            Q.QueueLocalUpdate( iLoc, i, Real(1) );
        else
            break;
    }
    Q.ProcessLocalQueues();

    // c := [-b;lambda]
    // =================
    Zeros( c, 2*n-1, 1 );
    c.Reserve( b.LocalHeight() );
    for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
        c.QueueUpdate( b.GlobalRow(iLoc), 0, -bLoc(iLoc) );
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
        if( c.GlobalRow(iLoc) > n )
            cLoc(iLoc) = lambda;
    c.ProcessQueues();

    // A := []
    // =======
    Zeros( A, 0, 2*n-1 );

    // bHat := []
    // ==========
    Zeros( bHat, 0, 1 );

    // G := |  D -I |
    //      | -D -I |
    // ==============
    Zeros( G, 2*(n-1), 2*n-1 );
    G.Reserve( 3*G.LocalHeight() );
    for( Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
    {
        const Int i = G.GlobalRow(iLoc);
        if( i < n-1 )
        {
            G.QueueLocalUpdate( iLoc, i,   Real( 1) );
            G.QueueLocalUpdate( iLoc, i+1, Real(-1) );
            G.QueueLocalUpdate( iLoc, i+n, Real(-1) );
        }
        else
        {
            G.QueueLocalUpdate( iLoc, i-(n-1),   Real(-1) );
            G.QueueLocalUpdate( iLoc, i+1-(n-1), Real( 1) );
            G.QueueLocalUpdate( iLoc, i+n-(n-1), Real(-1) );
        }
    }
    G.ProcessLocalQueues();

    // h := 0
    // ======
    Zeros( h, 2*(n-1), 1 );

    // Solve the affine QP
    // ===================
    DistMultiVec<Real> xHat(grid), y(grid), z(grid), s(grid);
    QP( Q, A, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    x = xHat( IR(0,n), ALL );
}

#define PROTO(Real) \
  template void TV \
  ( const AbstractDistMatrix<Real>& b, \
          Real lambda, \
          AbstractDistMatrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void TV \
  ( const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void TV \
  ( const DistMultiVec<Real>& b, \
          Real lambda, \
          DistMultiVec<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
