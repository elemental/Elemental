/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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
    DEBUG_ONLY(CallStackEntry cse("TV"))
    mpi::Comm comm = b.Grid().Comm();
    DistMultiVec<Real> bDMV(comm), xDMV(comm);
    bDMV = b;
    xDMV = x;
    TV( bDMV, lambda, xDMV, ctrl );
    x = xDMV;
}

template<typename Real>
void TV
( const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("TV"))
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
    Q.MakeConsistent();

    // c := [-b;lambda]
    // =================
    Zeros( c, 2*n-1, 1 );
    auto cx = c( xInd, IR(0,1) ); 
    cx = b; Scale( Real(-1), cx );
    auto ct = c( tInd, IR(0,1) );
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
    G.MakeConsistent();

    // h := 0
    // ======
    Zeros( h, 2*(n-1), 1 );

    // Solve the affine QP
    // ===================
    Matrix<Real> xHat, y, z, s;
    QP( Q, A, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    x = xHat( xInd, IR(0,1) );
}

template<typename Real>
void TV
( const DistMultiVec<Real>& b, 
        Real lambda,
        DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("TV"))
    const Int n = b.Height();
    mpi::Comm comm = b.Comm();
    const int commSize = mpi::Size(comm);

    DistSparseMatrix<Real> Q(comm), A(comm), G(comm);
    DistMultiVec<Real> c(comm), bHat(comm), h(comm);

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
    Q.MakeConsistent();

    // c := [-b;lambda]
    // =================
    Zeros( c, 2*n-1, 1 );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
        {
            const Int i = b.GlobalRow(iLoc);
            ++sendCounts[ c.RowOwner(i) ];
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets, recvOffsets;
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack -b
        // -------
        vector<Int> sSendBuf(totalSend);
        vector<Real> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
        {
            const Int i = b.GlobalRow(iLoc);
            const int owner = c.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            vSendBuf[offsets[owner]] = -b.GetLocal(iLoc,0);
            ++offsets[owner];
        }
        // Redistribute
        // ------------
        vector<Int> sRecvBuf(totalRecv);
        vector<Real> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        for( Int e=0; e<totalRecv; ++e )
            c.SetLocal( sRecvBuf[e]-c.FirstLocalRow(), 0, vRecvBuf[e] );
    }
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
    {
        const Int i = c.GlobalRow(iLoc);
        if( i > n )
            c.SetLocal( iLoc, 0, lambda );
    }

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
            G.QueueLocalUpdate( iLoc, i,   Real(1) );
            G.QueueLocalUpdate( iLoc, i+1, Real(-1) );
            G.QueueLocalUpdate( iLoc, i+n, Real(-1) );
        }
        else
        {
            G.QueueLocalUpdate( iLoc, i-(n-1),   Real(-1) );
            G.QueueLocalUpdate( iLoc, i+1-(n-1), Real(1) );
            G.QueueLocalUpdate( iLoc, i+n-(n-1), Real(-1) );
        }
    }
    G.MakeConsistent();

    // h := 0
    // ======
    Zeros( h, 2*(n-1), 1 );

    // Solve the affine QP
    // ===================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm), s(comm);
    QP( Q, A, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    {
        Zeros( x, n, 1 );
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            if( i < n )
                ++sendCounts[ x.RowOwner(i) ];
            else
                break;
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets, recvOffsets;
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack x
        // ------
        vector<Int> sSendBuf(totalSend);
        vector<Real> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            if( i < n )
            {
                const int owner = x.RowOwner(i);
                sSendBuf[offsets[owner]] = i;
                vSendBuf[offsets[owner]] = xHat.GetLocal(iLoc,0);
                ++offsets[owner];
            }
            else
                break;
        }
        // Redistribute
        // ------------
        vector<Int> sRecvBuf(totalRecv);
        vector<Real> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        for( Int e=0; e<totalRecv; ++e )
            x.SetLocal( sRecvBuf[e]-x.FirstLocalRow(), 0, vRecvBuf[e] );
    }
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
#include "El/macros/Instantiate.h"

} // namepace elem
