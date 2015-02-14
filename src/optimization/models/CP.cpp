/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// A Chebyshev point (CP) minimizes the supremum norm of A x - b, i.e.,
//
//   min || A x - b ||_oo.
//
// Real instances of the problem are expressable as a Linear Program via 
//
//   min t
//   s.t. -t <= A x - b <= t,
//
// which, in affine standard form, becomes
//
//   min [0; 1]^T [x; t]
//   s.t. |  A  -1 | | x | <= |  b |
//        | -A  -1 | | t |    | -b |
//
// NOTE: There is likely an appropriate citation, but the derivation is 
//       trivial. If one is found, it will be added.

namespace El {

template<typename Real>
void CP
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> c, AHat, bHat, G, h;

    // c := [zeros(n,1);1]
    // ===================
    Zeros( c, n+1, 1 );
    c.Set( n, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ====================== 
    Zeros( AHat, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( bHat, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( G, 2*m, n+1 );
    auto GTL = G( IR(0,m),   IR(0,n)   );
    auto GBL = G( IR(m,2*m), IR(0,n)   );
    auto gR =  G( IR(0,2*m), IR(n,n+1) );
    GTL = A;
    Axpy( Real(-1), GTL, GBL );
    Fill( gR, Real(-1) );

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( h, 2*m, 1 );
    auto hT = h( IR(0,m),   IR(0,1) );
    auto hB = h( IR(m,2*m), IR(0,1) );
    hT = b;
    Axpy( Real(-1), hT, hB );

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    x = xHat( IR(0,n), IR(0,1) );
}

template<typename Real>
void CP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<Real> c(g), AHat(g), bHat(g), G(g), h(g);

    // c := [zeros(n,1);1]
    // ===================
    Zeros( c, n+1, 1 );
    c.Set( n, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ====================== 
    Zeros( AHat, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( bHat, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( G, 2*m, n+1 );
    auto GTL = G( IR(0,m),   IR(0,n)   );
    auto GBL = G( IR(m,2*m), IR(0,n)   );
    auto gR =  G( IR(0,2*m), IR(n,n+1) );
    GTL = A;
    Axpy( Real(-1), GTL, GBL );
    Fill( gR, Real(-1) );

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( h, 2*m, 1 );
    auto hT = h( IR(0,m),   IR(0,1) );
    auto hB = h( IR(m,2*m), IR(0,1) );
    hT = b;
    Axpy( Real(-1), hT, hB );

    // Solve the affine LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    Copy( xHat( IR(0,n), IR(0,1) ), x );
}

template<typename Real>
void CP
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    const Int m = A.Height();
    const Int n = A.Width();
    SparseMatrix<Real> AHat, G;
    Matrix<Real> c, bHat, h;

    // c := [zeros(n,1);1]
    // ===================
    Zeros( c, n+1, 1 );
    c.Set( n, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ====================== 
    Zeros( AHat, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( bHat, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( G, 2*m, n+1 );
    const Int numEntriesA = A.NumEntries();
    G.Reserve( 2*numEntriesA + 2*m );
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Real value = A.Value(e);
        G.QueueUpdate( i,   j,  value );
        G.QueueUpdate( i+m, j, -value );
    }
    for( Int i=0; i<2*m; ++i )
        G.QueueUpdate( i, n, Real(-1) );
    G.MakeConsistent();

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( h, 2*m, 1 );
    auto hT = h( IR(0,m),   IR(0,1) );
    auto hB = h( IR(m,2*m), IR(0,1) );
    hT = b;
    Axpy( Real(-1), hT, hB );

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    x = xHat( IR(0,n), IR(0,1) );
}

template<typename Real>
void CP
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size(comm);
    DistSparseMatrix<Real> AHat(comm), G(comm);
    DistMultiVec<Real> c(comm), bHat(comm), h(comm);

    // c := [zeros(n,1);1]
    // ===================
    Zeros( c, n+1, 1 );
    // Have the process which owns the last entry set it to one
    if( c.GlobalRow(c.LocalHeight()-1) == n )
        c.SetLocal( c.LocalHeight()-1, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ====================== 
    Zeros( AHat, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( bHat, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( G, 2*m, n+1 );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int e=0; e<A.NumLocalEntries(); ++e )
        {
            const Int i = A.Row(e);
            ++sendCounts[ G.RowOwner(i)   ];
            ++sendCounts[ G.RowOwner(i+m) ];
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets, recvOffsets;
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack 
        // ----
        vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        vector<Real> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int e=0; e<A.NumLocalEntries(); ++e )
        {
            const Int i = A.Row(e);
            const Int j = A.Col(e);
            const Real value = A.Value(e);
            int owner = G.RowOwner(i);    
            sSendBuf[offsets[owner]] = i;
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = value;
            ++offsets[owner];
            owner = G.RowOwner(i+m);
            sSendBuf[offsets[owner]] = i+m;
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = -value;
            ++offsets[owner];
        }
        // Exchange
        // --------
        vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
        vector<Real> vRecvBuf(totalRecv);
        mpi::AllToAll
        ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        mpi::AllToAll
        ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
          vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
        // Unpack
        // ------
        G.Reserve( totalRecv+G.LocalHeight() );
        for( Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
            G.QueueLocalUpdate( iLoc, n, Real(-1) );
        for( Int e=0; e<totalRecv; ++e )
            G.QueueLocalUpdate
            ( sRecvBuf[e]-G.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] );
        G.MakeConsistent();
    }

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( h, 2*m, 1 );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
        {
            const Int i = b.GlobalRow(iLoc);
            ++sendCounts[ h.RowOwner(i)   ];
            ++sendCounts[ h.RowOwner(i+m) ];
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets, recvOffsets;
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        vector<Int> sSendBuf(totalSend);
        vector<Real> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
        {
            const Int i = b.GlobalRow(iLoc);
            const Real value = b.GetLocal(iLoc,0);
            int owner = h.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            vSendBuf[offsets[owner]] = value;
            ++offsets[owner];
            owner = h.RowOwner(i+m);
            sSendBuf[offsets[owner]] = i+m;
            vSendBuf[offsets[owner]] = -value;
            ++offsets[owner];
        }
        // Exchange
        // --------
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
            h.SetLocal( sRecvBuf[e]-h.FirstLocalRow(), 0, vRecvBuf[e] );
    }

    // Solve the affine LP
    // ===================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm), s(comm);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    Zeros( x, n, 1 );
    {
        // Compute the metadata
        // --------------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            if( i < n )
                ++sendCounts[ x.RowOwner(i) ];
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets, recvOffsets;
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        vector<Int> sSendBuf(totalSend);
        vector<Real> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            if( i < n )
            {
                const Real value = xHat.GetLocal(iLoc,0);
                int owner = x.RowOwner(i);
                sSendBuf[offsets[owner]] = i;
                vSendBuf[offsets[owner]] = value;
                ++offsets[owner];
            }
        }
        // Exchange
        // --------
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
  template void CP \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          DistMultiVec<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
