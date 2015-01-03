/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// Basis Pursuit is the search for a solution to A x = b which minimizes
// the one norm of x, i.e.,
//
//   min || x ||_1
//   s.t. A x = b.
//
// Real instances of the problem are expressable as a Linear Program [1] via 
// the transformation
//
//   min 1^T u,
//   s.t. A x = b, -u <= x <= u,
//
// which can be directly expressed in terms of the affine conic form LP
//
//   min c^T \hat x,
//   s.t. \hat A \hat x = b, G \hat x + s = h, s >= 0,
//
// by setting 
//
//   \hat x = [u; x], 
//        c = [1; 0],
//   \hat A = [0, A],
//        G = [-I, -I; I, -I], and
//        h = 0.
//
// Complex instances of Basis Pursuit can be solved via a Second-Order Cone
// Program where each second-order cone is of dimension three:
//
//     u >= sqrt(real(x)^2 + imag(x)^2).
//
// [1] Jiri Matousek and Bernd Gartner,
//     "Understanding and Using Linear Programming",
//     Section 8. More Applications

// TODO: Extend the existing LP control parameters
// TODO: Extend the (upcoming) SOCP control parameters

namespace El {

template<typename Real>
void BasisPursuit
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> c, xHat, AHat, G, h;

    // c := [1; 0]
    // ===========
    Zeros( c, 2*n, 1 );
    auto cT = c( IR(0,n),   IR(0,1) );
    Ones( cT, n, 1 );

    // \hat A := [0, A]
    // ================
    Zeros( AHat, m, 2*n );
    auto AHatR = AHat( IR(0,m), IR(n,2*n) );
    AHatR = A;

    // G := [-I, -I; I, -I]
    // ====================
    Zeros( G, 2*n, 2*n );
    auto GTL = G( IR(0,n),   IR(0,n) ); auto GTR = G( IR(0,n),   IR(n,2*n) );
    auto GBL = G( IR(n,2*n), IR(0,n) ); auto GBR = G( IR(n,2*n), IR(n,2*n) );
    Identity( GTL, n, n );
    Scale( Real(-1), GTL );
    Identity( GTR, n, n );
    Scale( Real(-1), GTR );
    Identity( GBL, n, n );
    Identity( GBR, n, n );
    Scale( Real(-1), GBR );

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine LP
    // ===================
    Matrix<Real> y, z, s;
    LP( AHat, G, b, c, h, xHat, y, z, s );

    // Extract x from \hat x = [u; x]
    // ==============================
    x = xHat( IR(n,2*n), IR(0,1) );
}

template<typename Real>
void BasisPursuit
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<Real> c(g), xHat(g), AHat(g), G(g), h(g);

    // c := [1; 0]
    // ===========
    Zeros( c, 2*n, 1 );
    auto cT = c( IR(0,n),   IR(0,1) );
    Ones( cT, n, 1 );

    // \hat A := [0, A]
    // ================
    Zeros( AHat, m, 2*n );
    auto AHatR = AHat( IR(0,m), IR(n,2*n) );
    AHatR = A;

    // G := [-I, -I; I, -I]
    // ====================
    Zeros( G, 2*n, 2*n );
    auto GTL = G( IR(0,n),   IR(0,n) ); auto GTR = G( IR(0,n),   IR(n,2*n) );
    auto GBL = G( IR(n,2*n), IR(0,n) ); auto GBR = G( IR(n,2*n), IR(n,2*n) );
    Identity( GTL, n, n );
    Scale( Real(-1), GTL );
    Identity( GTR, n, n );
    Scale( Real(-1), GTR );
    Identity( GBL, n, n );
    Identity( GBR, n, n );
    Scale( Real(-1), GBR );

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine LP
    // ===================
    DistMatrix<Real> y(g), z(g), s(g);
    LP( AHat, G, b, c, h, xHat, y, z, s );

    // Extract x from \hat x = [u; x]
    // ==============================
    x = xHat( IR(n,2*n), IR(0,1) );
}

template<typename Real>
void BasisPursuit
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    const Int m = A.Height();
    const Int n = A.Width();
    SparseMatrix<Real> AHat, G;
    Matrix<Real> c, xHat, h;

    // c := [1; 0]
    // ===========
    Zeros( c, 2*n, 1 );
    auto cT = c( IR(0,n),   IR(0,1) );
    Ones( cT, n, 1 );

    // \hat A := [0, A]
    // ================
    const Int numEntriesA = A.NumEntries();
    Zeros( AHat, m, 2*n );
    AHat.Reserve( numEntriesA );
    for( Int e=0; e<numEntriesA; ++e )
        AHat.QueueUpdate( A.Row(e), A.Col(e)+n, A.Value(e) );
    AHat.MakeConsistent();

    // G := [-I, -I; I, -I]
    // ====================
    Zeros( G, 2*n, 2*n );
    G.Reserve( 4*n );
    for( Int e=0; e<n; ++e )
    {
        G.QueueUpdate( e,   e,   Real(-1) );
        G.QueueUpdate( e+n, e,   Real(+1) );
        G.QueueUpdate( e,   e+n, Real(-1) ); 
        G.QueueUpdate( e+n, e+n, Real(-1) );
    }
    G.MakeConsistent();

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine LP
    // ===================
    Matrix<Real> y, z, s;
    LP( AHat, G, b, c, h, xHat, y, z, s );

    // Extract x from \hat x = [u, x]
    // ==============================
    Zeros( x, n, 1 );
    for( Int i=0; i<n; ++i )
        x.Set( i, 0, xHat.Get(i+n,0) );
}

template<typename Real>
void BasisPursuit
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("BasisPursuit"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> AHat(comm), G(comm);
    DistMultiVec<Real> c(comm), xHat(comm), h(comm);

    // c := [1; 0]
    // ===========
    Zeros( c, 2*n, 1 );
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
    {
        const Int i = c.GlobalRow(iLoc);
        if( i < n )
            c.SetLocal( iLoc, 0, Real(1) );
    }

    // \hat A := [0, A]
    // ================
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    const Int numLocalEntriesA = A.NumLocalEntries();
    Zeros( AHat, m, 2*n );
    AHat.Reserve( numLocalEntriesA );
    for( Int e=0; e<numLocalEntriesA; ++e )
        AHat.QueueLocalUpdate
        ( A.Row(e)-A.FirstLocalRow(), A.Col(e)+n, A.Value(e) );
    AHat.MakeConsistent();

    // G := [-I, -I; I, -I]
    // ====================
    Zeros( G, 2*n, 2*n );
    G.Reserve( 2*G.LocalHeight() );
    for( Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
    {
        const Int i = G.GlobalRow(iLoc);
        G.QueueLocalUpdate( iLoc, i,   Real(-1) );
        if( i < n )
            G.QueueLocalUpdate( iLoc, i+n, Real(-1) );
        else
            G.QueueLocalUpdate( iLoc, i-n, Real(+1) );
    }
    G.MakeConsistent();

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine LP
    // ===================
    DistMultiVec<Real> y(comm), z(comm), s(comm);
    LP( AHat, G, b, c, h, xHat, y, z, s );

    // Extract x from \hat x = [u; x]
    // ==============================
    x.Resize( n, 1 );
    // Determine the send and recv counts/offsets
    // ------------------------------------------
    const Int commSize = mpi::Size(comm);
    std::vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i >= n )
            ++sendCounts[ x.RowOwner(i-n) ];
    }
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    std::vector<int> sendOffsets, recvOffsets;
    const int totalSend = Scan( sendCounts, sendOffsets );
    const int totalRecv = Scan( recvCounts, recvOffsets );
    // Pack the data 
    // -------------
    std::vector<Int> sSendBuf(totalSend);
    std::vector<Real> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i >= n )
        {
            const int owner = x.RowOwner(i-n);
            sSendBuf[offsets[owner]] = i;
            vSendBuf[offsets[owner]] = xHat.GetLocal(iLoc,0);
            ++offsets[owner];
        }
    }
    // Exchange the data
    // -----------------
    std::vector<Int> sRecvBuf(totalRecv);
    std::vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    // Unpack the data
    // ---------------
    for( Int e=0; e<totalRecv; ++e )
        x.SetLocal( sRecvBuf[e]-x.FirstLocalRow(), 0, vRecvBuf[e] );
}

#define PROTO(Real) \
  template void BasisPursuit \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x ); \
  template void BasisPursuit \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x ); \
  template void BasisPursuit \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x ); \
  template void BasisPursuit \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          DistMultiVec<Real>& x );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
