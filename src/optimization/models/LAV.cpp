/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// Least Absolute Value (LAV) regression minimizes the one norm of the 
// residual of a system of equations, i.e.,
//
//   min || A x - b ||_1.
//
// Real instances of the problem are expressable as a Linear Program [1] via 
//
//   min 1^T (u + v)
//   s.t. [A, I, -I] [x; u; v] = b, u, v >= 0
//
// which, in affine standard form, becomes
//
//   min [0; 1; 1]^T [x; u; v]
//   s.t. [A, I, -I] | x | = b, | 0 -I  0 | | x | <= | 0 |.
//                   | u |      | 0  0 -I | | u |    | 0 |
//                   | v |                  | v |
//
// [1] 
//   A. Charnes, W. W. Cooper, and R. O. Ferguson,
//   "Optimal estimation of executive compensation by linear programming",
//   Management Science, Vol. 1, No. 2, pp. 138--151, 1955.
//

namespace El {

template<typename Real>
void LAV
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LAV"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> xInd(0,n), uInd(n,n+m), vInd(n+m,n+2*m);
    Matrix<Real> c, AHat, G, h;
    
    // c := [0;1;1]
    // ============
    Zeros( c, n+2*m, 1 );
    auto cuv = c( IR(n,n+2*m), IR(0,1) );
    Fill( cuv, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( AHat, m, n+2*m );
    auto AHatx = AHat( IR(0,m), xInd );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    AHatx = A;
    FillDiagonal( AHatu, Real( 1) );
    FillDiagonal( AHatv, Real(-1) );

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( G, 2*m, n+2*m );
    auto Guv = G( IR(0,2*m), IR(n,n+2*m) );
    FillDiagonal( Guv, Real(-1) );

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( h, 2*m, 1 );

    // Solve the affine linear program
    // ===============================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // Extract x
    // ==========
    x = xHat( xInd, IR(0,1) );
}

template<typename Real>
void LAV
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LAV"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> xInd(0,n), uInd(n,n+m), vInd(n+m,n+2*m);
    DistMatrix<Real> c(g), AHat(g), G(g), h(g);

    // c := [0;1;1]
    // ============
    Zeros( c, n+2*m, 1 );
    auto cuv = c( IR(n,n+2*m), IR(0,1) );
    Fill( cuv, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( AHat, m, n+2*m );
    auto AHatx = AHat( IR(0,m), xInd );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    AHatx = A;
    FillDiagonal( AHatu, Real( 1) );
    FillDiagonal( AHatv, Real(-1) );

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( G, 2*m, n+2*m );
    auto Guv = G( IR(0,2*m), IR(n,n+2*m) );
    FillDiagonal( Guv, Real(-1) );

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( h, 2*m, 1 );

    // Solve the affine linear program
    // ===============================
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    LP( AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // Extract x
    // =========
    Copy( xHat( xInd, IR(0,1) ), x );
}

template<typename Real>
void LAV
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LAV"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> xInd(0,n), uInd(n,n+m), vInd(n+m,n+2*m);
    SparseMatrix<Real> AHat, G;
    Matrix<Real> c, h;
    
    // c := [0;1;1]
    // ============
    Zeros( c, n+2*m, 1 );
    auto cuv = c( IR(n,n+2*m), IR(0,1) );
    Fill( cuv, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( AHat, m, n+2*m );
    const Int numEntriesA = A.NumEntries();
    AHat.Reserve( numEntriesA + 2*m );
    for( Int e=0; e<numEntriesA; ++e )
        AHat.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int i=0; i<m; ++i )
    {
        AHat.QueueUpdate( i, i+n,   Real( 1) );
        AHat.QueueUpdate( i, i+n+m, Real(-1) );
    }
    AHat.MakeConsistent();

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( G, 2*m, n+2*m );
    G.Reserve( G.Height() );
    for( Int i=0; i<2*m; ++i )
        G.QueueUpdate( i, i+n, Real(-1) );
    G.MakeConsistent();

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( h, 2*m, 1 );

    // Solve the affine linear program
    // ===============================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // Extract x
    // =========
    x = xHat( xInd, IR(0,1) );
}

template<typename Real>
void LAV
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LAV"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> AHat(comm), G(comm);
    DistMultiVec<Real> c(comm), h(comm);

    // c := [0;1;1]
    // ============
    Zeros( c, n+2*m, 1 );
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
        if( c.GlobalRow(iLoc) >= n )
            c.SetLocal( iLoc, 0, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( AHat, m, n+2*m );
    const Int numLocalEntriesA = A.NumLocalEntries();
    AHat.Reserve( numLocalEntriesA + 2*AHat.LocalHeight() );
    for( Int e=0; e<numLocalEntriesA; ++e )
        AHat.QueueLocalUpdate
        ( A.Row(e)-AHat.FirstLocalRow(), A.Col(e), A.Value(e) );
    for( Int iLoc=0; iLoc<AHat.LocalHeight(); ++iLoc )
    {
        const Int i = AHat.GlobalRow(iLoc);
        AHat.QueueLocalUpdate( iLoc, i+n,   Real( 1) );
        AHat.QueueLocalUpdate( iLoc, i+n+m, Real(-1) );
    }
    AHat.MakeConsistent();

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( G, 2*m, n+2*m );
    G.Reserve( G.LocalHeight() );
    for( Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
        G.QueueLocalUpdate( iLoc, G.GlobalRow(iLoc)+n, Real(-1) );
    G.MakeConsistent();

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( h, 2*m, 1 );

    // Solve the affine QP
    // ===================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm), s(comm);
    LP( AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // Extract x
    // =========
    Zeros( x, n, 1 );
    // Determine the send and recv counts/offsets
    // ------------------------------------------
    const Int commSize = mpi::Size(comm);
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
    // Pack the data 
    // -------------
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
    // Exchange the data
    // -----------------
    vector<Int> sRecvBuf(totalRecv);
    vector<Real> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    // Unpack the data
    // ---------------
    for( Int e=0; e<totalRecv; ++e )
        x.UpdateLocal( sRecvBuf[e]-x.FirstLocalRow(), 0, vRecvBuf[e] );
}

#define PROTO(Real) \
  template void LAV \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LAV \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LAV \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LAV \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          DistMultiVec<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
