/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// Basis pursuit denoising seeks the solution to the optimization problem
//
//   min (1/2) || b - A x ||_2^2 + lambda || x ||_1.
//
// Real instances of the problem are expressable as a Quadratic Program [1] via 
// the transformation
//
//   min lambda 1^T [u;v] + (1/2) r^T r
//   s.t. [A, -A] [u; v] + r = b, [u; v] >= 0.
//
// When expressed in affine conic form, the above expression becomes
//
//   min (1/2) [u;v;r]^T | 0 0 0 | | u | + lambda [1;1;0]^T | u |
//                       | 0 0 0 | | v |                    | v |
//                       | 0 0 I | | r |                    | r |
//
//   s.t. [A,-A,I] [u;v;r] = b, 
//
//        | -I  0 0 | | u | + s = | 0 |, s >= 0.
//        |  0 -I 0 | | v |       | 0 |
//                    | r | 
//
// Due to the linear transformation within the affine conic constraint,
// 
//   | -I  0 0 |
//   |  0 -I 0 |,
//
// being both sparse and exceedingly simple to analytically manipulate, 
// the dense variants of this algorithm will be unnecessarily slow relative
// to tailored algorithms (even without considering the use of iterative
// solvers for the KKT system exploiting fast algorithms for applying A).
//
// [1] Scott S. Chen, David L. Donoho, and Michael A. Saunders,
//     "Atomic Decomposition by Basis Pursuit",
//     SIAM Review, Vol. 43, No. 1, pp. 129--159, 2001

// TODO: Extend the existing QP control parameters
// TODO: Extend the (upcoming) SOCP control parameters

namespace El {

template<typename Real>
void BPDN
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BPDN"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m);

    Matrix<Real> Q, c, AHat, G, h;

    // Q := | 0 0 0 |
    //      | 0 0 0 |
    //      | 0 0 I |
    // ==============
    Zeros( Q, 2*n+m, 2*n+m );
    auto Qrr = Q( rInd, rInd );
    FillDiagonal( Qrr, Real(1) );

    // c := lambda*[1;1;0]
    // ===================
    Zeros( c, 2*n+m, 1 );
    auto cuv = c( IR(0,2*n), IR(0,1) );
    Fill( cuv, lambda );

    // \hat A := [A, -A, I]
    // ====================
    Zeros( AHat, m, 2*n+m );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    auto AHatr = AHat( IR(0,m), rInd );
    AHatu = A;
    Axpy( Real(-1), A, AHatv );
    FillDiagonal( AHatr, Real(1) );

    // G := | -I  0 0 |
    //      |  0 -I 0 |
    // ================
    Zeros( G, 2*n, 2*n+m );
    FillDiagonal( G, Real(-1) );

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine QP
    // ===================
    Matrix<Real> xHat, y, z, s;
    QP( Q, AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, IR(0,1) );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void BPDN
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BPDN"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m);
    DistMatrix<Real> Q(g), c(g), AHat(g), G(g), h(g);

    // Q := | 0 0 0 |
    //      | 0 0 0 |
    //      | 0 0 I |
    // ==============
    Zeros( Q, 2*n+m, 2*n+m );
    auto Qrr = Q( rInd, rInd );
    FillDiagonal( Qrr, Real(1) );

    // c := lambda*[1;1;0]
    // ===================
    Zeros( c, 2*n+m, 1 );
    auto cuv = c( IR(0,2*n), IR(0,1) );
    Fill( cuv, lambda );

    // \hat A := [A, -A, I]
    // ====================
    Zeros( AHat, m, 2*n+m );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    auto AHatr = AHat( IR(0,m), rInd );
    AHatu = A;
    Axpy( Real(-1), A, AHatv );
    FillDiagonal( AHatr, Real(1) );

    // G := | -I  0 0 |
    //      |  0 -I 0 |
    // ================
    Zeros( G, 2*n, 2*n+m );
    FillDiagonal( G, Real(-1) );

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine QP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    QP( Q, AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    Copy( xHat( uInd, IR(0,1) ), x );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void BPDN
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BPDN"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m);
    SparseMatrix<Real> Q, AHat, G;
    Matrix<Real> c, h;

    // Q := | 0 0 0 |
    //      | 0 0 0 |
    //      | 0 0 I |
    // ==============
    Zeros( Q, 2*n+m, 2*n+m );
    Q.Reserve( m );
    for( Int e=0; e<m; ++e )
        Q.QueueUpdate( 2*n+e, 2*n+e, Real(1) );
    Q.MakeConsistent();

    // c := lambda*[1;1;0]
    // ===================
    Zeros( c, 2*n+m, 1 );
    auto cuv = c( IR(0,2*n), IR(0,1) );
    Fill( cuv, lambda );

    // \hat A := [A, -A, I]
    // ====================
    const Int numEntriesA = A.NumEntries();
    Zeros( AHat, m, 2*n+m );
    AHat.Reserve( 2*numEntriesA+m );
    for( Int e=0; e<numEntriesA; ++e )
    {
        AHat.QueueUpdate( A.Row(e), A.Col(e),    A.Value(e) );
        AHat.QueueUpdate( A.Row(e), A.Col(e)+n, -A.Value(e) );
    }
    for( Int e=0; e<m; ++e )
        AHat.QueueUpdate( e, e+2*n, Real(1) );
    AHat.MakeConsistent();

    // G := | -I  0 0 |
    //      |  0 -I 0 |
    // ================
    Zeros( G, 2*n, 2*n+m );
    G.Reserve( 2*m );
    for( Int e=0; e<2*m; ++e )
        G.QueueUpdate( e, e, Real(-1) );
    G.MakeConsistent();

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine QP
    // ===================
    Matrix<Real> xHat, y, z, s;
    QP( Q, AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, IR(0,1) );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void BPDN
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        Real lambda,
        DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("BPDN"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> Q(comm), AHat(comm), G(comm);
    DistMultiVec<Real> c(comm), h(comm);

    // Q := | 0 0 0 |
    //      | 0 0 0 |
    //      | 0 0 I |
    // ==============
    Zeros( Q, 2*n+m, 2*n+m );
    {
        Int numLocalUpdates = 0;
        for( Int iLoc=0; iLoc<Q.LocalHeight(); ++iLoc )
            if( Q.GlobalRow(iLoc) >= 2*n )
                ++numLocalUpdates;
        Q.Reserve( numLocalUpdates );
        for( Int iLoc=0; iLoc<Q.LocalHeight(); ++iLoc )
            if( Q.GlobalRow(iLoc) >= 2*n )
                Q.QueueLocalUpdate( iLoc, Q.GlobalRow(iLoc), Real(1) );
        Q.MakeConsistent();
    }

    // c := lambda*[1;1;0]
    // ===================
    Zeros( c, 2*n+m, 1 );
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
        if( c.GlobalRow(iLoc) < 2*n )
            c.SetLocal( iLoc, 0, lambda );

    // \hat A := [A, -A, I]
    // ====================
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    const Int numLocalEntriesA = A.NumLocalEntries();
    Zeros( AHat, m, 2*n+m );
    AHat.Reserve( 2*numLocalEntriesA+AHat.LocalHeight() );
    for( Int e=0; e<numLocalEntriesA; ++e )
    {
        AHat.QueueLocalUpdate
        ( A.Row(e)-A.FirstLocalRow(), A.Col(e),    A.Value(e) );
        AHat.QueueLocalUpdate
        ( A.Row(e)-A.FirstLocalRow(), A.Col(e)+n, -A.Value(e) );
    }
    for( Int iLoc=0; iLoc<AHat.LocalHeight(); ++iLoc )
    {
        const Int i = AHat.GlobalRow(iLoc);
        AHat.QueueLocalUpdate( iLoc, i+2*n, Real(1) );
    }
    AHat.MakeConsistent();

    // G := | -I  0 0 |
    //      |  0 -I 0 |
    // ================
    Zeros( G, 2*n, 2*n+m );
    G.Reserve( G.LocalHeight() );
    for( Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
    {
        const Int i = G.GlobalRow(iLoc);
        G.QueueLocalUpdate( iLoc, i, Real(-1) );
    }
    G.MakeConsistent();

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine QP
    // ===================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm), s(comm);
    QP( Q, AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
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
        else if( i < 2*n )
            ++sendCounts[ x.RowOwner(i-n) ];
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
        else if( i < 2*n )
        {
            const int owner = x.RowOwner(i-n);
            sSendBuf[offsets[owner]] = i-n;
            vSendBuf[offsets[owner]] = -xHat.GetLocal(iLoc,0);
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
  template void BPDN \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void BPDN \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          Real lambda, \
          AbstractDistMatrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void BPDN \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void BPDN \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          Real lambda, \
          DistMultiVec<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
