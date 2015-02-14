/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// The Dantzig selector [1] seeks the solution to the problem
//
//   min || x ||_1 s.t. || A^T (b - A x) ||_oo <= lambda.
//
// One motivation is that it is an LP alternative to Basis pursuit denoising
// (BPDN), in that it seeks a sparse approximate minimizer of b - A x.
//
// Friendlander and Saunders propose [2] two LP formulations:
//
// (DS1)
//   min 1^T (u+v)
//   s.t. [A^T A, -A^T A, I] | u | = A^T b, u,v >= 0, || t ||_oo <= lambda,
//                           | v |
//                           | t |
//
// (DS2)
//   min 1^T (u+v)
//   s.t. | A, -A,  I,  0 | | u | = | b |, u,v >= 0, || t ||_oo <= lambda,
//        | 0,  0, A^T, I | | v |   | 0 |
//                          | r |
//                          | t |
//
//   which, in affine conic form, becomes
//
//   min [1;1;0;0]^T [u;v;r;t]
//   s.t. | A, -A,  I,  0 | | u | = | b |,
//        | 0,  0, A^T, I | | v |   | 0 |
//                          | r |
//                          | t |
//
//        | -I  0 0 0 | | u | + s = |    0     |, s >= 0.
//        |  0 -I 0 0 | | v |       |    0     |
//        |  0  0 0 I | | r |       | lambda e |
//        |  0  0 0-I | | t |       | lambda e |
//
// For dense and sparse matrices we respectively default to (DS1) and (DS2).
//
// TODO: 
//  Add the ability to switch between the (DS1) and (DS2) affine LP
//  formulations described in [2].
//
// [1] 
//   Emmanuel Candes and Terence Tao,
//   "The Dantzig selector: Statistical estimation when p is much 
//    larger than n",
//   The Annals of Statistics, pp. 2313--2351, 2007.
//
// [2] 
//   Michael Friedlander and Michael Saunders,
//  "Discussion: The Dantzig selector: Statistical estimation when p is much
//   larger than n",
//  The Annals of Statistics, Vol. 35, No. 6, pp. 2385--2391, 2007.

namespace El {

namespace ds {

template<typename Real>
void Var1
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ds::Var1"))
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), tInd(2*n,3*n);
    Matrix<Real> c, AHat, bHat, G, h;

    // c := [1;1;0]
    // ==============
    Zeros( c, 3*n, 1 );
    auto cuv = c( IR(0,2*n), IR(0,1) );
    Fill( cuv, Real(1) );

    // \hat A := [ A^T A, -A^T A,  I ]
    // ===============================
    Zeros( AHat, n, 3*n );
    auto AHatu = AHat( IR(0,n), uInd );
    auto AHatv = AHat( IR(0,n), vInd );
    auto AHatt = AHat( IR(0,n), tInd );
    Herk( LOWER, TRANSPOSE, Real(1), A, Real(0), AHatu );
    MakeSymmetric( LOWER, AHatu );
    Axpy( Real(-1), AHatu, AHatv );
    FillDiagonal( AHatt, Real(1) );

    // \hat b := A^T b
    // ===============
    Zeros( bHat, n, 1 );
    Gemv( TRANSPOSE, Real(1), A, b, Real(0), bHat );

    // G := | -I  0  0 |
    //      |  0 -I  0 |
    //      |  0  0  I |
    //      |  0  0 -I |
    // =================
    Zeros( G, 4*n, 3*n );
    auto Guv = G( IR(0,2*n), IR(0,2*n) );
    FillDiagonal( Guv, Real(-1) );
    auto G2t = G( IR(2*n,3*n), tInd );
    auto G3t = G( IR(3*n,4*n), tInd );
    FillDiagonal( G2t, Real( 1) );
    FillDiagonal( G3t, Real(-1) );

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto ht = h( IR(2*n,4*n), IR(0,1) );
    Fill( ht, lambda );

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, IR(0,1) );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void Var1
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& b, 
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ds::Var1"))

    auto APtr = ReadProxy<Real,MC,MR>(&APre);
    auto& A = *APtr;

    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> uInd(0,n), vInd(n,2*n), tInd(2*n,3*n);
    DistMatrix<Real> c(g), AHat(g), bHat(g), G(g), h(g);

    // c := [1;1;0]
    // ==============
    Zeros( c, 3*n, 1 );
    auto cuv = c( IR(0,2*n), IR(0,1) );
    Fill( cuv, Real(1) );

    // \hat A := [ A^T A, -A^T A,  I ]
    // ===============================
    Zeros( AHat, n, 3*n );
    auto AHatu = AHat( IR(0,n), uInd );
    auto AHatv = AHat( IR(0,n), vInd );
    auto AHatt = AHat( IR(0,n), tInd );
    Herk( LOWER, TRANSPOSE, Real(1), A, Real(0), AHatu );
    MakeSymmetric( LOWER, AHatu );
    Axpy( Real(-1), AHatu, AHatv );
    FillDiagonal( AHatt, Real(1) );

    // \hat b := A^T b
    // ===============
    Zeros( bHat, n, 1 );
    Gemv( TRANSPOSE, Real(1), A, b, Real(0), bHat );

    // G := | -I  0  0 |
    //      |  0 -I  0 |
    //      |  0  0  I |
    //      |  0  0 -I |
    // =================
    Zeros( G, 4*n, 3*n );
    auto Guv = G( IR(0,2*n), IR(0,2*n) );
    FillDiagonal( Guv, Real(-1) );
    auto G2t = G( IR(2*n,3*n), tInd );
    auto G3t = G( IR(3*n,4*n), tInd );
    FillDiagonal( G2t, Real( 1) );
    FillDiagonal( G3t, Real(-1) );

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto ht = h( IR(2*n,4*n), IR(0,1) );
    Fill( ht, lambda );

    // Solve the affine LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    Copy( xHat( uInd, IR(0,1) ), x );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void Var2
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ds::Var2"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m), tInd(2*n+m,3*n+m);
    Matrix<Real> c, AHat, bHat, G, h;

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n+m, 1 );
    auto c_u = c( uInd, IR(0,1) );
    auto c_v = c( vInd, IR(0,1) );
    Ones( c_u, n, 1 );
    Ones( c_v, n, 1 );

    // \hat A := | A, -A,  I,  0 |
    //           | 0,  0, A^T, I |
    // ===========================
    Zeros( AHat, m+n, 3*n+m );
    auto AHat_0u = AHat( IR(0,m),   uInd );
    auto AHat_0v = AHat( IR(0,m),   vInd );
    auto AHat_0r = AHat( IR(0,m),   rInd );
    auto AHat_1r = AHat( IR(m,m+n), rInd );
    auto AHat_1t = AHat( IR(m,m+n), tInd );
    AHat_0u = A;
    AHat_0v = A; Scale( Real(-1), AHat_0v );
    Identity( AHat_0r, m, m );
    Transpose( A, AHat_1r );
    Identity( AHat_1t, n, n );

    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    auto b0 = bHat( IR(0,m), IR(0,1) );
    b0 = b;

    // G := | -I  0 0  0 |
    //      |  0 -I 0  0 |
    //      |  0  0 0  I |
    //      |  0  0 0 -I |
    // ===================
    Zeros( G, 4*n, 3*n+m );
    auto G0u = G( IR(0,    n), uInd );
    auto G1v = G( IR(n,  2*n), vInd );
    auto G2t = G( IR(2*n,3*n), tInd );
    auto G3t = G( IR(3*n,4*n), tInd );
    Identity( G0u, n, n ); Scale( Real(-1), G0u );
    Identity( G1v, n, n ); Scale( Real(-1), G1v );
    Identity( G2t, n, n );
    Identity( G3t, n, n ); Scale( Real(-1), G3t );

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto h2 = h( IR(2*n,3*n), IR(0,1) );
    auto h3 = h( IR(3*n,4*n), IR(0,1) );
    Fill( h2, lambda );
    Fill( h3, lambda );

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, IR(0,1) );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void Var2
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ds::Var2"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m), tInd(2*n+m,3*n+m);
    DistMatrix<Real> c(g), AHat(g), bHat(g), G(g), h(g);

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n+m, 1 );
    auto c_u = c( uInd, IR(0,1) );
    auto c_v = c( vInd, IR(0,1) );
    Ones( c_u, n, 1 );
    Ones( c_v, n, 1 );

    // \hat A := | A, -A,  I,  0 |
    //           | 0,  0, A^T, I |
    // ===========================
    Zeros( AHat, m+n, 3*n+m );
    auto AHat_0u = AHat( IR(0,m),   uInd );
    auto AHat_0v = AHat( IR(0,m),   vInd );
    auto AHat_0r = AHat( IR(0,m),   rInd );
    auto AHat_1r = AHat( IR(m,m+n), rInd );
    auto AHat_1t = AHat( IR(m,m+n), tInd );
    AHat_0u = A;
    AHat_0v = A; Scale( Real(-1), AHat_0v );
    Identity( AHat_0r, m, m );
    Transpose( A, AHat_1r );
    Identity( AHat_1t, n, n );

    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    auto b0 = bHat( IR(0,m), IR(0,1) );
    b0 = b;

    // G := | -I  0 0  0 |
    //      |  0 -I 0  0 |
    //      |  0  0 0  I |
    //      |  0  0 0 -I |
    // ===================
    Zeros( G, 4*n, 3*n+m );
    auto G0u = G( IR(0,    n), uInd );
    auto G1v = G( IR(n,  2*n), vInd );
    auto G2t = G( IR(2*n,3*n), tInd );
    auto G3t = G( IR(3*n,4*n), tInd );
    Identity( G0u, n, n ); Scale( Real(-1), G0u );
    Identity( G1v, n, n ); Scale( Real(-1), G1v );
    Identity( G2t, n, n );
    Identity( G3t, n, n ); Scale( Real(-1), G3t );

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto h2 = h( IR(2*n,3*n), IR(0,1) );
    auto h3 = h( IR(3*n,4*n), IR(0,1) );
    Fill( h2, lambda );
    Fill( h3, lambda );

    // Solve the affine LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    Copy( xHat( uInd, IR(0,1) ), x );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void Var2
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ds::Var2"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntriesA = A.NumEntries();
    const Range<Int> uInd(0,n), vInd(n,2*n);
    SparseMatrix<Real> AHat, G;
    Matrix<Real> c, bHat, h;

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n, 1 );
    auto cuv = c( IR(0,2*n), IR(0,1) );
    Fill( cuv, Real(1) );

    // \hat A := | A, -A,  I,  0 |
    //           | 0,  0, A^T, I |
    // ===========================
    Zeros( AHat, m+n, 3*n+m );
    AHat.Reserve( 3*numEntriesA + m + n );
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Real value = A.Value(e);
        AHat.QueueUpdate( i,   j,     value );
        AHat.QueueUpdate( i,   j+n,  -value );
        AHat.QueueUpdate( j+m, i+2*n, value );
    }
    for( Int e=0; e<m+n; ++e )
        AHat.QueueUpdate( e, e+2*n, Real(1) );
    AHat.MakeConsistent();

    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    auto b0 = bHat( IR(0,m), IR(0,1) );
    b0 = b;

    // G := | -I  0 0  0 |
    //      |  0 -I 0  0 |
    //      |  0  0 0  I |
    //      |  0  0 0 -I |
    // ===================
    Zeros( G, 4*n, 3*n+m );
    G.Reserve( 4*n );
    for( Int i=0; i<4*n; ++i )
    {
        if( i < 2*n )
            G.QueueUpdate( i, i,       Real(-1) );
        else if( i < 3*n )
            G.QueueUpdate( i, i+m,     Real(+1) );
        else
            G.QueueUpdate( i, i+(m-n), Real(-1) );
    }
    G.MakeConsistent();

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto ht = h( IR(2*n,4*n), IR(0,1) );
    Fill( ht, lambda );

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, IR(0,1) );
    Axpy( Real(-1), xHat(vInd,IR(0,1)), x );
}

template<typename Real>
void Var2
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        Real lambda,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ds::Var2"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numLocalEntriesA = A.NumLocalEntries();
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size(comm);
    DistSparseMatrix<Real> AHat(comm), G(comm);
    DistMultiVec<Real> c(comm), bHat(comm), h(comm);

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n+m, 1 );
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
        if( c.GlobalRow(iLoc) < 2*n )
            c.SetLocal( iLoc, 0, Real(1) );

    // G := | -I  0 0  0 |
    //      |  0 -I 0  0 |
    //      |  0  0 0  I |
    //      |  0  0 0 -I |
    // ===================
    Zeros( G, 4*n, 3*n+m );
    G.Reserve( G.LocalHeight() );
    for( Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
    {
        const Int i = G.GlobalRow(iLoc);
        if( i < 2*n )
            G.QueueLocalUpdate( iLoc, i,       Real(-1) );
        else if( i < 3*n )
            G.QueueLocalUpdate( iLoc, i+m,     Real(+1) );
        else
            G.QueueLocalUpdate( iLoc, i+(m-n), Real(-1) );
    }
    G.MakeConsistent();

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    for( Int iLoc=0; iLoc<h.LocalHeight(); ++iLoc )
        if( h.GlobalRow(iLoc) >= 2*n )
            h.SetLocal( iLoc, 0, lambda );

    // \hat A := | A, -A,  I,  0 |
    //           | 0,  0, A^T, I |
    // ===========================
    Zeros( AHat, m+n, 3*n+m );
    {
        // Compute metadata
        // ----------------
        vector<int> sendCounts(commSize,0);
        for( Int e=0; e<numLocalEntriesA; ++e )
        {
            const Int i = A.Row(e);
            const Int j = A.Col(e);
            // Sending A
            ++sendCounts[ AHat.RowOwner(i) ];
            // Sending -A (yes, I know this is technically redundant)
            ++sendCounts[ AHat.RowOwner(i) ];
            // Sending A^T
            ++sendCounts[ AHat.RowOwner(j+m) ];
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets(commSize), recvOffsets(commSize);
        const int totalSend = Scan( sendCounts, sendOffsets );
        const int totalRecv = Scan( recvCounts, recvOffsets );
        // Pack
        // ----
        vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
        vector<Real> vSendBuf(totalSend);
        auto offsets = sendOffsets;
        for( Int e=0; e<numLocalEntriesA; ++e )
        {
            const Int i = A.Row(e); 
            const Int j = A.Col(e);
            const Real value = A.Value(e);
            // Sending A
            int owner = AHat.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            tSendBuf[offsets[owner]] = j;
            vSendBuf[offsets[owner]] = value;
            ++offsets[owner];
            // Sending -A
            sSendBuf[offsets[owner]] = i;
            tSendBuf[offsets[owner]] = j + n;
            vSendBuf[offsets[owner]] = -value;
            ++offsets[owner];
            // Sending A^T
            owner = AHat.RowOwner(j+m);
            sSendBuf[offsets[owner]] = j+m;
            tSendBuf[offsets[owner]] = i+2*n;
            vSendBuf[offsets[owner]] = value;
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
        // Unpack and add in I's
        // ---------------------
        AHat.Reserve( totalRecv + AHat.LocalHeight() );
        for( Int e=0; e<totalRecv; ++e )
            AHat.QueueLocalUpdate
            ( sRecvBuf[e]-AHat.FirstLocalRow(), tRecvBuf[e], vRecvBuf[e] ); 
        for( Int iLoc=0; iLoc<AHat.LocalHeight(); ++iLoc )
        {
            const Int i = AHat.GlobalRow(iLoc);
            AHat.QueueLocalUpdate( iLoc, i+2*n, Real(1) );
        }
        AHat.MakeConsistent();
    }
    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    {
        // Compute metadata
        // ----------------
        vector<int> sendCounts(commSize,0);
        for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
        {
            const Int i = b.GlobalRow(iLoc);
            ++sendCounts[ bHat.RowOwner(i) ];
        }
        vector<int> recvCounts(commSize);
        mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
        vector<int> sendOffsets(commSize), recvOffsets(commSize);
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
            const int owner = bHat.RowOwner(i);
            sSendBuf[offsets[owner]] = i;
            vSendBuf[offsets[owner]] = b.GetLocal(iLoc,0);
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
            bHat.SetLocal( sRecvBuf[e]-bHat.FirstLocalRow(), 0, vRecvBuf[e] );
    }

    // Solve the affine LP
    // ===================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm), s(comm);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    Zeros( x, n, 1 );
    {
        // Compute metadata
        // ----------------
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
        vector<int> sendOffsets(commSize), recvOffsets(commSize);
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
            x.UpdateLocal( sRecvBuf[e]-x.FirstLocalRow(), 0, vRecvBuf[e] );
    }
}

} // namespace ds

template<typename Real>
void DS
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("DS"))
    ds::Var1( A, b, lambda, x, ctrl );
}

template<typename Real>
void DS
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("DS"))
    ds::Var1( A, b, lambda, x, ctrl );
}

template<typename Real>
void DS
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("DS"))
    ds::Var2( A, b, lambda, x, ctrl );
}

template<typename Real>
void DS
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        Real lambda,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("DS"))
    ds::Var2( A, b, lambda, x, ctrl );
}

#define PROTO(Real) \
  template void DS \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void DS \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          Real lambda, \
          AbstractDistMatrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void DS \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void DS \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          Real lambda, \
          DistMultiVec<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
