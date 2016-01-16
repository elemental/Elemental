/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace nnls {

// Solve each problem 
//
//   min || A x - b ||_2 
//   s.t. x >= 0
//
// by transforming it into the SOCP
//
//   min t
//   s.t. || A x - b ||_2 <= t, x >= 0.
//
// In terms of the canonical affine SOCP,
//
//   min c^T xHat
//   s.t. AHat xHat = bHat, G xHat + s = h, s in K,
//
// we have that
//
//   xHat = [t; x],
//
//   c = [1; 0],
//
//   AHat = [], bHat = [], and
//
//   G = | -1  0 |, h = | 0 |.
//       |  0  A |      | b |
//       |  0 -I |      | 0 |
//
// Thus, the only input data that changes for each right-hand side is h.

template<typename Real>
void SOCP
( const Matrix<Real>& A, 
  const Matrix<Real>& B, 
        Matrix<Real>& X, 
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("nnls::SOCP"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();

    Matrix<Int> orders, firstInds;
    Zeros( orders, m+n+1, 1 );
    Zeros( firstInds, m+n+1, 1 );
    for( Int i=0; i<m+1; ++i )
    {
        orders.Set( i, 0, m+1 );
        firstInds.Set( i, 0, 0 );
    }
    for( Int i=0; i<n; ++i )
    {
        orders.Set( i+m+1, 0, 1 ); 
        firstInds.Set( i+m+1, 0, i+m+1 );
    }

    // G := | -1  0 |
    //      |  0  A |
    //      |  0 -I |
    Matrix<Real> G;
    {
        Zeros( G, m+n+1, n+1 );
        G.Set( 0, 0, -1 );
        auto GA = G( IR(1,m+1), IR(1,n+1) );
        GA = A;
        auto GI = G( IR(m+1,m+n+1), IR(1,n+1) );
        Identity( GI, n, n );
        GI *= -1;
    }

    // c := [1; 0]
    Matrix<Real> c;
    Zeros( c, n+1, 1 );
    c.Set( 0, 0, 1 );

    Matrix<Real> AHat, bHat;
    Zeros( AHat, 0, n+1 );
    Zeros( bHat, 0, 1 );

    Matrix<Real> h;
    Zeros( h, m+n+1, 1 ); 

    Zeros( X, n, k );
    Matrix<Real> xHat, y, z, s;
    for( Int j=0; j<k; ++j )
    {
        auto hb = h( IR(1,m+1), ALL );
        hb = B( ALL, IR(j) );

        El::SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );

        auto x = X( ALL, IR(j) );
        x = xHat( IR(1,END), ALL );
    }
}

template<typename Real>
void SOCP
( const ElementalMatrix<Real>& APre, 
  const ElementalMatrix<Real>& BPre, 
        ElementalMatrix<Real>& XPre,
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("nnls::SOCP"))

    DistMatrixReadProxy<Real,Real,MC,MR>
      AProx( APre ),
      BProx( BPre );
    DistMatrixWriteProxy<Real,Real,MC,MR>
      XProx( XPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& X = XProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    const Grid& g = A.Grid();

    DistMatrix<Int,VC,STAR> orders(g), firstInds(g);
    Zeros( orders, m+n+1, 1 );
    Zeros( firstInds, m+n+1, 1 );
    {
        const Int localHeight = orders.LocalHeight();    
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = orders.GlobalRow(iLoc);
            if( i < m+1 )
            {
                orders.SetLocal( iLoc, 0, m+1 );
                firstInds.SetLocal( iLoc, 0, 0 );
            }
            else
            {
                orders.SetLocal( iLoc, 0, 1 );
                firstInds.SetLocal( iLoc, 0, i );
            }
        }
    }

    // G := | -1  0 |
    //      |  0  A |
    //      |  0 -I |
    DistMatrix<Real> G(g);
    {
        Zeros( G, m+n+1, n+1 );
        G.Set( 0, 0, -1 );
        auto GA = G( IR(1,m+1), IR(1,n+1) );
        GA = A;
        auto GI = G( IR(m+1,m+n+1), IR(1,n+1) );
        Identity( GI, n, n );
        GI *= -1;
    }

    // c := [1; 0]
    DistMatrix<Real> c(g);
    Zeros( c, n+1, 1 );
    c.Set( 0, 0, 1 );

    DistMatrix<Real> AHat(g), bHat(g);
    Zeros( AHat, 0, n+1 );
    Zeros( bHat, 0, 1 );

    DistMatrix<Real> h(g);
    Zeros( h, m+n+1, 1 );    

    Zeros( X, n, k );
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    for( Int j=0; j<k; ++j )
    {
        auto hb = h( IR(1,m+1), ALL );
        hb = B( ALL, IR(j) );

        El::SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );
        
        auto x = X( ALL, IR(j) );
        x = xHat( IR(1,END), ALL );
    }
}

template<typename Real>
void SOCP
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& B, 
        Matrix<Real>& X, 
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("nnls::SOCP"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    const Int numEntriesA = A.NumEntries();

    Matrix<Int> orders, firstInds;
    Zeros( orders, m+n+1, 1 );
    Zeros( firstInds, m+n+1, 1 );
    for( Int i=0; i<m+1; ++i )
    {
        orders.Set( i, 0, m+1 );
        firstInds.Set( i, 0, 0 );
    }
    for( Int i=0; i<n; ++i )
    {
        orders.Set( i+m+1, 0, 1 );
        firstInds.Set( i+m+1, 0, i+m+1 );
    }

    // G := | -1  0 |
    //      |  0  A |
    //      |  0 -I |
    SparseMatrix<Real> G;
    {
        Zeros( G, m+n+1, n+1 );
        G.Reserve( n+1+numEntriesA );
        G.QueueUpdate( 0, 0, -1 );
        for( Int e=0; e<numEntriesA; ++e )
            G.QueueUpdate( A.Row(e)+1, A.Col(e)+1, A.Value(e) );
        for( Int j=0; j<n; ++j )
            G.QueueUpdate( j+m+1, j+1, Real(-1) );
        G.ProcessQueues();
    }

    // c := [1; 0]
    Matrix<Real> c;
    Zeros( c, n+1, 1 );
    c.Set( 0, 0, 1 );

    SparseMatrix<Real> AHat;
    Zeros( AHat, 0, n+1 );
    Matrix<Real> bHat;
    Zeros( bHat, 0, 1 );

    Matrix<Real> h;
    Zeros( h, m+n+1, 1 );

    Zeros( X, n, k ); 
    Matrix<Real> xHat, y, z, s;
    for( Int j=0; j<k; ++j )
    {
        auto hb = h( IR(1,m+1), ALL );
        hb = B( ALL, IR(j) );

        El::SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );

        auto x = X( ALL, IR(j) );
        x = xHat( IR(1,END), ALL );
    }
}

template<typename Real>
void SOCP
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& B, 
        DistMultiVec<Real>& X, 
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("nnls::SOCP"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    const Int numEntriesA = A.NumLocalEntries();
    mpi::Comm comm = A.Comm();

    DistMultiVec<Int> orders(comm), firstInds(comm);
    Zeros( orders, m+n+1, 1 );
    Zeros( firstInds, m+n+1, 1 );
    {
        const Int localHeight = orders.LocalHeight();    
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = orders.GlobalRow(iLoc);
            if( i < m+1 )
            {
                orders.SetLocal( iLoc, 0, m+1 );
                firstInds.SetLocal( iLoc, 0, 0 );
            }
            else
            {
                orders.SetLocal( iLoc, 0, 1 );
                firstInds.SetLocal( iLoc, 0, i );
            }
        }
    }

    // G := | -1  0 |
    //      |  0  A |
    //      |  0 -I |
    DistSparseMatrix<Real> G(comm);
    {
        Zeros( G, m+n+1, n+1 );

        // Count the number of entries of G to reserve
        // -------------------------------------------
        Int numLocalUpdates = 0;
        const Int localHeight = G.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = G.GlobalRow(iLoc);
            if( i == 0 || i > m )
                ++numLocalUpdates;
        }
        G.Reserve( numLocalUpdates+numEntriesA, numEntriesA );

        // Queue the local updates
        // -----------------------
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = G.GlobalRow(iLoc);
            if( i == 0 )
                G.QueueLocalUpdate( iLoc, 0, -1 );
            else if( i > m )
                G.QueueLocalUpdate( iLoc, i-m, -1 );
        }

        // Queue the remote updates
        // ------------------------
        for( Int e=0; e<numEntriesA; ++e )
            G.QueueUpdate( A.Row(e)+1, A.Col(e)+1, A.Value(e) );

        G.ProcessQueues();
    }

    // c := [1; 0]
    DistMultiVec<Real> c(comm);
    Zeros( c, n+1, 1 );
    c.Set( 0, 0, 1 );

    DistSparseMatrix<Real> AHat(comm);
    Zeros( AHat, 0, n+1 );

    DistMultiVec<Real> bHat(comm);
    Zeros( bHat, 0, 1 );

    DistMultiVec<Real> h(comm);
    Zeros( h, m+n+1, 1 );

    X.SetComm( A.Comm() );
    Zeros( X, n, k ); 
    DistMultiVec<Real> x(comm), xHat(comm), y(comm), z(comm), s(comm);
    for( Int j=0; j<k; ++j )
    {
        Zeros( h, m+n+1, 1 );
        const Int bLocalHeight = B.LocalHeight();
        h.Reserve( bLocalHeight );
        for( Int iLoc=0; iLoc<bLocalHeight; ++iLoc )
            h.QueueUpdate( B.GlobalRow(iLoc)+1, 0, B.GetLocal(iLoc,j) );
        h.ProcessQueues();

        El::SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );

        Zeros( x, n, 1 );
        const Int xHatLocalHeight = xHat.LocalHeight();
        X.Reserve( xHatLocalHeight );
        for( Int iLoc=0; iLoc<xHatLocalHeight; ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            if( i > 0 )
                X.QueueUpdate( i-1, j, xHat.GetLocal(iLoc,0) );
        }
        X.ProcessQueues();
    }
}

} // namespace nnls
} // namespace El
