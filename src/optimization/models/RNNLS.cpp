/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Given || [dA, db] ||_2 <= rho, minimize the worst-case error of
//
//     || (A+dA) x - (b+db) ||_2,
//
// which can be shown to be equal to
//
//     || A x - b ||_2 + rho || [x; 1] ||_2,
//
// subject to x >= 0, which can be formulated as the SOCP
//
//     min t + rho s 
//     s.t. || A x - b ||_2 <= t, || [x; 1] ||_2 <= s, x >= 0.
//
// (See [1] or Subsection 2.7 of [2].)
//
// [1] L. El Ghaoui and H. Lebret, "Robust solutions to least-squares problems
//     with uncertain data", SIAM J. Matrix Anal. and Appl., Vol. 18, No. 4,
//     1997. DOI: http://epubs.siam.org/doi/abs/10.1137/S0895479896298130
//
// [2] M.S. Lobo, L. Vandenberghe, S. Boyd, and H. Lebret, 
//     "Applications of second-order cone programming", 
//     Linear Algebra and its Applications, Vol. 284, Issues 1-3, 1998. 
//     DOI: http://www.sciencedirect.com/science/article/pii/S0024379598100320
// 

template<typename Real>
void RNNLS
( const Matrix<Real>& A, 
  const Matrix<Real>& b,
        Real rho,
        Matrix<Real>& x, 
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("RNNLS"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds;
    Zeros( orders, m+2*n+3, 1 );
    Zeros( firstInds, m+2*n+3, 1 );
    for( Int i=0; i<m+1; ++i )
    {
        orders.Set( i, 0, m+1 );
        firstInds.Set( i, 0, 0 );
    }
    for( Int i=0; i<n+2; ++i )
    {
        orders.Set( i+m+1, 0, n+2 ); 
        firstInds.Set( i+m+1, 0, m+1 );
    }
    for( Int i=0; i<n; ++i )
    {
        orders.Set( i+m+n+3, 0, 1 );
        firstInds.Set( i+m+n+3, 0, i+m+n+3 );
    }

    // G := | -1  0  0 |
    //      |  0  0  A |
    //      |  0 -1  0 |
    //      |  0  0 -I |
    //      |  0  0  0 |
    //      |  0  0 -I |
    Matrix<Real> G;
    {
        Zeros( G, m+2*n+3, n+2 );
        G.Set( 0, 0, -1 );
        auto GA = G( IR(1,m+1), IR(2,n+2) );
        GA = A;
        G.Set( m+1, 1, -1 );
        {
            auto GI = G( IR(m+2,m+n+2), IR(2,END) );
            Identity( GI, n, n );
            GI *= -1;
        }
        {
            auto GI = G( IR(m+n+3,END), IR(2,END) );
            Identity( GI, n, n );
            GI *= -1;
        }
    }

    // h := | 0 |
    //      | b |
    //      | 0 |
    //      | 0 |
    //      | 1 |
    //      | 0 |
    Matrix<Real> h;
    Zeros( h, m+2*n+3, 1 ); 
    auto hb = h( IR(1,m+1), ALL );
    hb = b;
    h.Set( m+n+2, 0, 1 );

    // c := [1; rho; 0]
    Matrix<Real> c;
    Zeros( c, n+2, 1 );
    c.Set( 0, 0, 1 );
    c.Set( 1, 0, rho );

    Matrix<Real> AHat, bHat;
    Zeros( AHat, 0, n+2 );
    Zeros( bHat, 0, 1 );

    Matrix<Real> xHat, y, z, s;
    SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );
    x = xHat( IR(2,END), ALL );
}

template<typename Real>
void RNNLS
( const ElementalMatrix<Real>& APre, 
  const ElementalMatrix<Real>& bPre, 
        Real rho,
        ElementalMatrix<Real>& xPre,
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("RNNLS"))

    DistMatrixReadProxy<Real,Real,MC,MR>
      AProx( APre ),
      bProx( bPre );
    DistMatrixWriteProxy<Real,Real,MC,MR>
      xProx( xPre );
    auto& A = AProx.GetLocked();
    auto& b = bProx.GetLocked();
    auto& x = xProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();

    DistMatrix<Int,VC,STAR> orders(g), firstInds(g);
    Zeros( orders, m+2*n+3, 1 );
    Zeros( firstInds, m+2*n+3, 1 );
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
            else if( i < m+n+3 )
            {
                orders.SetLocal( iLoc, 0, n+2 );
                firstInds.SetLocal( iLoc, 0, m+1 );
            }
            else
            {
                orders.SetLocal( iLoc, 0, 1 );
                firstInds.SetLocal( iLoc, 0, i );
            }
        }
    }

    // G := | -1  0  0 |
    //      |  0  0  A |
    //      |  0 -1  0 |
    //      |  0  0 -I |
    //      |  0  0  0 |
    //      |  0  0 -I |
    DistMatrix<Real> G(g);
    {
        Zeros( G, m+2*n+3, n+2 );
        G.Set( 0, 0, -1 );
        auto GA = G( IR(1,m+1), IR(2,END) );
        GA = A;
        G.Set( m+1, 1, -1 );
        {
            auto GI = G( IR(m+2,m+n+2), IR(2,END) );
            Identity( GI, n, n );
            GI *= -1;
        }
        {
            auto GI = G( IR(m+n+3,END), IR(2,END) );
            Identity( GI, n, n );
            GI *= -1;
        }
    }

    // h := | 0 |
    //      | b |
    //      | 0 |
    //      | 0 |
    //      | 1 |
    //      | 0 |
    DistMatrix<Real> h(g);
    Zeros( h, m+2*n+3, 1 );
    auto hb = h( IR(1,m+1), ALL );
    hb = b;
    h.Set( m+n+2, 0, 1 );

    // c := [1; rho; 0]
    DistMatrix<Real> c(g);
    Zeros( c, n+2, 1 );
    c.Set( 0, 0, 1 );
    c.Set( 1, 0, rho );

    DistMatrix<Real> AHat(g), bHat(g);
    Zeros( AHat, 0, n+2 );
    Zeros( bHat, 0, 1 );

    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );
    x = xHat( IR(2,END), ALL );
}

template<typename Real>
void RNNLS
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b, 
        Real rho,
        Matrix<Real>& x, 
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("RNNLS"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds;
    Zeros( orders, m+2*n+3, 1 );
    Zeros( firstInds, m+2*n+3, 1 );
    for( Int i=0; i<m+1; ++i )
    {
        orders.Set( i, 0, m+1 );
        firstInds.Set( i, 0, 0 );
    }
    for( Int i=0; i<n+2; ++i )
    {
        orders.Set( i+m+1, 0, n+2 ); 
        firstInds.Set( i+m+1, 0, m+1 );
    }
    for( Int i=0; i<n; ++i )
    {
        orders.Set( i+m+n+3, 0, 1 );
        firstInds.Set( i+m+n+3, 0, i+m+n+3 );
    }

    // G := | -1  0  0 |
    //      |  0  0  A |
    //      |  0 -1  0 |
    //      |  0  0 -I |
    //      |  0  0  0 |
    //      |  0  0 -I |
    SparseMatrix<Real> G;
    {
        const Int numEntriesA = A.NumEntries();
        Zeros( G, m+2*n+3, n+2 );
        G.Reserve( numEntriesA+2*n+2 );
        G.QueueUpdate( 0, 0, -1 );
        for( Int e=0; e<numEntriesA; ++e )
            G.QueueUpdate( A.Row(e)+1, A.Col(e)+2, A.Value(e) );
        G.QueueUpdate( m+1, 1, -1 );
        for( Int j=0; j<n; ++j )
            G.QueueUpdate( j+m+2, j+2, -1 );
        for( Int j=0; j<n; ++j )
            G.QueueUpdate( j+m+n+3, j+2, -1 );
        G.ProcessQueues();
    }

    // h := | 0 |
    //      | b |
    //      | 0 |
    //      | 0 |
    //      | 1 |
    //      | 0 |
    Matrix<Real> h;
    Zeros( h, m+2*n+3, 1 ); 
    auto hb = h( IR(1,m+1), ALL );
    hb = b;
    h.Set( m+n+2, 0, 1 );

    // c := [1; rho; 0]
    Matrix<Real> c;
    Zeros( c, n+2, 1 );
    c.Set( 0, 0, 1 );
    c.Set( 1, 0, rho );

    SparseMatrix<Real> AHat;
    Zeros( AHat, 0, n+2 );
    Matrix<Real> bHat;
    Zeros( bHat, 0, 1 );

    Matrix<Real> xHat, y, z, s;
    SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );
    x = xHat( IR(2,END), ALL );
}

template<typename Real>
void RNNLS
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, 
        Real rho,
        DistMultiVec<Real>& x, 
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("RNNLS"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();

    DistMultiVec<Int> orders(comm), firstInds(comm);
    Zeros( orders, m+2*n+3, 1 );
    Zeros( firstInds, m+2*n+3, 1 );
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
            else if( i < m+n+3 )
            {
                orders.SetLocal( iLoc, 0, n+2 ); 
                firstInds.SetLocal( iLoc, 0, m+1 );
            }
            else
            {
                orders.SetLocal( iLoc, 0, 1 );
                firstInds.SetLocal( iLoc, 0, i );
            }
        }
    }

    // G := | -1  0  0 |
    //      |  0  0  A |
    //      |  0 -1  0 |
    //      |  0  0 -I |
    //      |  0  0  0 |
    //      |  0  0 -I |
    DistSparseMatrix<Real> G(comm);
    {
        Zeros( G, m+2*n+3, n+2 );

        // Count the number of entries of G to reserve
        // -------------------------------------------
        Int numLocalUpdates = 0;
        const Int localHeight = G.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = G.GlobalRow(iLoc);
            if( i == 0 || i == m+1 || (i>m+1 && i<m+n+2) || i>=m+n+3 )
                ++numLocalUpdates;
        }
        const Int numEntriesA = A.NumLocalEntries();
        G.Reserve( numLocalUpdates+numEntriesA, numEntriesA );

        // Queue the local updates
        // ----------------------- 
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = G.GlobalRow(iLoc);
            if( i == 0 )
                G.QueueLocalUpdate( iLoc, 0, -1 );
            else if( i == m+1 )
                G.QueueLocalUpdate( iLoc, 1, -1 );
            else if( i > m+1 && i < m+n+2 )
                G.QueueLocalUpdate( iLoc, i-m, -1 );
            else if( i >= m+n+3 )
                G.QueueLocalUpdate( iLoc, i-(m+n+1), -1 );
        }

        // Queue the remote updates
        // ------------------------
        for( Int e=0; e<numEntriesA; ++e )
            G.QueueUpdate( A.Row(e)+1, A.Col(e)+2, A.Value(e) );

        G.ProcessQueues();
    }

    // h := | 0 |
    //      | b |
    //      | 0 |
    //      | 0 |
    //      | 1 |
    //      | 0 |
    DistMultiVec<Real> h(comm);
    Zeros( h, m+2*n+3, 1 ); 
    {
        const Int bLocalHeight = b.LocalHeight();
        h.Reserve( bLocalHeight );
        for( Int iLoc=0; iLoc<bLocalHeight; ++iLoc )
            h.QueueUpdate( b.GlobalRow(iLoc)+1, 0, b.GetLocal(iLoc,0) );
        h.ProcessQueues();
    }
    h.Set( m+n+2, 0, 1 );

    // c := [1; rho; 0]
    DistMultiVec<Real> c(comm);
    Zeros( c, n+2, 1 );
    c.Set( 0, 0, 1 );
    c.Set( 1, 0, rho );

    DistSparseMatrix<Real> AHat(comm);
    Zeros( AHat, 0, n+2 );
    DistMultiVec<Real> bHat(comm);
    Zeros( bHat, 0, 1 );

    DistMultiVec<Real> xHat(comm), y(comm), z(comm), s(comm);
    SOCP( AHat, G, bHat, c, h, orders, firstInds, xHat, y, z, s, ctrl );
    x = xHat( IR(2,END), ALL );
}

#define PROTO(Real) \
  template void RNNLS \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
          Real rho, \
          Matrix<Real>& x, \
    const socp::affine::Ctrl<Real>& ctrl ); \
  template void RNNLS \
  ( const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& b, \
          Real rho, \
          ElementalMatrix<Real>& x, \
    const socp::affine::Ctrl<Real>& ctrl ); \
  template void RNNLS \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
          Real rho, \
          Matrix<Real>& x, \
    const socp::affine::Ctrl<Real>& ctrl ); \
  template void RNNLS \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
          Real rho, \
          DistMultiVec<Real>& x, \
    const socp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
