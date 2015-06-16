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
// decomposing x into its positive and negative parts, say (u,v), and posing
//
//   min 1^T [u;v]
//   s.t. [A, -A] [u; v] = b, [u; v] >= 0.
//
// After solving this LP, the solution is set to x := u - v.
//
// Complex instances of Basis Pursuit require Second-Order Cone Programming.
//
// [1] Scott S. Chen, David L. Donoho, and Michael A. Saunders,
//     "Atomic Decomposition by Basis Pursuit",
//     SIAM Review, Vol. 43, No. 1, pp. 129--159, 2001

// TODO: Extend the existing LP control parameters
// TODO: Extend the SOCP control parameters

namespace El {
namespace bp {

template<typename Real>
void LPIPM
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::LPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n);
    Matrix<Real> c, AHat;

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    Zeros( AHat, m, 2*n );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    AHatu = A;
    Axpy( Real(-1), A, AHatv );

    // Solve the direct LP
    // ===================
    Matrix<Real> xHat, y, z;
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    Axpy( Real(-1), xHat(vInd,ALL), x );
}

template<typename Real>
void LPIPM
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::LPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> uInd(0,n), vInd(n,2*n);
    DistMatrix<Real> c(g), AHat(g);

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    Zeros( AHat, m, 2*n );
    auto AHatu = AHat( IR(0,m), uInd );
    auto AHatv = AHat( IR(0,m), vInd );
    AHatu = A;
    Axpy( Real(-1), A, AHatv );

    // Solve the direct LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g);
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    Copy( xHat( uInd, ALL ), x );
    Axpy( Real(-1), xHat(vInd,ALL), x );
}

template<typename Real>
void LPIPM
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::LPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n);
    SparseMatrix<Real> AHat;
    Matrix<Real> c;

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // =================
    const Int numEntriesA = A.NumEntries();
    Zeros( AHat, m, 2*n );
    AHat.Reserve( 2*numEntriesA );
    for( Int e=0; e<numEntriesA; ++e )
    {
        AHat.QueueUpdate( A.Row(e), A.Col(e),    A.Value(e) );
        AHat.QueueUpdate( A.Row(e), A.Col(e)+n, -A.Value(e) );
    }
    AHat.ProcessQueues();

    // Solve the direct LP
    // ===================
    Matrix<Real> xHat, y, z;
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    Axpy( Real(-1), xHat(vInd,ALL), x );
}

template<typename Real>
void LPIPM
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::LPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> AHat(comm);
    DistMultiVec<Real> c(comm);

    // c := ones(2*n,1)
    // ================
    Ones( c, 2*n, 1 );

    // \hat A := [A, -A]
    // ================
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    const Int numLocalEntriesA = A.NumLocalEntries();
    Zeros( AHat, m, 2*n );
    AHat.Reserve( 2*numLocalEntriesA );
    for( Int e=0; e<numLocalEntriesA; ++e )
    {
        AHat.QueueUpdate( A.Row(e), A.Col(e),    A.Value(e) );
        AHat.QueueUpdate( A.Row(e), A.Col(e)+n, -A.Value(e) );
    }
    AHat.ProcessLocalQueues();

    // Solve the direct LP
    // ===================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm);
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    Zeros( x, n, 1 );
    x.Reserve( 2*xHat.LocalHeight() );
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i < n )
            x.QueueUpdate( i,   0,  xHat.GetLocal(iLoc,0) );
        else
            x.QueueUpdate( i-n, 0, -xHat.GetLocal(iLoc,0) );
    }
    x.ProcessQueues();
}

//
// The Basis Pursuit (BP) problem
//
//   min_x || x ||_1 s.t. A x = b
//
// can be reformulated as the Second-Order 
//
//   min_x e^T xHat s.t. (A E) xHat = b, xHat in K,
//
// where 'E' is the extraction operator, which, in this case, extracts the
// odd entries of the input vector.
//

template<typename Real>
void SOCPIPM
( const Matrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::SOCPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds, labels;
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    Zeros( labels,    2*n, 1 );
    for( Int i=0; i<orders.Height(); ++i )
    {
        orders.Set( i, 0, 2 );
        firstInds.Set( i, 0, i-(i%2) ); 
        labels.Set( i, 0, i/2 );
    }

    Matrix<Real> c;
    SOCIdentity( c, orders, firstInds );

    // \hat A := A E
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    Matrix<Real> AHat;
    Zeros( AHat, m, 2*n ); 
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            AHat.Set( i, 2*j+1, A.Get(i,j) );

    // Solve the direct SOCP
    // =====================
    Matrix<Real> xHat, y, z;
    SOCP( AHat, b, c, xHat, y, z, orders, firstInds, labels, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    for( Int i=0; i<xHat.Height(); ++i )
        if( i % 2 == 1 )
            x.Set( (i-1)/2, 0, xHat.Get(i,0) );
}

template<typename Real>
void SOCPIPM
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::SOCPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds, labels;
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    Zeros( labels,    2*n, 1 );
    for( Int i=0; i<orders.Height(); ++i )
    {
        orders.Set( i, 0, 2 );
        firstInds.Set( i, 0, i-(i%2) ); 
        labels.Set( i, 0, i/2 );
    }

    Matrix<Real> c;
    SOCIdentity( c, orders, firstInds );

    // \hat A := A E
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    SparseMatrix<Real> AHat;
    Zeros( AHat, m, 2*n ); 
    const Int numEntriesA = A.NumEntries();
    AHat.Reserve( numEntriesA );
    for( Int e=0; e<numEntriesA; ++e )
        AHat.QueueUpdate( A.Row(e), 2*A.Col(e)+1, A.Value(e) );
    AHat.ProcessQueues();

    // Solve the direct SOCP
    // =====================
    Matrix<Real> xHat, y, z;
    SOCP( AHat, b, c, xHat, y, z, orders, firstInds, labels, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    for( Int i=0; i<xHat.Height(); ++i )
        if( i % 2 == 1 )
            x.Set( (i-1)/2, 0, xHat.Get(i,0) );
}

template<typename Real>
void SOCPIPM
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::SOCPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();

    DistMatrix<Int,VC,STAR> orders(grid), firstInds(grid), labels(grid);
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    Zeros( labels,    2*n, 1 );
    for( Int iLoc=0; iLoc<orders.LocalHeight(); ++iLoc )
    {
        const Int i = orders.GlobalRow(iLoc);
        orders.SetLocal( iLoc, 0, 2 );
        firstInds.SetLocal( iLoc, 0, i-(i%2) ); 
        labels.SetLocal( iLoc, 0, i/2 );
    }

    DistMatrix<Real> c(grid);
    SOCIdentity( c, orders, firstInds );

    // \hat A := A E
    DistMatrix<Real> AHat(grid);
    Zeros( AHat, m, 2*n ); 
    if( A.RedundantRank() == 0 )
    {
        AHat.Reserve( A.LocalHeight()*A.LocalWidth() );
        for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
            for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
                AHat.QueueUpdate
                ( A.GlobalRow(iLoc), 2*A.GlobalCol(jLoc)+1, 
                  A.GetLocal(iLoc,jLoc) );
    }
    AHat.ProcessQueues();

    // Solve the direct SOCP
    // =====================
    DistMatrix<Real> xHat(grid), y(grid), z(grid);
    SOCP( AHat, b, c, xHat, y, z, orders, firstInds, labels, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    x.Reserve( xHat.LocalHeight() );
    if( xHat.IsLocalCol(0) )
    {
        for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            if( i % 2 == 1 )
                x.QueueUpdate( (i-1)/2, 0, xHat.GetLocal(iLoc,0) );
        }
    }
    x.ProcessQueues();
}

template<typename Real>
void SOCPIPM
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bp::SOCPIPM"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();

    DistMultiVec<Int> orders(comm), firstInds(comm), labels(comm);
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    Zeros( labels,    2*n, 1 );
    for( Int iLoc=0; iLoc<orders.LocalHeight(); ++iLoc )
    {
        const Int i = orders.GlobalRow(iLoc);
        orders.SetLocal( iLoc, 0, 2 );
        firstInds.SetLocal( iLoc, 0, i-(i%2) ); 
        labels.SetLocal( iLoc, 0, i/2 );
    }

    DistMultiVec<Real> c(comm);
    SOCIdentity( c, orders, firstInds );
    
    // \hat A := A E
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    DistSparseMatrix<Real> AHat(comm);
    Zeros( AHat, m, 2*n ); 
    const Int numLocalEntriesA = A.NumLocalEntries();
    AHat.Reserve( numLocalEntriesA );
    for( Int e=0; e<numLocalEntriesA; ++e )
        AHat.QueueUpdate( A.Row(e), 2*A.Col(e)+1, A.Value(e) );
    AHat.ProcessLocalQueues();

    // Solve the direct SOCP
    // =====================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm);
    SOCP( AHat, b, c, xHat, y, z, orders, firstInds, labels, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    x.Reserve( xHat.LocalHeight() );
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i % 2 == 1 )
            x.QueueUpdate( (i-1)/2, 0, xHat.GetLocal(iLoc,0) );
    }
    x.ProcessQueues();
}

} // namespace bp
} // namespace El
