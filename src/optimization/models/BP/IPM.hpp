/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

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
( const Matrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
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
    AHatv -= A;

    // Solve the direct LP
    // ===================
    Matrix<Real> xHat, y, z;
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void LPIPM
( const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& b, 
        ElementalMatrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
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
    AHatv -= A;

    // Solve the direct LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g);
    LP( AHat, b, c, xHat, y, z, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void LPIPM
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
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
    x -= xHat( vInd, ALL );
}

template<typename Real>
void LPIPM
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
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
// The real Basis Pursuit (BP) problem
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
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds;
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    for( Int i=0; i<orders.Height(); ++i )
    {
        orders.Set( i, 0, 2 );
        firstInds.Set( i, 0, i-(i%2) ); 
    }

    Matrix<Real> c;
    soc::Identity( c, orders, firstInds );

    // \hat A := A E
    // =============
    Matrix<Real> AHat;
    Zeros( AHat, m, 2*n ); 
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            AHat(i,2*j+1) = A(i,j);

    // Solve the direct SOCP
    // =====================
    Matrix<Real> xHat, y, z;
    SOCP( AHat, b, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    for( Int i=0; i<xHat.Height(); ++i )
        if( i % 2 == 1 )
            x((i-1)/2) = xHat(i);
}

template<typename Real>
void SOCPIPM
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds;
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    for( Int i=0; i<orders.Height(); ++i )
    {
        orders(i) = 2;
        firstInds(i) = i-(i%2); 
    }

    Matrix<Real> c;
    soc::Identity( c, orders, firstInds );

    // \hat A := A E
    // =============
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
    SOCP( AHat, b, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    for( Int i=0; i<xHat.Height(); ++i )
        if( i % 2 == 1 )
            x((i-1)/2) = xHat(i);
}

template<typename Real>
void SOCPIPM
( const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& b, 
        ElementalMatrix<Real>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();

    DistMatrix<Int,VC,STAR> orders(grid), firstInds(grid);
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    auto& ordersLoc = orders.Matrix();
    auto& firstIndsLoc = firstInds.Matrix();
    for( Int iLoc=0; iLoc<orders.LocalHeight(); ++iLoc )
    {
        const Int i = orders.GlobalRow(iLoc);
        ordersLoc(iLoc) = 2;
        firstIndsLoc(iLoc) = i-(i%2); 
    }

    DistMatrix<Real> c(grid);
    soc::Identity( c, orders, firstInds );

    // \hat A := A E
    DistMatrix<Real> AHat(grid);
    Zeros( AHat, m, 2*n ); 
    auto& ALoc = A.LockedMatrix();
    if( A.RedundantRank() == 0 )
    {
        AHat.Reserve( A.LocalHeight()*A.LocalWidth() );
        for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
            for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
                AHat.QueueUpdate
                ( A.GlobalRow(iLoc), 2*A.GlobalCol(jLoc)+1, ALoc(iLoc,jLoc) );
    }
    AHat.ProcessQueues();

    // Solve the direct SOCP
    // =====================
    DistMatrix<Real> xHat(grid), y(grid), z(grid);
    SOCP( AHat, b, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    x.Reserve( xHat.LocalHeight() );
    auto& xHatLoc = xHat.LockedMatrix();
    if( xHat.IsLocalCol(0) )
    {
        for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            if( i % 2 == 1 )
                x.QueueUpdate( (i-1)/2, 0, xHatLoc(iLoc) );
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
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();

    DistMultiVec<Int> orders(comm), firstInds(comm);
    Zeros( orders,    2*n, 1 );
    Zeros( firstInds, 2*n, 1 );
    auto& ordersLoc = orders.Matrix();
    auto& firstIndsLoc = firstInds.Matrix();
    for( Int iLoc=0; iLoc<orders.LocalHeight(); ++iLoc )
    {
        const Int i = orders.GlobalRow(iLoc);
        ordersLoc(iLoc) = 2;
        firstIndsLoc(iLoc) = i-(i%2);
    }

    DistMultiVec<Real> c(comm);
    soc::Identity( c, orders, firstInds );
    
    // \hat A := A E
    // =============
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
    SOCP( AHat, b, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E xHat
    // ===========
    Zeros( x, n, 1 );
    x.Reserve( xHat.LocalHeight() );
    auto& xHatLoc = xHat.LockedMatrix();
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i % 2 == 1 )
            x.QueueUpdate( (i-1)/2, 0, xHatLoc(iLoc) );
    }
    x.ProcessQueues();
}

//
// The complex Basis Pursuit (BP) problem
//
//   min_x || x ||_1 s.t. A x = b
//
// can be reformulated as the Second-Order 
//
//   min_x e^T xHat 
//   s.t. | Real(A) E_R, -Imag(A) E_I | xHat = | Real(b) |, xHat in K,
//        | Imag(A) E_R,  Real(A) E_I |        | Imag(b) |
//
// where K is a product of n second-order cones of dimension 3. The operator
// E_R extracts the real components by selecting every integer which is 
// equal to 1 modulo 3, while E_I extracts the imaginary components by 
// selecting indices which are equal to 2 modulo 3.
//

template<typename Real>
void SOCPIPM
( const Matrix<Complex<Real>>& A, 
  const Matrix<Complex<Real>>& b, 
        Matrix<Complex<Real>>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds;
    Zeros( orders,    3*n, 1 );
    Zeros( firstInds, 3*n, 1 );
    for( Int i=0; i<orders.Height(); ++i )
    {
        orders(i) = 3;
        firstInds(i) = i-(i%3);
    }

    Matrix<Real> c;
    soc::Identity( c, orders, firstInds );

    // \hat A := |  Real(A) E_R - Imag(A) E_I |
    //           |  Imag(A) E_R + Real(A) E_I |
    // ========================================
    Matrix<Real> AHat;
    Zeros( AHat, 2*m, 3*n ); 
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real alphaReal = RealPart(A(i,j));
            const Real alphaImag = ImagPart(A(i,j));
            AHat(i,  3*j+1) =  alphaReal;
            AHat(i,  3*j+2) = -alphaImag;
            AHat(i+m,3*j+1) =  alphaImag;
            AHat(i+m,3*j+2) =  alphaReal;
        }
    }

    // \hat b := | Real(b) |
    //           | Imag(b) |
    // =====================
    Matrix<Real> bHat;
    Zeros( bHat, 2*m, 1 );
    for( Int i=0; i<m; ++i )
    {
        bHat(i) = RealPart(b(i));
        bHat(i+m) = ImagPart(b(i));
    }

    // Solve the direct SOCP
    // =====================
    Matrix<Real> xHat, y, z;
    SOCP( AHat, bHat, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E_R xHat + i E_I xHat
    // ==========================
    Zeros( x, n, 1 );
    for( Int i=0; i<xHat.Height(); ++i )
    {
        if( i % 3 == 1 )
            x.UpdateRealPart( (i-1)/3, 0, xHat(i) );
        else if( i % 3 == 2 )
            x.UpdateImagPart( (i-2)/3, 0, xHat(i) );
    }
}

template<typename Real>
void SOCPIPM
( const SparseMatrix<Complex<Real>>& A, 
  const Matrix<Complex<Real>>& b, 
        Matrix<Complex<Real>>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Int> orders, firstInds;
    Zeros( orders,    3*n, 1 );
    Zeros( firstInds, 3*n, 1 );
    for( Int i=0; i<orders.Height(); ++i )
    {
        orders(i) = 3;
        firstInds(i) = i-(i%3);
    }

    Matrix<Real> c;
    soc::Identity( c, orders, firstInds );

    // \hat A := |  Real(A) E_R - Imag(A) E_I |
    //           |  Imag(A) E_R + Real(A) E_I |
    // ========================================
    SparseMatrix<Real> AHat;
    Zeros( AHat, 2*m, 3*n ); 
    const Int numEntriesA = A.NumEntries();
    AHat.Reserve( 4*numEntriesA );
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Real alphaReal = RealPart(A.Value(e));
        const Real alphaImag = ImagPart(A.Value(e));
        AHat.QueueUpdate( i,   3*j+1,  alphaReal );
        AHat.QueueUpdate( i,   3*j+2, -alphaImag );
        AHat.QueueUpdate( i+m, 3*j+1,  alphaImag );
        AHat.QueueUpdate( i+m, 3*j+2,  alphaReal );
    }
    AHat.ProcessQueues();

    // \hat b := | Real(b) |
    //           | Imag(b) |
    // =====================
    Matrix<Real> bHat;
    Zeros( bHat, 2*m, 1 );
    for( Int i=0; i<m; ++i )
    {
        bHat(i) = RealPart(b(i));
        bHat(i+m) = ImagPart(b(i));
    }

    // Solve the direct SOCP
    // =====================
    Matrix<Real> xHat, y, z;
    SOCP( AHat, bHat, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E_R xHat + i E_I xHat
    // ==========================
    Zeros( x, n, 1 );
    for( Int i=0; i<xHat.Height(); ++i )
    {
        if( i % 3 == 1 )
            x.UpdateRealPart( (i-1)/3, 0, xHat(i) );
        else if( i % 3 == 2 )
            x.UpdateImagPart( (i-2)/3, 0, xHat(i) );
    }
}

template<typename Real>
void SOCPIPM
( const ElementalMatrix<Complex<Real>>& A, 
  const ElementalMatrix<Complex<Real>>& b, 
        ElementalMatrix<Complex<Real>>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();

    DistMatrix<Int,VC,STAR> orders(grid), firstInds(grid);
    Zeros( orders,    3*n, 1 );
    Zeros( firstInds, 3*n, 1 );
    auto& ordersLoc = orders.Matrix();
    auto& firstIndsLoc = firstInds.Matrix();
    for( Int iLoc=0; iLoc<orders.LocalHeight(); ++iLoc )
    {
        const Int i = orders.GlobalRow(iLoc);
        ordersLoc(iLoc) = 3;
        firstIndsLoc(iLoc) = i-(i%3);
    }

    DistMatrix<Real> c(grid);
    soc::Identity( c, orders, firstInds );

    // \hat A := |  Real(A) E_R - Imag(A) E_I |
    //           |  Imag(A) E_R + Real(A) E_I |
    // ========================================
    DistMatrix<Real> AHat(grid);
    Zeros( AHat, 2*m, 3*n ); 
    auto& ALoc = A.LockedMatrix();
    const Int ALocHeight = A.LocalHeight();
    const Int ALocWidth = A.LocalWidth();
    AHat.Reserve( 4*ALocHeight*ALocWidth );
    for( Int jLoc=0; jLoc<ALocWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<ALocHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const Real alphaReal = RealPart(ALoc(iLoc,jLoc));
            const Real alphaImag = ImagPart(ALoc(iLoc,jLoc));
            AHat.QueueUpdate( i,   3*j+1,  alphaReal );
            AHat.QueueUpdate( i,   3*j+2, -alphaImag );
            AHat.QueueUpdate( i+m, 3*j+1,  alphaImag );
            AHat.QueueUpdate( i+m, 3*j+2,  alphaReal );
        }
    }
    AHat.ProcessQueues();

    // \hat b := | Real(b) |
    //           | Imag(b) |
    // =====================
    DistMatrix<Real> bHat(grid);
    Zeros( bHat, 2*m, 1 );
    bHat.Reserve( 2*b.LocalHeight() );
    auto& bLoc = b.LockedMatrix();
    for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
    {
        const Int i = b.GlobalRow(iLoc);
        bHat.QueueUpdate( i,   0, RealPart(bLoc(i)) );
        bHat.QueueUpdate( i+m, 0, ImagPart(bLoc(i)) );
    }
    bHat.ProcessQueues();

    // Solve the direct SOCP
    // =====================
    DistMatrix<Real> xHat(grid), y(grid), z(grid);
    SOCP( AHat, bHat, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E_R xHat + i E_I xHat
    // ==========================
    Zeros( x, n, 1 );
    x.Reserve( xHat.LocalHeight() );
    if( xHat.IsLocalCol(0) )
    {
        for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        {
            const Int i = xHat.GlobalRow(iLoc);
            const Real val = xHat.GetLocal(iLoc,0);
            if( i % 3 == 1 )
                x.QueueUpdate( (i-1)/3, 0, Complex<Real>(val,0) );
            else if( i % 3 == 2 )
                x.QueueUpdate( (i-2)/3, 0, Complex<Real>(0,val) );
        }
    }
    x.ProcessQueues();
}

template<typename Real>
void SOCPIPM
( const DistSparseMatrix<Complex<Real>>& A, 
  const DistMultiVec<Complex<Real>>& b, 
        DistMultiVec<Complex<Real>>& x,
  const socp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();

    DistMultiVec<Int> orders(comm), firstInds(comm);
    Zeros( orders,    3*n, 1 );
    Zeros( firstInds, 3*n, 1 );
    auto& ordersLoc = orders.Matrix();
    auto& firstIndsLoc = firstInds.Matrix();
    for( Int iLoc=0; iLoc<orders.LocalHeight(); ++iLoc )
    {
        const Int i = orders.GlobalRow(iLoc);
        ordersLoc(iLoc) = 3;
        firstIndsLoc(iLoc) = i-(i%3);
    }

    DistMultiVec<Real> c(comm);
    soc::Identity( c, orders, firstInds );
    
    // \hat A := |  Real(A) E_R - Imag(A) E_I |
    //           |  Imag(A) E_R + Real(A) E_I |
    // ========================================
    DistSparseMatrix<Real> AHat(comm);
    const Int numLocalEntriesA = A.NumLocalEntries();
    Zeros( AHat, 2*m, 3*n ); 
    AHat.Reserve( 4*numLocalEntriesA, 4*numLocalEntriesA );
    for( Int e=0; e<numLocalEntriesA; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Real alphaReal = RealPart(A.Value(e));
        const Real alphaImag = ImagPart(A.Value(e));
        AHat.QueueUpdate( i,   3*j+1,  alphaReal );
        AHat.QueueUpdate( i,   3*j+2, -alphaImag );
        AHat.QueueUpdate( i+m, 3*j+1,  alphaImag );
        AHat.QueueUpdate( i+m, 3*j+2,  alphaReal );
    }
    AHat.ProcessQueues();

    // \hat b := | Real(b) |
    //           | Imag(b) |
    // =====================
    DistMultiVec<Real> bHat(comm);
    Zeros( bHat, 2*m, 1 );
    bHat.Reserve( 2*b.LocalHeight() );
    auto& bLoc = b.LockedMatrix();
    for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
    {
        const Int i = b.GlobalRow(iLoc);
        bHat.QueueUpdate( i,   0, RealPart(bLoc(iLoc)) );
        bHat.QueueUpdate( i+m, 0, ImagPart(bLoc(iLoc)) );
    }
    bHat.ProcessQueues();

    // Solve the direct SOCP
    // =====================
    DistMultiVec<Real> xHat(comm), y(comm), z(comm);
    SOCP( AHat, bHat, c, orders, firstInds, xHat, y, z, ctrl );

    // x := E_R xHat + i E_I xHat
    // ==========================
    Zeros( x, n, 1 );
    x.Reserve( xHat.LocalHeight() );
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        const Real val = xHat.GetLocal(iLoc,0);
        if( i % 3 == 1 )
            x.QueueUpdate( (i-1)/3, 0, Complex<Real>(val,0) );
        else if( i % 3 == 2 )
            x.QueueUpdate( (i-2)/3, 0, Complex<Real>(0,val) );
    }
    x.ProcessQueues();
}

} // namespace bp
} // namespace El
