/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

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
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), tInd(2*n,3*n);
    Matrix<Real> c, AHat, bHat, G, h;

    // c := [1;1;0]
    // ==============
    Zeros( c, 3*n, 1 );
    auto cuv = c( IR(0,2*n), ALL );
    Fill( cuv, Real(1) );

    // \hat A := [ A^T A, -A^T A,  I ]
    // ===============================
    Zeros( AHat, n, 3*n );
    auto AHatu = AHat( ALL, uInd );
    auto AHatv = AHat( ALL, vInd );
    auto AHatt = AHat( ALL, tInd );
    Herk( LOWER, TRANSPOSE, Real(1), A, Real(0), AHatu );
    MakeSymmetric( LOWER, AHatu );
    AHatv -= AHatu;
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
    auto ht = h( IR(2*n,4*n), ALL );
    Fill( ht, lambda );

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void Var1
( const AbstractDistMatrix<Real>& APre,
  const AbstractDistMatrix<Real>& b,
        Real lambda,
        AbstractDistMatrix<Real>& xPre,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Real,Real,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    DistMatrixWriteProxy<Real,Real,MC,MR> xProx( xPre );
    auto& x = xProx.Get();

    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> uInd(0,n), vInd(n,2*n), tInd(2*n,3*n);
    DistMatrix<Real> c(g), AHat(g), bHat(g), G(g), h(g);

    // c := [1;1;0]
    // ==============
    Zeros( c, 3*n, 1 );
    auto cuv = c( IR(0,2*n), ALL );
    Fill( cuv, Real(1) );

    // \hat A := [ A^T A, -A^T A,  I ]
    // ===============================
    Zeros( AHat, n, 3*n );
    auto AHatu = AHat( ALL, uInd );
    auto AHatv = AHat( ALL, vInd );
    auto AHatt = AHat( ALL, tInd );
    Herk( LOWER, TRANSPOSE, Real(1), A, Real(0), AHatu );
    MakeSymmetric( LOWER, AHatu );
    AHatv -= AHatu;
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
    auto ht = h( IR(2*n,4*n), ALL );
    Fill( ht, lambda );

    // Solve the affine LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void Var2
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m), tInd(2*n+m,3*n+m);
    Matrix<Real> c, AHat, bHat, G, h;

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n+m, 1 );
    auto c_u = c( uInd, ALL );
    auto c_v = c( vInd, ALL );
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
    AHat_0v = A; AHat_0v *= -1;
    Identity( AHat_0r, m, m );
    Transpose( A, AHat_1r );
    Identity( AHat_1t, n, n );

    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    auto b0 = bHat( IR(0,m), ALL );
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
    Identity( G0u, n, n ); G0u *= -1;
    Identity( G1v, n, n ); G1v *= -1;
    Identity( G2t, n, n );
    Identity( G3t, n, n ); G3t *= -1;

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto h2 = h( IR(2*n,3*n), ALL );
    auto h3 = h( IR(3*n,4*n), ALL );
    Fill( h2, lambda );
    Fill( h3, lambda );

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void Var2
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
        Real lambda,
        AbstractDistMatrix<Real>& xPre,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixWriteProxy<Real,Real,MC,MR> xProx( xPre );
    auto& x = xProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m), tInd(2*n+m,3*n+m);
    DistMatrix<Real> c(g), AHat(g), bHat(g), G(g), h(g);

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n+m, 1 );
    auto c_u = c( uInd, ALL );
    auto c_v = c( vInd, ALL );
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
    AHat_0v = A; AHat_0v *= -1;
    Identity( AHat_0r, m, m );
    Transpose( A, AHat_1r );
    Identity( AHat_1t, n, n );

    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    auto b0 = bHat( IR(0,m), ALL );
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
    Identity( G0u, n, n ); G0u *= -1;
    Identity( G1v, n, n ); G1v *= -1;
    Identity( G2t, n, n );
    Identity( G3t, n, n ); G3t *= -1;

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto h2 = h( IR(2*n,3*n), ALL );
    auto h3 = h( IR(3*n,4*n), ALL );
    Fill( h2, lambda );
    Fill( h3, lambda );

    // Solve the affine LP
    // ===================
    DistMatrix<Real> xHat(g), y(g), z(g), s(g);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void Var2
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntriesA = A.NumEntries();
    const Range<Int> uInd(0,n), vInd(n,2*n);
    SparseMatrix<Real> AHat, G;
    Matrix<Real> c, bHat, h;

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n+m, 1 );
    auto cuv = c( IR(0,2*n), ALL );
    Fill( cuv, Real(1) );

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
    G.ProcessQueues();

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto ht = h( IR(2*n,4*n), ALL );
    Fill( ht, lambda );

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
    AHat.ProcessQueues();

    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    auto b0 = bHat( IR(0,m), ALL );
    b0 = b;

    // Solve the affine LP
    // ===================
    Matrix<Real> xHat, y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void Var2
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
        Real lambda,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numLocalEntriesA = A.NumLocalEntries();
    const Grid& grid = A.Grid();

    DistSparseMatrix<Real> AHat(grid), G(grid);
    DistMultiVec<Real> c(grid), bHat(grid), h(grid);

    auto& bLoc = b.LockedMatrix();

    // c := [1;1;0;0]
    // ==============
    Zeros( c, 3*n+m, 1 );
    auto& cLoc = c.Matrix();
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
        if( c.GlobalRow(iLoc) < 2*n )
            cLoc(iLoc) = Real(1);

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
    G.ProcessLocalQueues();

    // h := [0;0;lambda e;lambda e]
    // ============================
    Zeros( h, 4*n, 1 );
    auto& hLoc = h.Matrix();
    for( Int iLoc=0; iLoc<h.LocalHeight(); ++iLoc )
        if( h.GlobalRow(iLoc) >= 2*n )
            hLoc(iLoc) = lambda;

    // \hat A := | A, -A,  I,  0 |
    //           | 0,  0, A^T, I |
    // ===========================
    Zeros( AHat, m+n, 3*n+m );
    AHat.Reserve( 3*numLocalEntriesA+AHat.LocalHeight(), 3*numLocalEntriesA );
    for( Int e=0; e<numLocalEntriesA; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Real value = A.Value(e);
        AHat.QueueUpdate( i,   j,     value );
        AHat.QueueUpdate( i,   j+n,  -value );
        AHat.QueueUpdate( j+m, i+2*n, value );
    }
    for( Int iLoc=0; iLoc<AHat.LocalHeight(); ++iLoc )
    {
        const Int i = AHat.GlobalRow(iLoc);
        AHat.QueueUpdate( i, i+2*n, Real(1) );
    }
    AHat.ProcessQueues();

    // \hat b := | b |
    //           | 0 |
    // ===============
    Zeros( bHat, m+n, 1 );
    bHat.Reserve( b.LocalHeight() );
    for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
    {
        const Int i = b.GlobalRow(iLoc);
        bHat.QueueUpdate( i, 0, bLoc(iLoc) );
    }
    bHat.ProcessQueues();

    // Solve the affine LP
    // ===================
    DistMultiVec<Real> xHat(grid), y(grid), z(grid), s(grid);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    auto& xHatLoc = xHat.LockedMatrix();
    Zeros( x, n, 1 );
    Int numRemoteUpdates = 0;
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
        if( xHat.GlobalRow(iLoc) < 2*n )
            ++numRemoteUpdates;
        else
            break;
    x.Reserve( numRemoteUpdates );
    for( Int iLoc=0; iLoc<xHat.LocalHeight(); ++iLoc )
    {
        const Int i = xHat.GlobalRow(iLoc);
        if( i < n )
            x.QueueUpdate( i, 0, xHatLoc(iLoc) );
        else if( i < 2*n )
            x.QueueUpdate( i-n, 0, -xHatLoc(iLoc) );
        else
            break;
    }
    x.ProcessQueues();
}

} // namespace ds

// TODO(poulson): Add the ability to choose the variant (while preserving the
// Var1 dense default and Var2 sparse default).

template<typename Real>
void DS
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    ds::Var1( A, b, lambda, x, ctrl );
}

template<typename Real>
void DS
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    ds::Var1( A, b, lambda, x, ctrl );
}

template<typename Real>
void DS
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
        Real lambda,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    ds::Var2( A, b, lambda, x, ctrl );
}

template<typename Real>
void DS
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
        Real lambda,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    ds::Var2( A, b, lambda, x, ctrl );
}

#define PROTO(Real) \
  template void DS \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void DS \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
          Real lambda, \
          AbstractDistMatrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void DS \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void DS \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
          Real lambda, \
          DistMultiVec<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
