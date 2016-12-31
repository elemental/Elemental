/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// An elastic net seeks the solution to the optimization problem
//
//   min || b - A x ||_2^2 + lambda_1 || x ||_1.+ lambda_2 || x ||_2^2.
//
// Real instances of the problem are expressable as a Quadratic Program [1] via
// the transformation
//
//   min r^T r + lambda_1 1^T [u;v] + lambda_2 (u^T u + v^T v)
//   s.t. [A, -A] [u; v] + r = b, [u; v] >= 0.
//
// When expressed in affine conic form, the above expression becomes
//
//   min (1/2) [u;v;r]^T | 2*lambda_2     0      0 | | u | +
//                       |     0      2*lambda_2 0 | | v |
//                       |     0          0      2 | | r |
//       lambda_1 [1;1;0]^T [u;v;r]
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

namespace El {

template<typename Real>
void EN
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Real lambda1,
        Real lambda2,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m);

    Matrix<Real> Q, c, AHat, G, h;

    // Q := | 2*lambda_2     0      0 |
    //      |     0      2*lambda_2 0 |
    //      |     0          0      2 |
    // ================================
    Zeros( Q, 2*n+m, 2*n+m );
    auto QTL = Q( IR(0,2*n), IR(0,2*n) );
    FillDiagonal( QTL, 2*lambda2 );
    auto Qrr = Q( rInd, rInd );
    FillDiagonal( Qrr, Real(1) );

    // c := lambda_1*[1;1;0]
    // =====================
    Zeros( c, 2*n+m, 1 );
    auto cuv = c( IR(0,2*n), ALL );
    Fill( cuv, lambda1 );

    // \hat A := [A, -A, I]
    // ====================
    Zeros( AHat, m, 2*n+m );
    auto AHatu = AHat( ALL, uInd );
    auto AHatv = AHat( ALL, vInd );
    auto AHatr = AHat( ALL, rInd );
    AHatu = A;
    AHatv -= A;
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
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void EN
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
        Real lambda1,
        Real lambda2,
        AbstractDistMatrix<Real>& xPre,
  const qp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixWriteProxy<Real,Real,MC,MR> xProx( xPre );
    auto& x = xProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m);
    DistMatrix<Real> Q(g), c(g), AHat(g), G(g), h(g);

    // Q := | 2*lambda_2     0      0 |
    //      |     0      2*lambda_2 0 |
    //      |     0          0      2 |
    // ================================
    Zeros( Q, 2*n+m, 2*n+m );
    auto QTL = Q( IR(0,2*n), IR(0,2*n) );
    FillDiagonal( QTL, 2*lambda2 );
    auto Qrr = Q( rInd, rInd );
    FillDiagonal( Qrr, Real(1) );

    // c := lambda_1*[1;1;0]
    // =====================
    Zeros( c, 2*n+m, 1 );
    auto cuv = c( IR(0,2*n), ALL );
    Fill( cuv, lambda1 );

    // \hat A := [A, -A, I]
    // ====================
    Zeros( AHat, m, 2*n+m );
    auto AHatu = AHat( ALL, uInd );
    auto AHatv = AHat( ALL, vInd );
    auto AHatr = AHat( ALL, rInd );
    AHatu = A;
    AHatv -= A;
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
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void EN
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
        Real lambda1,
        Real lambda2,
        Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> uInd(0,n), vInd(n,2*n), rInd(2*n,2*n+m);
    SparseMatrix<Real> Q, AHat, G;
    Matrix<Real> c, h;

    // Q := | 2*lambda_2     0      0 |
    //      |     0      2*lambda_2 0 |
    //      |     0          0      2 |
    // ================================
    Zeros( Q, 2*n+m, 2*n+m );
    Q.Reserve( 2*n+m );
    for( Int e=0; e<2*n; ++e )
        Q.QueueUpdate( e, e, 2*lambda2 );
    for( Int e=0; e<m; ++e )
        Q.QueueUpdate( 2*n+e, 2*n+e, Real(1) );
    Q.ProcessQueues();

    // c := lambda_1*[1;1;0]
    // =====================
    Zeros( c, 2*n+m, 1 );
    auto cuv = c( IR(0,2*n), ALL );
    Fill( cuv, lambda1 );

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
    AHat.ProcessQueues();

    // G := | -I  0 0 |
    //      |  0 -I 0 |
    // ================
    Zeros( G, 2*n, 2*n+m );
    G.Reserve( 2*m );
    for( Int e=0; e<2*m; ++e )
        G.QueueUpdate( e, e, Real(-1) );
    G.ProcessQueues();

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine QP
    // ===================
    Matrix<Real> xHat, y, z, s;
    QP( Q, AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
    x = xHat( uInd, ALL );
    x -= xHat( vInd, ALL );
}

template<typename Real>
void EN
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
        Real lambda1,
        Real lambda2,
        DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();

    DistSparseMatrix<Real> Q(grid), AHat(grid), G(grid);
    DistMultiVec<Real> c(grid), h(grid);

    // Q := | 2*lambda_2     0      0 |
    //      |     0      2*lambda_2 0 |
    //      |     0          0      2 |
    // ================================
    Zeros( Q, 2*n+m, 2*n+m );
    Q.Reserve( Q.LocalHeight() );
    for( Int iLoc=0; iLoc<Q.LocalHeight(); ++iLoc )
    {
        const Int i = Q.GlobalRow(iLoc);
        if( i < 2*n )
            Q.QueueLocalUpdate( iLoc, i, 2*lambda2 );
        else
            Q.QueueLocalUpdate( iLoc, i, Real(2) );
    }
    Q.ProcessLocalQueues();

    // c := lambda_1*[1;1;0]
    // =====================
    Zeros( c, 2*n+m, 1 );
    for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
        if( c.GlobalRow(iLoc) < 2*n )
            c.SetLocal( iLoc, 0, lambda1 );

    // \hat A := [A, -A, I]
    // ====================
    // NOTE: Since A and \hat A are the same height and each distributed within
    //       columns, it is possible to form \hat A from A without communication
    const Int numLocalEntriesA = A.NumLocalEntries();
    Zeros( AHat, m, 2*n+m );
    AHat.Reserve( 2*numLocalEntriesA+AHat.LocalHeight() );
    for( Int e=0; e<numLocalEntriesA; ++e )
    {
        AHat.QueueUpdate( A.Row(e), A.Col(e),    A.Value(e) );
        AHat.QueueUpdate( A.Row(e), A.Col(e)+n, -A.Value(e) );
    }
    for( Int iLoc=0; iLoc<AHat.LocalHeight(); ++iLoc )
    {
        const Int i = AHat.GlobalRow(iLoc);
        AHat.QueueLocalUpdate( iLoc, i+2*n, Real(1) );
    }
    AHat.ProcessLocalQueues();

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
    G.ProcessLocalQueues();

    // h := 0
    // ======
    Zeros( h, 2*n, 1 );

    // Solve the affine QP
    // ===================
    DistMultiVec<Real> xHat(grid), y(grid), z(grid), s(grid);
    QP( Q, AHat, G, b, c, h, xHat, y, z, s, ctrl );

    // x := u - v
    // ==========
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
            x.QueueUpdate( i, 0, xHat.GetLocal(iLoc,0) );
        else if( i < 2*n )
            x.QueueUpdate( i-n, 0, -xHat.GetLocal(iLoc,0) );
        else
            break;
    }
    x.ProcessQueues();
}

#define PROTO(Real) \
  template void EN \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
          Real lambda1, \
          Real lambda2, \
          Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void EN \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
          Real lambda1, \
          Real lambda2, \
          AbstractDistMatrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void EN \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
          Real lambda1, \
          Real lambda2, \
          Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void EN \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
          Real lambda1, \
          Real lambda2, \
          DistMultiVec<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
