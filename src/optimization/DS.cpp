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
// We begin with support for (DS2), due to its avoidance of explicitly forming
// A^T A.
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
    Matrix<Real> c, xHat, AHat, bHat, G, h;

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
    Matrix<Real> y, z, s;
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
    DistMatrix<Real> c(g), xHat(g), AHat(g), bHat(g), G(g), h(g);

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
    DistMatrix<Real> y(g), z(g), s(g);
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
    LogicError("This routine is not yet written");
}

template<typename Real>
void Var2
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        Real lambda,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("ds::Var2"))
    LogicError("This routine is not yet written");
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
    ds::Var2( A, b, lambda, x, ctrl );
}

template<typename Real>
void DS
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("DS"))
    ds::Var2( A, b, lambda, x, ctrl );
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
