/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GLM_HPP
#define ELEM_GLM_HPP

#include ELEM_GQR_INC

// This driver solves a sequence of General (Gauss-Markov) Linear Model (GLM)
// problems using a Generalized QR factorization. 
//
// The problem formulation is 
//     min || y ||_2 subject to d = A x + B y,
//     x,y
// where A is m x n, B is m x p, and n <= m <= n+p. It is assumed that 
// A has full column rank, n, and [A B] has full row rank, m.
//
// A Generalized QR factorization of (A,B),
//     A = Q R = Q | R11 |, B = Q T Z = Q | T11 T12 | Z,
//                 | 0   |                |   0 T22 |
// where Q and Z are unitary and R and T are upper-trapezoidal, allows us to
// re-express the constraint as 
//     (Q^H d) = | R11 | x + | T11 T12 | (Z y).
//               |   0 |     |   0 T22 |
// which is re-written as
//      | g1 | = | R11 x + T11 c1 + T12 c2 |
//      | g2 |   |                  T22 c2 |.
// Since || c ||_2 == || Z y ||_2 = || y ||_2 is to be minimized, and c2 is 
// fixed, our only freedom is in the choice of c1, which we set to zero.
// Then all that is left is to solve
//      R11 x = g1 - T12 c2
// for x.
//
// On exit, A and B are overwritten with their implicit Generalized QR 
// factorization, D is overwritten with X, and Y is set to the solution Y.
//
// Note that essentially the same scheme is used in LAPACK's {S,D,C,Z}GGGLM.

namespace elem {

template<typename F> 
inline void
GLM( Matrix<F>& A, Matrix<F>& B, Matrix<F>& D, Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("GLM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Width();
    const Int numRhs = D.Width();
    if( m != B.Height() || m != D.Height() )
        LogicError("A, B, and D must be the same height");
    if( m < n )
        LogicError("GLM requires height(A) >= width(A)");
    if( n+p < m )
        LogicError("GLM requires width(A)+width(B) >= height(A)");
    const bool checkIfSingular = true;

    // Compute the implicit Generalized QR decomposition of (A,B)
    Matrix<F> tA, tB;
    Matrix<Base<F>> dA, dB;
    GQR( A, tA, dA, B, tB, dB );

    // G := Q^H D
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, D );

    // Partition the relevant matrices
    Matrix<F> G1, G2;
    PartitionDown( D, G1, G2, n );
    Matrix<F> R11, R21;
    PartitionDown( A, R11, R21, n );
    Matrix<F> T11, T12, T21, T22;
    PartitionUpOffsetDiagonal
    ( p-m,
      B, T11, T12,
         T21, T22, m-n );
    Zeros( Y, p, numRhs );
    Matrix<F> C1, C2;
    PartitionDown( Y, C1, C2, n+p-m );

    // Solve T22 C2 = G2
    C2 = G2;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T22, C2, checkIfSingular );

    // G1 := G1 - T12 C2
    Gemm( NORMAL, NORMAL, F(-1), T12, C2, F(1), G1 );
    
    // Solve R11 X = G1 
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, G1, checkIfSingular );
    D.Resize( n, numRhs );

    // Y := Z^H C
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, Y );
}

template<typename F> 
inline void
GLM( DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<F>& D, DistMatrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("GLM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Width();
    const Int numRhs = D.Width();
    if( m != B.Height() || m != D.Height() )
        LogicError("A, B, and D must be the same height");
    if( m < n )
        LogicError("GLM requires height(A) >= width(A)");
    if( n+p < m )
        LogicError("GLM requires width(A)+width(B) >= height(A)");
    const Grid& g = A.Grid();
    if( g != B.Grid() || g != D.Grid() )
        LogicError("All matrices must have the same grid");
    Y.SetGrid( g );
    const bool checkIfSingular = true;

    // Compute the implicit Generalized QR decomposition of (A,B)
    DistMatrix<F,MD,STAR> tA(g), tB(g);
    DistMatrix<Base<F>,MD,STAR> dA(g), dB(g);
    GQR( A, tA, dA, B, tB, dB );

    // G := Q^H D
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, D );

    // Partition the relevant matrices
    DistMatrix<F> G1(g), G2(g);
    PartitionDown( D, G1, G2, n );
    DistMatrix<F> R11(g), R21(g);
    PartitionDown( A, R11, R21, n );
    DistMatrix<F> T11(g), T12(g), T21(g), T22(g);
    PartitionUpOffsetDiagonal
    ( p-m,
      B, T11, T12,
         T21, T22, m-n );
    Zeros( Y, p, numRhs );
    DistMatrix<F> C1(g), C2(g);
    PartitionDown( Y, C1, C2, n+p-m );

    // Solve T22 C2 = G2
    C2 = G2;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T22, C2, checkIfSingular );

    // G1 := G1 - T12 C2
    Gemm( NORMAL, NORMAL, F(-1), T12, C2, F(1), G1 );
    
    // Solve R11 X = G1
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, G1, checkIfSingular );
    D.Resize( n, numRhs );

    // Y := Z^H C
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, Y );
}

} // namespace elem

#endif // ifndef ELEM_GLM_HPP
