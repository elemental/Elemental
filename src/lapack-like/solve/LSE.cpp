/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

// This driver solves a sequence of Equality-constrained Least Squares (LSE)
// problems using a Generalized RQ factorization. 
//
// The problem formulation is 
//     min || A x - c ||_2 subject to B x = d,
//      x
// where A is m x n, B is p x n, and p <= n <= m+p. It is assumed that 
// B has full row rank, p, and [A;B] has full column rank, n.
//
// A Generalized RQ factorization of (B,A),
//    B = T Q = | 0 T12 | Q,  A = Z | R11 R12 | Q,
//                                  |   0 R22 |
// where Q and Z are unitary and R and T are upper-trapezoidal, allows us to
// re-express the constraint
//     T Q x = d,
// as
//     | 0 T12 | | y1 | = d,
//               | y2 |   
// where y = Q x, which only requires the solution of the upper-triangular 
// system
//     T12 y2 = d.
//
// The objective can be rewritten as
//     || A x - c ||_2 = || Z^H A x - Z^H c ||_2
//                     = ||   R Q x - Z^H c ||_2
// which, defining g = Z^H c, can be partitioned as
//     | R11 R12 | | y1 | - | g1 | = | R11 y1 + R12 y2 - g1 |.
//     |   0 R22 | | y2 |   | g2 |   |          R22 y2 - g2 |
// Since y2 is fixed by the constraint, the norm is minimized by setting the
// top term to zero, which involves solving the upper-triangular system
//     R11 y1 = g1 - R12 y2.
//       
// On exit, A and B are overwritten with their implicit Generalized RQ 
// factorization of (B,A), and, optionally, C is overwritten with the rotated
// residual matrix
//     Z^H (A X - C) = (R Q X - Z^H C) = |           0 |,
//                                       | R22 Y2 - G1 |
// where R22 is an upper-trapezoidal (not necessarily triangular) matrix.
// D is overwritten with arbitrary values.
//
// Note that essentially the same scheme is used in LAPACK's {S,D,C,Z}GGLSE.

namespace El {

template<typename F> 
void LSE
( Matrix<F>& A, Matrix<F>& B, Matrix<F>& C, Matrix<F>& D, Matrix<F>& X, 
  bool computeResidual )
{
    DEBUG_ONLY(CallStackEntry cse("LSE"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Height();
    const Int numRhs = D.Width();
    if( m != C.Height() )
        LogicError("A and C must be the same height");
    if( p != D.Height() )
        LogicError("B and D must be the same height");
    if( numRhs != C.Width() )
        LogicError("C and D must be the same width");
    if( n < p )
        LogicError("LSE requires width(A) >= height(B)");
    if( m+p < n )
        LogicError("LSE requires height(A)+height(B) >= width(A)");
    const bool checkIfSingular = true;

    // Compute the implicit Generalized RQ decomposition of (B,A)
    Matrix<F> tA, tB;
    Matrix<Base<F>> dA, dB;
    GRQ( B, tB, dB, A, tA, dA );

    // G := Z^H C
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, C );

    // Partition the relevant matrices
    Zeros( X, n, numRhs );
    Matrix<F> Y1, Y2;
    PartitionUp( X, Y1, Y2, p );
    Matrix<F> T11, T12;
    PartitionLeft( B, T11, T12, p );
    Matrix<F> R11, R12, R21, R22;
    PartitionDownDiagonal( A, R11, R12, R21, R22, n-p );
    Matrix<F> G1, G2;
    PartitionDown( C, G1, G2, n-p );

    // Solve T12 Y2 = D
    Y2 = D; 
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T12, Y2, checkIfSingular );

    // G1 := G1 - R12 Y2
    Gemm( NORMAL, NORMAL, F(-1), R12, Y2, F(1), G1 );

    // Solve R11 Y1 = G1
    Y1 = G1;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, Y1, checkIfSingular );

    if( computeResidual )
    {
        // R22 is upper-trapezoidal, and so it is best to decompose it in terms
        // of its upper-left triangular block and either its bottom zero 
        // block or right non-zero block. Putting k=Min(p,m-(n-p)), then
        // the k x k upper-left block is upper-triangular. If m >= n, the
        // bottom m-(n-p) - k = m-n rows are zero, otherwise the right 
        // p - k = n-m.columns are nonzero.
        if( m < n )
        {
            Matrix<F> R22L, R22R;
            PartitionLeft( R22, R22L, R22R, n-m );
            Matrix<F> DT, DB;
            PartitionUp( D, DT, DB, n-m );
            Gemm( NORMAL, NORMAL, F(-1), R22R, DB, F(1), G2 );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22L, DT );
            Axpy( F(-1), DT, G2 );
        }
        else
        {
            Matrix<F> R22T, R22B;
            PartitionUp( R22, R22T, R22B, m-n );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22T, D );
            Matrix<F> G2T, G2B;
            PartitionUp( G2, G2T, G2B, m-n );
            Axpy( F(-1), D, G2T );
        }
        Zero( G1 );
    }

    // X := Q^H Y
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, X );
}

template<typename F> 
void LSE
( DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<F>& C, DistMatrix<F>& D, 
  DistMatrix<F>& X, bool computeResidual )
{
    DEBUG_ONLY(CallStackEntry cse("LSE"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int p = B.Height();
    const Int numRhs = D.Width();
    if( m != C.Height() )
        LogicError("A and C must be the same height");
    if( p != D.Height() )
        LogicError("B and D must be the same height");
    if( numRhs != C.Width() )
        LogicError("C and D must be the same width");
    if( n < p )
        LogicError("LSE requires width(A) >= height(B)");
    if( m+p < n )
        LogicError("LSE requires height(A)+height(B) >= width(A)");
    const Grid& g = A.Grid();
    if( g != B.Grid() || g != C.Grid() || g != D.Grid() )
        LogicError("All matrices must be distributed over the same grid");
    X.SetGrid( g );
    const bool checkIfSingular = true;

    // Compute the implicit Generalized RQ decomposition of (B,A)
    DistMatrix<F,MD,STAR> tA(g), tB(g);
    DistMatrix<Base<F>,MD,STAR> dA(g), dB(g);
    GRQ( B, tB, dB, A, tA, dA );

    // G := Z^H C
    qr::ApplyQ( LEFT, ADJOINT, A, tA, dA, C );

    // Partition the relevant matrices
    Zeros( X, n, numRhs );
    DistMatrix<F> Y1(g), Y2(g);
    PartitionUp( X, Y1, Y2, p );
    DistMatrix<F> T11(g), T12(g);
    PartitionLeft( B, T11, T12, p );
    DistMatrix<F> R11(g), R12(g), R21(g), R22(g);
    PartitionDownDiagonal( A, R11, R12, R21, R22, n-p );
    DistMatrix<F> G1(g), G2(g);
    PartitionDown( C, G1, G2, n-p );

    // Solve T12 Y2 = D
    Y2 = D; 
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), T12, Y2, checkIfSingular );

    // G1 := G1 - R12 Y2
    Gemm( NORMAL, NORMAL, F(-1), R12, Y2, F(1), G1 );

    // Solve R11 Y1 = G1
    Y1 = G1;
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R11, Y1, checkIfSingular );

    if( computeResidual )
    {
        // R22 is upper-trapezoidal, and so it is best to decompose it in terms
        // of its upper-left triangular block and either its bottom zero 
        // block or right non-zero block. Putting k=Min(p,m-(n-p)), then
        // the k x k upper-left block is upper-triangular. If m >= n, the
        // bottom m-(n-p) - k = m-n rows are zero, otherwise the right 
        // p - k = n-m.columns are nonzero.
        if( m < n )
        {
            DistMatrix<F> R22L(g), R22R(g);
            PartitionLeft( R22, R22L, R22R, n-m );
            DistMatrix<F> DT(g), DB(g);
            PartitionUp( D, DT, DB, n-m );
            Gemm( NORMAL, NORMAL, F(-1), R22R, DB, F(1), G2 );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22L, DT );
            Axpy( F(-1), DT, G2 );
        }
        else
        {
            DistMatrix<F> R22T(g), R22B(g);
            PartitionUp( R22, R22T, R22B, m-n );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22T, D );
            DistMatrix<F> G2T(g), G2B(g);
            PartitionUp( G2, G2T, G2B, m-n );
            Axpy( F(-1), D, G2T );
        }
        Zero( G1 );
    }

    // X := Q^H Y
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, X );
}

#define PROTO(F) \
  template void LSE \
  ( Matrix<F>& A, Matrix<F>& B, Matrix<F>& C, Matrix<F>& D, \
    Matrix<F>& X, bool computeResidual ); \
  template void LSE \
  ( DistMatrix<F>& A, DistMatrix<F>& B, DistMatrix<F>& C, DistMatrix<F>& D, \
    DistMatrix<F>& X, bool computeResidual ); \

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
