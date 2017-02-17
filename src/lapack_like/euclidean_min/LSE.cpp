/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "El/core/FlamePart.hpp"

// This file implements both dense and sparse-direct solutions of
// Equality-constrained Least Squares (LSE):
//
//     min_x || A x - c ||_2 subject to B x = d.
//
// For dense instances of the problem, a Generalized RQ factorization can be
// employed as long as A is m x n, B is p x n, and p <= n <= m+p. It is
// assumed that B has full row rank, p, and [A;B] has full column rank, n.
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
// On exit of the internal dense implementation, A and B are overwritten with
// their implicit Generalized RQ factorization of (B,A), and, optionally, C is
// overwritten with the rotated residual matrix
//     Z^H (A X - C) = (R Q X - Z^H C) = |           0 |,
//                                       | R22 Y2 - G1 |
// where R22 is an upper-trapezoidal (not necessarily triangular) matrix.
// D is overwritten with arbitrary values.
//
// Note that essentially the same scheme is used in LAPACK's {S,D,C,Z}GGLSE.
//
// For sparse instances of the LSE problem, the symmetric quasi-semidefinite
// augmented system
//
//     | 0    A^H    B^H | |     x    |   | 0 |
//     | A -alpha*I   0  | | -r/alpha | = | c |
//     | B     0      0  | |  y/alpha |   | d |
//
// is formed, equilibrated, and then a priori regularization is added in order
// to make the system sufficiently quasi-definite. A Cholesky-like factorization
// of this regularized system is then used as a preconditioner for FGMRES(k).
//

namespace El {

namespace lse {

template<typename F>
void Overwrite
( Matrix<F>& A,
  Matrix<F>& B,
  Matrix<F>& C,
  Matrix<F>& D,
  Matrix<F>& X, bool computeResidual )
{
    EL_DEBUG_CSE
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
    auto ind1 = IR(0,n-p);
    auto ind2 = IR(n-p,END);
    auto Y1 = X( ind1, ALL );
    auto Y2 = X( ind2, ALL );
    auto T11 = B( ALL, ind1 );
    auto T12 = B( ALL, ind2 );
    auto R11 = A( ind1, ind1 );
    auto R12 = A( ind1, ind2 );
    auto R21 = A( ind2, ind1 );
    auto R22 = A( ind2, ind2 );
    auto G1 = C( ind1, ALL );
    auto G2 = C( ind2, ALL );

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
            G2 -= DT;
        }
        else
        {
            Matrix<F> R22T, R22B;
            PartitionUp( R22, R22T, R22B, m-n );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22T, D );
            Matrix<F> G2T, G2B;
            PartitionUp( G2, G2T, G2B, m-n );
            G2T -= D;
        }
        Zero( G1 );
    }

    // X := Q^H Y
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, X );
}

template<typename F>
void Overwrite
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& BPre,
  AbstractDistMatrix<F>& CPre,
  AbstractDistMatrix<F>& DPre,
  AbstractDistMatrix<F>& XPre,
  bool computeResidual )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR>
      AProx( APre ),
      BProx( BPre ),
      CProx( CPre ),
      DProx( DPre );
    DistMatrixWriteProxy<F,F,MC,MR>
      XProx( XPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();
    auto& C = CProx.Get();
    auto& D = DProx.Get();
    auto& X = XProx.Get();

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
            G2 -= DT;
        }
        else
        {
            DistMatrix<F> R22T(g), R22B(g);
            PartitionUp( R22, R22T, R22B, m-n );
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R22T, D );
            DistMatrix<F> G2T(g), G2B(g);
            PartitionUp( G2, G2T, G2B, m-n );
            G2T -= D;
        }
        Zero( G1 );
    }

    // X := Q^H Y
    rq::ApplyQ( LEFT, ADJOINT, B, tB, dB, X );
}

} // namespace lse

template<typename F>
void LSE
( const Matrix<F>& A,
  const Matrix<F>& B,
  const Matrix<F>& C,
  const Matrix<F>& D,
        Matrix<F>& X )
{
    EL_DEBUG_CSE
    Matrix<F> ACopy( A ), BCopy( B ), CCopy( C ), DCopy( D );
    lse::Overwrite( ACopy, BCopy, CCopy, DCopy, X );
}

template<typename F>
void LSE
( const AbstractDistMatrix<F>& A,
  const AbstractDistMatrix<F>& B,
  const AbstractDistMatrix<F>& C,
  const AbstractDistMatrix<F>& D,
        AbstractDistMatrix<F>& X )
{
    EL_DEBUG_CSE
    DistMatrix<F> ACopy( A ), BCopy( B ), CCopy( C ), DCopy( D );
    lse::Overwrite( ACopy, BCopy, CCopy, DCopy, X );
}

template<typename F>
void LSE
( const SparseMatrix<F>& A,
  const SparseMatrix<F>& B,
  const Matrix<F>& C,
  const Matrix<F>& D,
        Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Height();
    const Int numRHS = C.Width();

    // Rescale W = [ A; B ]
    // ====================
    Matrix<Real> dC;
    SparseMatrix<F> W;
    VCat( A, B, W );
    if( ctrl.equilibrate )
    {
        ColumnTwoNorms( W, dC );
        auto normMap = []( const Real& beta )
          { return beta < Sqrt(limits::Epsilon<Real>()) ? Real(1) : beta; };
        EntrywiseMap( dC, MakeFunction(normMap) );
        DiagonalSolve( RIGHT, NORMAL, dC, W );
    }
    else
    {
        Ones( dC, n, 1 );
    }
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        normScale = TwoNormEstimate( W, ctrl.basisSize );
        if( ctrl.progress )
            cout << "Estimated || [ A; B ] ||_2 ~= " << normScale << endl;
        W *= F(1)/normScale;
        dC *= normScale;
    }

    // Form the augmented RHS
    // ======================
    //   G = [ 0; C; D ]
    Matrix<F> G;
    Zeros( G, n+m+k, numRHS );
    {
        auto Gr = G( IR(n,n+m),     ALL );
        auto Gy = G( IR(n+m,n+m+k), ALL );
        Gr = C;
        Gy = D;
    }

    // Form the augmented matrix
    // =========================
    //
    //         | 0     A^H    B^H |
    //     J = | A  -alpha*I   0  |
    //         | B      0      0  |
    //
    const Int numEntriesW = W.NumEntries();
    SparseMatrix<F> J;
    Zeros( J, n+m+k, n+m+k );
    J.Reserve( 2*numEntriesW+m );
    for( Int e=0; e<numEntriesW; ++e )
    {
        J.QueueUpdate( W.Row(e)+n, W.Col(e),        W.Value(e)  );
        J.QueueUpdate( W.Col(e),   W.Row(e)+n, Conj(W.Value(e)) );
    }
    for( Int e=0; e<m; ++e )
        J.QueueUpdate( e+n, e+n, -ctrl.alpha );
    J.ProcessQueues();

    // Solve the Symmetric Quasi-SemiDefinite system
    // =============================================
    SQSDSolve( n, J, G, ctrl.sqsdCtrl );

    // Extract X from G = [ Dc*X; -R/alpha; Y/alpha ]
    // ==============================================
    X = G( IR(0,n), ALL );
    DiagonalSolve( LEFT, NORMAL, dC, X );
}

template<typename F>
void LSE
( const DistSparseMatrix<F>& A,
  const DistSparseMatrix<F>& B,
  const DistMultiVec<F>& C,
  const DistMultiVec<F>& D,
        DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Height();
    const Int numRHS = C.Width();
    const Grid& grid = A.Grid();
    const int commRank = grid.Rank();

    // Rescale W = [ A; B ]
    // ====================
    DistMultiVec<Real> dC(grid);
    DistSparseMatrix<F> W(grid);
    VCat( A, B, W );
    if( ctrl.equilibrate )
    {
        ColumnTwoNorms( W, dC );
        auto normMap = []( const Real& beta )
          { return beta < Sqrt(limits::Epsilon<Real>()) ? Real(1) : beta; };
        EntrywiseMap( dC, MakeFunction(normMap) );
        DiagonalSolve( RIGHT, NORMAL, dC, W );
    }
    else
    {
        Ones( dC, n, 1 );
    }
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        normScale = TwoNormEstimate( W, ctrl.basisSize );
        if( ctrl.progress && commRank == 0 )
            cout << "Estimated || [ A; B ] ||_2 ~= " << normScale << endl;
        W *= F(1)/normScale;
        dC *= normScale;
    }

    // Form the augmented RHS
    // ======================
    //   G = [ 0; C; D ]
    DistMultiVec<F> G(grid);
    Zeros( G, n+m+k, numRHS );
    {
        const Int CLocalHeight = C.LocalHeight();
        const Int DLocalHeight = D.LocalHeight();
        const Int sendCount = (CLocalHeight+DLocalHeight)*numRHS;
        G.Reserve( sendCount );
        for( Int iLoc=0; iLoc<CLocalHeight; ++iLoc )
        {
            const Int i = C.GlobalRow(iLoc);
            for( Int j=0; j<numRHS; ++j )
                G.QueueUpdate( i+n, j, C.GetLocal(iLoc,j) );
        }
        for( Int iLoc=0; iLoc<DLocalHeight; ++iLoc )
        {
            const Int i = D.GlobalRow(iLoc);
            for( Int j=0; j<numRHS; ++j )
                G.QueueUpdate( i+n+m, j, D.GetLocal(iLoc,j) );
        }
        G.ProcessQueues();
    }

    // Form the augmented matrix
    // =========================
    //
    //         | 0     A^H    B^H |
    //     J = | A  -alpha*I   0  |
    //         | B      0      0  |
    //
    const Int numEntriesW = W.NumLocalEntries();
    DistSparseMatrix<F> J(grid);
    Zeros( J, n+m+k, n+m+k );
    {
        const Int JLocalHeight = J.LocalHeight();
        Int numNegAlphaUpdates = 0;
        for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i >= n+m )
                break;
            else if( i >= n )
                ++numNegAlphaUpdates;
        }
        const Int numSend = 2*numEntriesW;
        J.Reserve( numSend+numNegAlphaUpdates, numSend );

        for( Int e=0; e<numEntriesW; ++e )
        {
            const Int i = W.Row(e);
            const Int j = W.Col(e);
            const F value = W.Value(e);

            J.QueueUpdate( i+n, j, value );
            J.QueueUpdate( j, i+n, Conj(value) );
        }
        for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i >= n+m )
                break;
            else if( i >= n )
                J.QueueLocalUpdate( iLoc, i, -ctrl.alpha );
        }

        J.ProcessQueues();
    }

    // Solve the Symmetric Quasi-SemiDefinite system
    // =============================================
    SQSDSolve( n, J, G, ctrl.sqsdCtrl );

    // Extract X from G = [ Dc*X; -R/alpha; Y/alpha ]
    // ==============================================
    X = G( IR(0,n), ALL );
    DiagonalSolve( LEFT, NORMAL, dC, X );
}

#define PROTO(F) \
  template void lse::Overwrite \
  ( Matrix<F>& A, \
    Matrix<F>& B, \
    Matrix<F>& C, \
    Matrix<F>& D, \
    Matrix<F>& X, bool computeResidual ); \
  template void lse::Overwrite \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<F>& C, \
    AbstractDistMatrix<F>& D, \
    AbstractDistMatrix<F>& X, bool computeResidual ); \
  template void LSE \
  ( const Matrix<F>& A, \
    const Matrix<F>& B, \
    const Matrix<F>& C, \
    const Matrix<F>& D, \
          Matrix<F>& X ); \
  template void LSE \
  ( const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& B, \
    const AbstractDistMatrix<F>& C, \
    const AbstractDistMatrix<F>& D, \
          AbstractDistMatrix<F>& X ); \
  template void LSE \
  ( const SparseMatrix<F>& A, \
    const SparseMatrix<F>& B, \
    const Matrix<F>& C, \
    const Matrix<F>& D, \
          Matrix<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl ); \
  template void LSE \
  ( const DistSparseMatrix<F>& A, \
    const DistSparseMatrix<F>& B, \
    const DistMultiVec<F>& C, \
    const DistMultiVec<F>& D, \
          DistMultiVec<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
