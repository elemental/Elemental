/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace ls {

template<typename F>
void Overwrite
( Orientation orientation,
        Matrix<F>& A,
  const Matrix<F>& B,
        Matrix<F>& X )
{
    EL_DEBUG_CSE

    Matrix<F> phase;
    Matrix<Base<F>> signature;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, phase, signature );
        qr::SolveAfter( orientation, A, phase, signature, B, X );
    }
    else
    {
        LQ( A, phase, signature );
        lq::SolveAfter( orientation, A, phase, signature, B, X );
    }
}

template<typename F>
void Overwrite
( Orientation orientation,
        AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& B,
        AbstractDistMatrix<F>& X )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<F,MD,STAR> phase(A.Grid());
    DistMatrix<Base<F>,MD,STAR> signature(A.Grid());

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, phase, signature );
        qr::SolveAfter( orientation, A, phase, signature, B, X );
    }
    else
    {
        LQ( A, phase, signature );
        lq::SolveAfter( orientation, A, phase, signature, B, X );
    }
}

} // namespace ls

template<typename F>
void LeastSquares
( Orientation orientation,
  const Matrix<F>& A,
  const Matrix<F>& B,
        Matrix<F>& X )
{
    EL_DEBUG_CSE
    Matrix<F> ACopy( A );
    ls::Overwrite( orientation, ACopy, B, X );
}

template<typename F>
void LeastSquares
( Orientation orientation,
  const AbstractDistMatrix<F>& A,
  const AbstractDistMatrix<F>& B,
        AbstractDistMatrix<F>& X )
{
    EL_DEBUG_CSE
    DistMatrix<F> ACopy( A );
    ls::Overwrite( orientation, ACopy, B, X );
}

// The following routines solve either
//
//   Minimum length:
//     min_X || X ||_F
//     s.t. W X = B, or
//
//   Least squares:
//     min_X || W X - B ||_F,
//
// where W=op(A) is either A, A^T, or A^H, via forming a Hermitian
// quasi-semidefinite system
//
//    | alpha*I  W | | R/alpha | = | B |,
//    |   W^H    0 | | X       |   | 0 |
//
// when height(W) >= width(W), or
//
//    | alpha*I  W^H | |     X   | = | 0 |,
//    |   W       0  | | alpha*Y |   | B |
//
// when height(W) < width(W).
//
// The latter guarantees that W X = B and X in range(W^H), which shows that
// X solves the minimum length problem. The former defines R = B - W X and
// ensures that R is in the null-space of W^H (therefore solving the least
// squares problem).
//
// Note that, ideally, alpha is roughly the minimum (nonzero) singular value
// of W, which implies that the condition number of the quasi-semidefinite
// system is roughly equal to the condition number of W (see the analysis of
// Bjorck). If it is too expensive to estimate the minimum singular value, and
// W is equilibrated to have a unit two-norm, a typical choice for alpha is
// epsilon^0.25.
//
// The Hermitian quasi-semidefinite systems are solved by converting them into
// Hermitian quasi-definite form via a priori regularization, applying an
// LDL^H factorization with static pivoting to the regularized system, and
// using the iteratively-refined solution of with the regularized factorization
// as a preconditioner for the original problem (defaulting to Flexible GMRES
// for now).
//
// This approach originated within
//
//    Michael Saunders,
//   "Chapter 8, Cholesky-based Methods for Sparse Least Squares:
//    The Benefits of Regularization",
//    in L. Adams and J.L. Nazareth (eds.), Linear and Nonlinear Conjugate
//    Gradient-Related Methods, SIAM, Philadelphia, 92--100 (1996).
//
// But note that SymmLQ and LSQR were used rather than flexible GMRES, and
// iteratively refining *within* the preconditioner was not discussed.
//

// NOTE: The following routines are implemented as a special case of Tikhonov
//       regularization with either an m x 0 or 0 x n regularization matrix.

namespace ls {

template<typename F>
void Equilibrated
( const SparseMatrix<F>& A,
  const Matrix<F>& B,
        Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numRHS = B.Width();
    const Int numEntriesA = A.NumEntries();

    SparseMatrix<F> J;
    Zeros( J, m+n, m+n );
    J.Reserve( 2*numEntriesA + Max(m,n) );
    if( m >= n )
    {
        // Form J = [alpha*I, A; A^H, 0]
        // =============================
        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( A.Row(e),   A.Col(e)+m,      A.Value(e)  );
            J.QueueUpdate( A.Col(e)+m, A.Row(e),   Conj(A.Value(e)) );
        }
        for( Int e=0; e<m; ++e )
            J.QueueUpdate( e, e, ctrl.alpha );
    }
    else
    {
        // Form J = [alpha*I, A^H; A, 0]
        // =============================
        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, Conj(A.Value(e)) );
            J.QueueUpdate( A.Row(e)+n, A.Col(e),        A.Value(e)  );
        }
        for( Int e=0; e<n; ++e )
            J.QueueUpdate( e, e, ctrl.alpha );
    }
    J.ProcessQueues();

    Matrix<F> D;
    Zeros( D, m+n, numRHS );
    if( m >= n )
    {
        // Form D = [B; 0]
        // ==================
        auto DT = D( IR(0,m), ALL );
        DT = B;
    }
    else
    {
        // Form D = [0; B]
        // ===============
        auto DB = D( IR(n,m+n), ALL );
        DB = B;
    }

    // Solve the Symmetric Quasi-SemiDefinite system
    // =============================================
    SQSDSolve( Max(m,n), J, D, ctrl.sqsdCtrl );

    Zeros( X, n, numRHS );
    if( m >= n )
    {
        // Extract X from [R/alpha; X]
        // ===========================
        auto DB = D( IR(m,m+n), ALL );
        X = DB;
    }
    else
    {
        // Extract X from [X; alpha*Y]
        // ===========================
        auto DT = D( IR(0,n), ALL );
        X = DT;
    }
}

} // namespace ls

template<typename F>
void LeastSquares
( Orientation orientation,
  const SparseMatrix<F>& A,
  const Matrix<F>& B,
        Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    SparseMatrix<F> ABar;
    if( orientation == NORMAL )
        ABar = A;
    else if( orientation == TRANSPOSE )
        Transpose( A, ABar );
    else
        Adjoint( A, ABar );
    auto BBar = B;
    const Int m = ABar.Height();
    const Int n = ABar.Width();

    // Equilibrate the matrix
    // ======================
    Matrix<Real> dR, dC;
    if( ctrl.equilibrate )
    {
        auto normMap = []( const Real& beta )
          { return beta < Sqrt(limits::Epsilon<Real>()) ? Real(1) : beta; };
        if( m >= n )
        {
            ColumnTwoNorms( ABar, dC );
            EntrywiseMap( dC, MakeFunction(normMap) );
            DiagonalSolve( RIGHT, NORMAL, dC, ABar );
            Ones( dR, m, 1 );
        }
        else
        {
            RowTwoNorms( ABar, dR );
            EntrywiseMap( dR, MakeFunction(normMap) );
            DiagonalSolve( LEFT, NORMAL, dR, ABar );
            Ones( dC, n, 1 );
        }
    }
    else
    {
        Ones( dR, m, 1 );
        Ones( dC, n, 1 );
    }
    if( ctrl.scaleTwoNorm )
    {
        // Scale ABar down to roughly unit two-norm
        const Real normScale = TwoNormEstimate( ABar, ctrl.basisSize );
        if( ctrl.progress )
            Output("Estimated || A ||_2 ~= ",normScale);
        ABar *= F(1)/normScale;
        dR *= normScale;
    }

    // Equilibrate the RHS
    // ===================
    DiagonalSolve( LEFT, NORMAL, dR, BBar );

    // Solve the equilibrated least squares problem
    // ============================================
    ls::Equilibrated( ABar, BBar, X, ctrl );

    // Unequilibrate the solution
    // ==========================
    DiagonalSolve( LEFT, NORMAL, dC, X );
}

namespace ls {

template<typename F>
void Equilibrated
( const DistSparseMatrix<F>& A,
  const DistMultiVec<F>& B,
        DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )
    const Grid& grid = A.Grid();

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numRHS = B.Width();

    // J := [alpha*I,A;A^H,0] or [alpha*I,A^H;A,0]
    // ===========================================
    DistSparseMatrix<F> J(grid);
    Zeros( J, m+n, m+n );
    {
        const Int JLocalHeight = J.LocalHeight();
        const Int numLocalEntriesA = A.NumLocalEntries();
        Int numAlphaUpdates = 0;
        for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i < Max(m,n) )
                ++numAlphaUpdates;
            else
                break;
        }
        const Int numSend = 2*numLocalEntriesA;
        J.Reserve( numSend+numAlphaUpdates, numSend );
        for( Int e=0; e<numLocalEntriesA; ++e )
        {
            const Int i = A.Row(e);
            const Int j = A.Col(e);
            const F value = A.Value(e);
            if( m >= n )
            {
                J.QueueUpdate( i, j+m, value );
                J.QueueUpdate( j+m, i, Conj(value) );
            }
            else
            {
                J.QueueUpdate( i+n, j, value );
                J.QueueUpdate( j, i+n, Conj(value) );
            }
        }
        for( Int iLoc=0; iLoc<JLocalHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i < Max(m,n) )
                J.QueueLocalUpdate( iLoc, i, ctrl.alpha );
            else
                break;
        }
        J.ProcessQueues();
    }

    // Set D to [B; 0] or [0; B]
    // =========================
    DistMultiVec<F> D(grid);
    Zeros( D, m+n, numRHS );
    {
        const Int BLocalHeight = B.LocalHeight();
        const Int numSend = BLocalHeight*numRHS;
        D.Reserve( numSend );
        for( Int iLoc=0; iLoc<BLocalHeight; ++iLoc )
        {
            const Int i = B.GlobalRow(iLoc);
            if( m >= n )
            {
                for( Int j=0; j<numRHS; ++j )
                    D.QueueUpdate( i, j, B.GetLocal(iLoc,j) );
            }
            else
            {
                for( Int j=0; j<numRHS; ++j )
                    D.QueueUpdate( i+n, j, B.GetLocal(iLoc,j) );
            }
        }
        D.ProcessQueues();
    }

    // Solve the Symmetric Quasi-SemiDefinite system
    // =============================================
    SQSDSolve( Max(m,n), J, D, ctrl.sqsdCtrl );

    // Extract X from [R/alpha; X] or [X; alpha*Y] and then rescale
    // ============================================================
    if( m >= n )
        X = D( IR(m,m+n), ALL );
    else
        X = D( IR(0,n),   ALL );
}

} // namespace ls

template<typename F>
void LeastSquares
( Orientation orientation,
  const DistSparseMatrix<F>& A,
  const DistMultiVec<F>& B,
        DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Grid& grid = A.Grid();
    const int commRank = grid.Rank();

    DistSparseMatrix<F> ABar(grid);
    if( orientation == NORMAL )
        ABar = A;
    else if( orientation == TRANSPOSE )
        Transpose( A, ABar );
    else
        Adjoint( A, ABar );
    auto BBar = B;
    const Int m = ABar.Height();
    const Int n = ABar.Width();

    // Equilibrate the matrix
    // ======================
    DistMultiVec<Real> dR(grid), dC(grid);
    if( ctrl.equilibrate )
    {
        auto normMap = []( const Real& beta )
          { return beta < Sqrt(limits::Epsilon<Real>()) ? Real(1) : beta; };
        if( m >= n )
        {
            ColumnTwoNorms( ABar, dC );
            EntrywiseMap( dC, MakeFunction(normMap) );
            DiagonalSolve( RIGHT, NORMAL, dC, ABar );
            Ones( dR, m, 1 );
        }
        else
        {
            RowTwoNorms( ABar, dR );
            EntrywiseMap( dR, MakeFunction(normMap) );
            DiagonalSolve( LEFT, NORMAL, dR, ABar );
            Ones( dC, n, 1 );
        }
    }
    else
    {
        Ones( dR, m, 1 );
        Ones( dC, n, 1 );
    }
    if( ctrl.scaleTwoNorm )
    {
        // Scale ABar down to roughly unit two-norm
        const Real normScale = TwoNormEstimate( ABar, ctrl.basisSize );
        if( ctrl.progress && commRank == 0 )
            Output("Estimated || A ||_2 ~= ",normScale);
        ABar *= F(1)/normScale;
        dR *= normScale;
    }

    // Equilibrate the RHS
    // ===================
    DiagonalSolve( LEFT, NORMAL, dR, BBar );

    // Solve the equilibrated least squares problem
    // ============================================
    ls::Equilibrated( ABar, BBar, X, ctrl );

    // Unequilibrate the solution
    // ==========================
    DiagonalSolve( LEFT, NORMAL, dC, X );
}

#define PROTO(F) \
  template void ls::Overwrite \
  ( Orientation orientation, \
          Matrix<F>& A, \
    const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void ls::Overwrite \
  ( Orientation orientation, \
          AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& B, \
          AbstractDistMatrix<F>& X ); \
  template void LeastSquares \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void LeastSquares \
  ( Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& B, \
          AbstractDistMatrix<F>& X ); \
  template void LeastSquares \
  ( Orientation orientation, \
    const SparseMatrix<F>& A, \
    const Matrix<F>& B, \
          Matrix<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl ); \
  template void LeastSquares \
  ( Orientation orientation, \
    const DistSparseMatrix<F>& A, \
    const DistMultiVec<F>& B, \
          DistMultiVec<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
