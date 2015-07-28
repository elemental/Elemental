/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace ls {

template<typename F> 
void Overwrite
( Orientation orientation, 
        Matrix<F>& A, 
  const Matrix<F>& B, 
        Matrix<F>& X )
{
    DEBUG_ONLY(CSE cse("ls::Overwrite"))

    Matrix<F> t;
    Matrix<Base<F>> d;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, t, d );
        qr::SolveAfter( orientation, A, t, d, B, X );
    }
    else
    {
        LQ( A, t, d );
        lq::SolveAfter( orientation, A, t, d, B, X );
    }
}

template<typename F> 
void Overwrite
( Orientation orientation, 
        AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& B, 
        AbstractDistMatrix<F>& X )
{
    DEBUG_ONLY(CSE cse("ls::Overwrite"))

    auto APtr = ReadProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());

    const Int m = A.Height();
    const Int n = A.Width();
    if( m >= n )
    {
        QR( A, t, d );
        qr::SolveAfter( orientation, A, t, d, B, X );
    }
    else
    {
        LQ( A, t, d );
        lq::SolveAfter( orientation, A, t, d, B, X );
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
    DEBUG_ONLY(CSE cse("LeastSquares"))
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
    DEBUG_ONLY(CSE cse("LeastSquares"))
    DistMatrix<F> ACopy( A );
    ls::Overwrite( orientation, ACopy, B, X ); 
}

// The following routines approximately solve either
//
//   Minimum length: 
//
//     min_X { || X ||_F : W X = B },
//
// or
//
//   Least squares:  
//
//     min_X || W X - B ||_F,
//
// which are both special cases of 
//
//   min_{X,R} { || gamma X ||_F^2 + || R ||_F^2 : W X + delta R = B },
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

// TODO: Also add support for normal equations
template<typename F>
inline void Equilibrated
( const SparseMatrix<F>& A,
  const Matrix<F>& B, 
        Matrix<F>& X,
  Base<F> alpha, 
  Base<F> damp,
  Base<F> dampTmp,
  const RegLDLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("ls::Equilibrated");
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )
    typedef Base<F> Real;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numRHS = B.Width();
    const Int numEntriesA = A.NumEntries();

    Matrix<Real> reg;
    SparseMatrix<F> J;
    Zeros( reg, m+n, 1 );
    Zeros( J, m+n, m+n );
    J.Reserve( 2*numEntriesA + m+n );
    if( m >= n )
    {
        // Form J = [alpha*I, A; A^H, -damp^2*I]
        // =====================================
        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( A.Row(e),   A.Col(e)+m,      A.Value(e)  );
            J.QueueUpdate( A.Col(e)+m, A.Row(e),   Conj(A.Value(e)) );
        }
        for( Int i=0; i<m+n; ++i )
        {
            if( i < m )         
            {
                J.QueueUpdate( i, i, alpha );
            }
            else
            {
                J.QueueUpdate( i, i, -damp*damp );
                reg.Set( i, 0, -dampTmp*dampTmp );
            }
        }
    }
    else
    {
        // Form J = [alpha*I, A^H; A, -damp^2*I]
        // =====================================
        for( Int e=0; e<numEntriesA; ++e )
        {
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, Conj(A.Value(e)) );
            J.QueueUpdate( A.Row(e)+n, A.Col(e),        A.Value(e)  );
        }
        for( Int i=0; i<m+n; ++i )
        {
            if( i < n )
            {
                J.QueueUpdate( i, i, alpha );
            }
            else
            {
                J.QueueUpdate( i, i, -damp*damp );
                reg.Set( i, 0, -dampTmp*dampTmp );
            }
        }
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

    // Compute the regularized quasi-semidefinite fact of J
    // ====================================================
    SparseMatrix<F> JOrig;
    JOrig = J;
    UpdateRealPartOfDiagonal( J, Real(1), reg, 0, true );

    vector<Int> map, invMap;
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    ldl::NestedDissection( J.LockedGraph(), map, rootSep, info );
    InvertMap( map, invMap );
    ldl::Front<F> JFront( J, map, info );
    LDL( info, JFront );

    // Successively solve each of the linear systems
    // =============================================
    // TODO: Extend the iterative refinement to handle multiple RHS
    Matrix<F> u;
    Zeros( u, m+n, 1 );
    for( Int j=0; j<numRHS; ++j )
    {
        auto d = D( ALL, IR(j) );
        u = d;
        ldl::SolveAfter( JOrig, reg, invMap, info, JFront, u, ctrl );
        d = u;
    }

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
    DEBUG_ONLY(CSE cse("LeastSquares"))
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
        auto normMap = []( Real beta ) 
          { return ( beta < Sqrt(Epsilon<Real>()) ? Real(1) : beta ); };
        if( m >= n )
        {
            ColumnTwoNorms( ABar, dC );
            EntrywiseMap( dC, function<Real(Real)>(normMap) );
            DiagonalSolve( RIGHT, NORMAL, dC, ABar );
            Ones( dR, m, 1 );
        }
        else
        {
            RowTwoNorms( ABar, dR );
            EntrywiseMap( dR, function<Real(Real)>(normMap) );
            DiagonalSolve( LEFT, NORMAL, dR, ABar );
            Ones( dC, n, 1 );
        }
    }
    else
    {
        Ones( dR, m, 1 );
        Ones( dC, n, 1 );
    }
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        // Scale ABar down to roughly unit two-norm
        normScale = TwoNormEstimate( ABar, ctrl.basisSize ); 
        if( ctrl.progress )
            cout << "Estimated || A ||_2 ~= " << normScale << endl;
        ABar *= F(1)/normScale;
        dR *= normScale;
    }

    // Equilibrate the RHS
    // ===================
    DiagonalSolve( LEFT, NORMAL, dR, BBar );

    // Solve the equilibrated least squares problem
    // ============================================
    ls::Equilibrated
    ( ABar, BBar, X, ctrl.alpha, ctrl.damp, ctrl.dampTmp, ctrl.regLDLCtrl );

    // Unequilibrate the solution
    // ==========================
    DiagonalSolve( LEFT, NORMAL, dC, X );
}

namespace ls {

// TODO: Also add support for normal equations
template<typename F>
void Equilibrated
( const DistSparseMatrix<F>& A,
  const DistMultiVec<F>& B, 
        DistMultiVec<F>& X,
  Base<F> alpha,
  Base<F> damp,
  Base<F> dampTmp,
  const RegLDLCtrl<Base<F>>& ctrl,
  bool time )
{
    DEBUG_ONLY(
      CSE cse("ls::Equilibrated");
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )
    typedef Base<F> Real;
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);
    Timer timer;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int numRHS = B.Width();

    // J := [alpha*I,A;A^H,-damp^2] or [alpha*I,A^H;A,-damp^2]
    // =======================================================
    DistMultiVec<Real> reg(comm);
    DistSparseMatrix<F> J(comm);
    Zeros( reg, m+n, 1 );
    Zeros( J, m+n, m+n );
    const Int localHeight = J.LocalHeight();
    const Int numLocalEntriesA = A.NumLocalEntries();
    J.Reserve( localHeight+2*numLocalEntriesA, 2*numLocalEntriesA );
    // Queue the diagonal
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = J.GlobalRow(iLoc); 
        if( i < Max(m,n) )
            J.QueueLocalUpdate( iLoc, i, alpha );
        else
        {
            J.QueueLocalUpdate( iLoc, i, -damp*damp );
            reg.SetLocal( iLoc, 0, -dampTmp*dampTmp );
        }
    }
    // Queue A and A^T
    for( Int e=0; e<numLocalEntriesA; ++e )
    {
        if( m >= n )
        {
            J.QueueUpdate( A.Row(e),   A.Col(e)+m,      A.Value(e),  false );
            J.QueueUpdate( A.Col(e)+m, A.Row(e),   Conj(A.Value(e)), false );
        }
        else
        {
            J.QueueUpdate( A.Row(e)+n, A.Col(e),        A.Value(e),  false  );
            J.QueueUpdate( A.Col(e),   A.Row(e)+n, Conj(A.Value(e)), false );
        }
    }
    J.ProcessQueues();

    // Set D to [B; 0] or [0; B]
    // =========================
    DistMultiVec<F> D(comm);
    Zeros( D, m+n, numRHS );
    const Int BLocalHeight = B.LocalHeight();
    D.Reserve( BLocalHeight*numRHS );
    for( Int iLoc=0; iLoc<BLocalHeight; ++iLoc )
    {
        const Int i = B.GlobalRow(iLoc);
        for( Int j=0; j<numRHS; ++j )
        {
            if( m >= n )
                D.QueueUpdate( i, j, B.GetLocal(iLoc,j) );
            else
                D.QueueUpdate( i+n, j, B.GetLocal(iLoc,j) );
        }
    }      
    D.ProcessQueues();

    // Compute the regularized quasi-semidefinite fact of J
    // ====================================================
    DistSparseMatrix<F> JOrig(comm);
    JOrig = J;
    UpdateRealPartOfDiagonal( J, Real(1), reg, 0, true );

    DistMap map, invMap;
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    if( commRank == 0 && time )
        timer.Start();
    ldl::NestedDissection( J.LockedDistGraph(), map, rootSep, info );
    if( commRank == 0 && time )
        cout << "  ND: " << timer.Stop() << " secs" << endl;
    InvertMap( map, invMap );
    ldl::DistFront<F> JFront( J, map, rootSep, info );

    if( commRank == 0 && time )
        timer.Start();
    LDL( info, JFront, LDL_2D );
    if( commRank == 0 && time )
        cout << "  LDL: " << timer.Stop() << " secs" << endl;

    // Successively solve each of the k linear systems
    // ===============================================
    // TODO: Extend the iterative refinement to handle multiple right-hand sides
    DistMultiVec<F> u(comm);
    Zeros( u, m+n, 1 );
    auto& DLoc = D.Matrix();
    auto& uLoc = u.Matrix();
    if( commRank == 0 && time )
        timer.Start();
    for( Int j=0; j<numRHS; ++j )
    {
        auto dLoc = DLoc( ALL, IR(j) );
        Copy( dLoc, uLoc );
        ldl::SolveAfter( JOrig, reg, invMap, info, JFront, u, ctrl );
        Copy( uLoc, dLoc );
    }
    if( commRank == 0 && time )
        cout << "  Solve: " << timer.Stop() << " secs" << endl;

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
    DEBUG_ONLY(CSE cse("LeastSquares"))
    typedef Base<F> Real;
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);

    DistSparseMatrix<F> ABar(comm);
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
    DistMultiVec<Real> dR(comm), dC(comm);
    if( ctrl.equilibrate )
    {
        auto normMap = []( Real beta ) 
          { return ( beta < Sqrt(Epsilon<Real>()) ? Real(1) : beta ); };
        if( m >= n )
        {
            ColumnTwoNorms( ABar, dC );
            EntrywiseMap( dC, function<Real(Real)>(normMap) );
            DiagonalSolve( RIGHT, NORMAL, dC, ABar );
            Ones( dR, m, 1 );
        }
        else
        {
            RowTwoNorms( ABar, dR );
            EntrywiseMap( dR, function<Real(Real)>(normMap) );
            DiagonalSolve( LEFT, NORMAL, dR, ABar );
            Ones( dC, n, 1 );
        }
    }
    else
    {
        Ones( dR, m, 1 );
        Ones( dC, n, 1 );
    } 
    Real normScale = 1;
    if( ctrl.scaleTwoNorm )
    {
        // Scale ABar down to roughly unit two-norm
        normScale = TwoNormEstimate( ABar, ctrl.basisSize );
        if( ctrl.progress && commRank == 0 )
            cout << "Estimated || A ||_2 ~= " << normScale << endl;
        ABar *= F(1)/normScale;
        dR *= normScale;
    }

    // Equilibrate the RHS
    // ===================
    DiagonalSolve( LEFT, NORMAL, dR, BBar );

    // Solve the equilibrated least squares problem
    // ============================================
    ls::Equilibrated
    ( ABar, BBar, X, ctrl.alpha, ctrl.damp, ctrl.dampTmp, ctrl.regLDLCtrl, 
      ctrl.time );

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
#include "El/macros/Instantiate.h"

} // namespace El
