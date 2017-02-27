/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// Least Absolute Value (LAV) regression minimizes the one norm of the
// residual of a system of equations, i.e.,
//
//   min || A x - b ||_1.
//
// Real instances of the problem are expressable as a Linear Program [1] via
//
//   min 1^T (u + v)
//   s.t. [A, I, -I] [x; u; v] = b, u, v >= 0
//
// which, in affine standard form, becomes
//
//   min [0; 1; 1]^T [x; u; v]
//   s.t. [A, I, -I] | x | = b, | 0 -I  0 | | x | <= | 0 |.
//                   | u |      | 0  0 -I | | u |    | 0 |
//                   | v |                  | v |
//
// [1]
//   A. Charnes, W. W. Cooper, and R. O. Ferguson,
//   "Optimal estimation of executive compensation by linear programming",
//   Management Science, Vol. 1, No. 2, pp. 138--151, 1955.
//

namespace El {

template<typename Real>
void LAV
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> xInd(0,n), uInd(n,n+m), vInd(n+m,n+2*m);

    AffineLPProblem<Matrix<Real>,Matrix<Real>> problem;

    // c := [0;1;1]
    // ============
    Zeros( problem.c, n+2*m, 1 );
    auto cuv = problem.c( IR(n,n+2*m), ALL );
    Fill( cuv, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( problem.A, m, n+2*m );
    auto AHatx = problem.A( IR(0,m), xInd );
    auto AHatu = problem.A( IR(0,m), uInd );
    auto AHatv = problem.A( IR(0,m), vInd );
    AHatx = A;
    FillDiagonal( AHatu, Real( 1) );
    FillDiagonal( AHatv, Real(-1) );

    // \hat b := b
    // ===========
    problem.b = b;

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( problem.G, 2*m, n+2*m );
    auto Guv = problem.G( IR(0,2*m), IR(n,n+2*m) );
    FillDiagonal( Guv, Real(-1) );

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( problem.h, 2*m, 1 );

    // Solve the affine linear program
    // ===============================
    AffineLPSolution<Matrix<Real>> solution;
    LP( problem, solution, ctrl );

    // Extract x
    // ==========
    x = solution.x( xInd, ALL );
}

template<typename Real>
void LAV
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
        AbstractDistMatrix<Real>& xPre,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixWriteProxy<Real,Real,MC,MR> xProx( xPre );
    auto& x = xProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const Range<Int> xInd(0,n), uInd(n,n+m), vInd(n+m,n+2*m);

    AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>> problem;
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.G.SetGrid( grid );
    problem.h.SetGrid( grid );

    // c := [0;1;1]
    // ============
    Zeros( problem.c, n+2*m, 1 );
    auto cuv = problem.c( IR(n,n+2*m), ALL );
    Fill( cuv, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( problem.A, m, n+2*m );
    auto AHatx = problem.A( IR(0,m), xInd );
    auto AHatu = problem.A( IR(0,m), uInd );
    auto AHatv = problem.A( IR(0,m), vInd );
    AHatx = A;
    FillDiagonal( AHatu, Real( 1) );
    FillDiagonal( AHatv, Real(-1) );

    // \hat b := b
    // ===========
    problem.b = b;

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( problem.G, 2*m, n+2*m );
    auto Guv = problem.G( IR(0,2*m), IR(n,n+2*m) );
    FillDiagonal( Guv, Real(-1) );

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( problem.h, 2*m, 1 );

    // Solve the affine linear program
    // ===============================
    AffineLPSolution<DistMatrix<Real>> solution;
    solution.x.SetGrid( grid );
    solution.s.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
    LP( problem, solution, ctrl );

    // Extract x
    // =========
    x = solution.x( xInd, ALL );
}

template<typename Real>
void LAV
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> xInd(0,n), uInd(n,n+m), vInd(n+m,n+2*m);

    AffineLPProblem<SparseMatrix<Real>,Matrix<Real>> problem;

    // c := [0;1;1]
    // ============
    Zeros( problem.c, n+2*m, 1 );
    auto cuv = problem.c( IR(n,n+2*m), ALL );
    Fill( cuv, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( problem.A, m, n+2*m );
    const Int numEntriesA = A.NumEntries();
    problem.A.Reserve( numEntriesA + 2*m );
    for( Int e=0; e<numEntriesA; ++e )
        problem.A.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int i=0; i<m; ++i )
    {
        problem.A.QueueUpdate( i, i+n,   Real( 1) );
        problem.A.QueueUpdate( i, i+n+m, Real(-1) );
    }
    problem.A.ProcessQueues();

    // \hat b := b
    // ===========
    problem.b = b;

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( problem.G, 2*m, n+2*m );
    problem.G.Reserve( problem.G.Height() );
    for( Int i=0; i<2*m; ++i )
        problem.G.QueueUpdate( i, i+n, Real(-1) );
    problem.G.ProcessQueues();

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( problem.h, 2*m, 1 );

    // Solve the affine linear program
    // ===============================
    AffineLPSolution<Matrix<Real>> solution;
    LP( problem, solution, ctrl );

    // Extract x
    // =========
    x = solution.x( xInd, ALL );
}

template<typename Real>
void LAV
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();

    AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>> problem;
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.G.SetGrid( grid );
    problem.h.SetGrid( grid );

    // c := [0;1;1]
    // ============
    Zeros( problem.c, n+2*m, 1 );
    for( Int iLoc=0; iLoc<problem.c.LocalHeight(); ++iLoc )
        if( problem.c.GlobalRow(iLoc) >= n )
            problem.c.SetLocal( iLoc, 0, Real(1) );

    // \hat A := [A, I, -I]
    // ====================
    Zeros( problem.A, m, n+2*m );
    const Int numLocalEntriesA = A.NumLocalEntries();
    problem.A.Reserve( numLocalEntriesA + 2*problem.A.LocalHeight() );
    for( Int e=0; e<numLocalEntriesA; ++e )
        problem.A.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int iLoc=0; iLoc<problem.A.LocalHeight(); ++iLoc )
    {
        const Int i = problem.A.GlobalRow(iLoc);
        problem.A.QueueLocalUpdate( iLoc, i+n,   Real( 1) );
        problem.A.QueueLocalUpdate( iLoc, i+n+m, Real(-1) );
    }
    problem.A.ProcessLocalQueues();

    // \hat b := b
    // ===========
    problem.b = b;

    // G := | 0 -I  0 |
    //      | 0  0 -I |
    // ================
    Zeros( problem.G, 2*m, n+2*m );
    problem.G.Reserve( problem.G.LocalHeight() );
    for( Int iLoc=0; iLoc<problem.G.LocalHeight(); ++iLoc )
        problem.G.QueueLocalUpdate
        ( iLoc, problem.G.GlobalRow(iLoc)+n, Real(-1) );
    problem.G.ProcessLocalQueues();

    // h := | 0 |
    //      | 0 |
    // ==========
    Zeros( problem.h, 2*m, 1 );

    // Solve the affine QP
    // ===================
    AffineLPSolution<DistMultiVec<Real>> solution;
    solution.x.SetGrid( grid );
    solution.s.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
    LP( problem, solution, ctrl );

    // Extract x
    // =========
    x = solution.x( IR(0,n), ALL );
}

#define PROTO(Real) \
  template void LAV \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LAV \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LAV \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LAV \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
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
