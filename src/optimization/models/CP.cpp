/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// A Chebyshev point (CP) minimizes the supremum norm of A x - b, i.e.,
//
//   min || A x - b ||_oo.
//
// Real instances of the problem are expressable as a Linear Program via
//
//   min t
//   s.t. -t <= A x - b <= t,
//
// which, in affine standard form, becomes
//
//   min [0; 1]^T [x; t]
//   s.t. |  A  -1 | | x | <= |  b |
//        | -A  -1 | | t |    | -b |
//
// NOTE: There is likely an appropriate citation, but the derivation is
//       trivial. If one is found, it will be added.

namespace El {

template<typename Real>
void CP
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    AffineLPProblem<Matrix<Real>,Matrix<Real>> problem;

    // c := [zeros(n,1);1]
    // ===================
    Zeros( problem.c, n+1, 1 );
    problem.c(n) = Real(1);

    // \hat A := zeros(0,n+1)
    // ======================
    Zeros( problem.A, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( problem.b, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( problem.G, 2*m, n+1 );
    auto GTL = problem.G( IR(0,m),   IR(0,n)   );
    auto GBL = problem.G( IR(m,2*m), IR(0,n)   );
    auto gR =  problem.G( IR(0,2*m), IR(n,n+1) );
    GTL = A;
    GBL -= GTL;
    Fill( gR, Real(-1) );

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( problem.h, 2*m, 1 );
    auto hT = problem.h( IR(0,m),   ALL );
    auto hB = problem.h( IR(m,2*m), ALL );
    hT = b;
    hB -= hT;

    // Solve the affine LP
    // ===================
    AffineLPSolution<Matrix<Real>> solution;
    LP( problem, solution, ctrl );

    // Extract x from [x;t]
    // ====================
    x = solution.x( IR(0,n), ALL );
}

template<typename Real>
void CP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>> problem;
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.G.SetGrid( grid );
    problem.h.SetGrid( grid );

    // c := [zeros(n,1);1]
    // ===================
    Zeros( problem.c, n+1, 1 );
    problem.c.Set( n, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ======================
    Zeros( problem.A, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( problem.b, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( problem.G, 2*m, n+1 );
    auto GTL = problem.G( IR(0,m),   IR(0,n)   );
    auto GBL = problem.G( IR(m,2*m), IR(0,n)   );
    auto gR =  problem.G( IR(0,2*m), IR(n,n+1) );
    GTL = A;
    GBL -= GTL;
    Fill( gR, Real(-1) );

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( problem.h, 2*m, 1 );
    auto hT = problem.h( IR(0,m),   ALL );
    auto hB = problem.h( IR(m,2*m), ALL );
    hT = b;
    hB -= hT;

    // Solve the affine LP
    // ===================
    AffineLPSolution<DistMatrix<Real>> solution;
    LP( problem, solution, ctrl );

    // Extract x from [x;t]
    // ====================
    x = solution.x( IR(0,n), ALL );
}

template<typename Real>
void CP
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    AffineLPProblem<SparseMatrix<Real>,Matrix<Real>> problem;

    // c := [zeros(n,1);1]
    // ===================
    Zeros( problem.c, n+1, 1 );
    problem.c(n) = Real(1);

    // \hat A := zeros(0,n+1)
    // ======================
    Zeros( problem.A, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( problem.b, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( problem.G, 2*m, n+1 );
    const Int numEntriesA = A.NumEntries();
    problem.G.Reserve( 2*numEntriesA + 2*m );
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Real value = A.Value(e);
        problem.G.QueueUpdate( i,   j,  value );
        problem.G.QueueUpdate( i+m, j, -value );
    }
    for( Int i=0; i<2*m; ++i )
        problem.G.QueueUpdate( i, n, Real(-1) );
    problem.G.ProcessQueues();

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( problem.h, 2*m, 1 );
    auto hT = problem.h( IR(0,m),   ALL );
    auto hB = problem.h( IR(m,2*m), ALL );
    hT = b;
    hB -= hT;

    // Solve the affine LP
    // ===================
    AffineLPSolution<Matrix<Real>> solution;
    LP( problem, solution, ctrl );

    // Extract x from [x;t]
    // ====================
    x = solution.x( IR(0,n), ALL );
}

template<typename Real>
void CP
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

    // c := [zeros(n,1);1]
    // ===================
    Zeros( problem.c, n+1, 1 );
    problem.c.Set( n, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ======================
    Zeros( problem.A, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( problem.b, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( problem.G, 2*m, n+1 );
    problem.G.Reserve
    ( 2*A.NumLocalEntries()+problem.G.LocalHeight(), 2*A.NumLocalEntries() );
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        problem.G.QueueUpdate( A.Row(e),   A.Col(e),  A.Value(e) );
        problem.G.QueueUpdate( A.Row(e)+m, A.Col(e), -A.Value(e) );
    }
    for( Int iLoc=0; iLoc<problem.G.LocalHeight(); ++iLoc )
        problem.G.QueueLocalUpdate( iLoc, n, Real(-1) );
    problem.G.ProcessQueues();

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( problem.h, 2*m, 1 );
    problem.h.Reserve( 2*b.LocalHeight() );
    for( Int iLoc=0; iLoc<b.LocalHeight(); ++iLoc )
    {
        const Int i = b.GlobalRow(iLoc);
        problem.h.QueueUpdate( i,   0,  b.GetLocal(iLoc,0) );
        problem.h.QueueUpdate( i+m, 0, -b.GetLocal(iLoc,0) );
    }
    problem.h.ProcessQueues();

    // Solve the affine LP
    // ===================
    AffineLPSolution<DistMultiVec<Real>> solution;
    LP( problem, solution, ctrl );

    // Extract x from [x;t]
    // ====================
    x = solution.x( IR(0,n), ALL );
}

#define PROTO(Real) \
  template void CP \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
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
