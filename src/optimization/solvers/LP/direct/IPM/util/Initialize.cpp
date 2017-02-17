/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "../util.hpp"

namespace El {
namespace lp {
namespace direct {

//
// Despite the fact that the CVXOPT documentation [1] suggests a single-stage
// procedure for initializing (x,y,z,s), a post-processed two-stage procedure
// is currently used by the code [2], which, in the case that G = -I and h = 0,
//
// 1) Minimize || x ||^2, s.t. A x = b  by solving
//
//    | 0 A^T -I | |  x |   | 0 |
//    | A  0   0 | |  u | = | b |,
//    |-I  0  -I | | -s |   | 0 |
//
//   where 'u' is an unused dummy variable. A Schur-complement manipulation
//   yields
//
//    | I A^T | | x |   | 0 |
//    | A  0  | | u | = | b |.
//
// 2) Minimize || z ||^2, s.t. A^T y - z + c = 0 by solving
//
//    | 0 A^T -I | | u |   | -c |
//    | A  0   0 | | y | = |  0 |,
//    |-I  0  -I | | z |   |  0 |
//
//    where 'u = -z' is an unused dummy variable. A Schur-complement
//    manipulation yields
//
//    | I A^T | | -z |   | -c |
//    | A  0  | |  y | = |  0 |.
//
// 3) Set
//
//      alpha_p := -min(x), and
//      alpha_d := -min(z).
//
//    Then shift x and z according to the rules:
//
//      x := ( alpha_p > -sqrt(eps)*Max(1,||x||_2) ? x + (1+alpha_p)e : x )
//      z := ( alpha_d > -sqrt(eps)*Max(1,||z||_2) ? z + (1+alpha_d)e : z ),
//
//    where 'eps' is the machine precision, 'e' is a vector of all ones
//    (for more general conic optimization problems, it is the product of
//    identity elements from the Jordan algebras whose squares yield the
//    relevant cone.
//
//    Since the post-processing in step (3) has a large discontinuity as the
//    minimum entry approaches sqrt(eps)*Max(1,||q||_2), we also provide
//    the ability to instead use an entrywise lower clip.
//
// [1] L. Vandenberghe
//     "The CVXOPT linear and quadratic cone program solvers"
//     <http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf>
//
// [2] L. Vandenberghe
//     CVXOPT's source file, "src/python/coneprog.py"
//     <https://github.com/cvxopt/cvxopt/blob/f3ca94fb997979a54b913f95b816132f7fd44820/src/python/coneprog.py>
//

template<typename Real>
void Initialize
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  bool primalInit,
  bool dualInit,
  bool standardShift )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    if( primalInit )
        if( solution.x.Height() != n || solution.x.Width() != 1 )
            LogicError("x was of the wrong size");
    if( dualInit )
    {
        if( solution.y.Height() != m || solution.y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( solution.z.Height() != n || solution.z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO(poulson): Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    Matrix<Real> J, ones;
    Ones( ones, n, 1 );
    AugmentedKKT( problem.A, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    Matrix<Real> dSub;
    Permutation p;
    LDL( J, dSub, p, false );

    DirectLPResidual<Matrix<Real>> residual;
    Matrix<Real> u, v, d;
    Zeros( residual.dualConic, n, 1 );
    if( !primalInit )
    {
        // Minimize || x ||^2, s.t. A x = b  by solving
        //
        //    | I A^T | | x |   | 0 |
        //    | A  0  | | u | = | b |,
        //
        // where 'u' is an unused dummy variable.
        Zeros( residual.dualEquality, n, 1 );
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Zeros( residual.dualConic, n, 1 );
        AugmentedKKTRHS
        ( ones, residual.dualEquality, residual.primalEquality,
          residual.dualConic, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution
        ( ones, ones, residual.dualConic, d, solution.x, u, v );
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y - z + c = 0 by solving
        //
        //    | I A^T | | -z |   | -c |
        //    | A  0  | |  y | = |  0 |.
        residual.dualEquality = problem.c;
        Zeros( residual.primalEquality, m, 1 );
        AugmentedKKTRHS
        ( ones, residual.dualEquality, residual.primalEquality,
          residual.dualConic, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution
        ( ones, ones, residual.dualConic, d, solution.z, solution.y, u );
        solution.z *= -1;
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real xNorm = Nrm2( solution.x );
    const Real zNorm = Nrm2( solution.z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(xNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : x + alpha*e >= 0 }
        // -------------------------------------------
        const auto xMinPair = VectorMinLoc( solution.x );
        const Real alphaPrimal = -xMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized x was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( solution.z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( solution.x, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( solution.z, alphaDual+1 );
    }
    else
    {
        LowerClip( solution.x, gammaPrimal );
        LowerClip( solution.z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  bool primalInit,
  bool dualInit,
  bool standardShift )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Grid& grid = problem.A.Grid();
    if( primalInit )
        if( solution.x.Height() != n || solution.x.Width() != 1 )
            LogicError("x was of the wrong size");
    if( dualInit )
    {
        if( solution.y.Height() != m || solution.y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( solution.z.Height() != n || solution.z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO(poulson): Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    DistMatrix<Real> J(grid), ones(grid);
    Ones( ones, n, 1 );
    AugmentedKKT( problem.A, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    DistMatrix<Real> dSub(grid);
    DistPermutation p(grid);
    LDL( J, dSub, p, false );

    DirectLPResidual<DistMatrix<Real>> residual;
    ForceSimpleAlignments( residual, grid );
    Zeros( residual.dualConic, n, 1 );
    DistMatrix<Real> u(grid), v(grid), d(grid);
    if( !primalInit )
    {
        // Minimize || x ||^2, s.t. A x = b  by solving
        //
        //    | I A^T | | x |   | 0 |
        //    | A  0  | | u | = | b |,
        //
        // where 'u' is an unused dummy variable.
        Zeros( residual.dualEquality, n, 1 );
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Zeros( residual.dualConic, n, 1 );
        AugmentedKKTRHS
        ( ones, residual.dualEquality, residual.primalEquality,
          residual.dualConic, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution
        ( ones, ones, residual.dualConic, d, solution.x, u, v );
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y - z + c = 0 by solving
        //
        //    | I A^T | | -z |   | -c |
        //    | A  0  | |  y | = |  0 |.
        residual.dualEquality = problem.c;
        Zeros( residual.primalEquality, m, 1 );
        AugmentedKKTRHS
        ( ones, residual.dualEquality, residual.primalEquality,
          residual.dualConic, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution
        ( ones, ones, residual.dualConic, d, solution.z, solution.y, u );
        solution.z *= -1;
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real xNorm = Nrm2( solution.x );
    const Real zNorm = Nrm2( solution.z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(xNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : x + alpha*e >= 0 }
        // -------------------------------------------
        const auto xMinPair = VectorMinLoc( solution.x );
        const Real alphaPrimal = -xMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized x was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( solution.z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( solution.x, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( solution.z, alphaDual+1 );
    }
    else
    {
        LowerClip( solution.x, gammaPrimal );
        LowerClip( solution.z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
        SparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    const Int n = problem.A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::direct::Initialize
    ( Q, problem.A, problem.b, problem.c, solution.x, solution.y, solution.z,
      sparseLDLFact,
      primalInit, dualInit, standardShift, solveCtrl );
}

template<typename Real>
void Initialize
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
        DistSparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    const Int n = problem.A.Width();
    DistSparseMatrix<Real> Q(problem.A.Grid());
    Q.Resize( n, n );
    qp::direct::Initialize
    ( Q, problem.A, problem.b, problem.c, solution.x, solution.y, solution.z,
      sparseLDLFact,
      primalInit, dualInit, standardShift, solveCtrl );
}

#define PROTO(Real) \
  template void Initialize \
  ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
    bool primalInit, \
    bool dualInit, \
    bool standardShift ); \
  template void Initialize \
  ( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          DirectLPSolution<DistMatrix<Real>>& solution, \
    bool primalInit, \
    bool dualInit, \
    bool standardShift ); \
  template void Initialize \
  ( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
          SparseLDLFactorization<Real>& sparseLDLFact, \
    bool primalInit, \
    bool dualInit, \
    bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl ); \
  template void Initialize \
  ( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          DirectLPSolution<DistMultiVec<Real>>& solution, \
          DistSparseLDLFactorization<Real>& sparseLDLFact, \
    bool primalInit, \
    bool dualInit, \
    bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace direct
} // namespace lp
} // namespace El
