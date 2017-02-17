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
namespace affine {

//
// Despite the fact that the CVXOPT documentation [1] suggests a single-stage
// procedure for initializing (x,y,z,s), a post-processed two-stage procedure
// is currently used by the code [2]:
//
// 1) Minimize || G x - h ||^2, s.t. A x = b  by solving
//
//    | 0 A^T G^T | |  x |   | 0 |
//    | A  0   0  | |  u | = | b |,
//    | G  0  -I  | | -s |   | h |
//
//   where 'u' is an unused dummy variable.
//
// 2) Minimize || z ||^2, s.t. A^T y + G^T z + c = 0 by solving
//
//    | 0 A^T G^T | | u |   | -c |
//    | A  0   0  | | y | = |  0 |,
//    | G  0  -I  | | z |   |  0 |
//
//    where 'u' is an unused dummy variable.
//
// 3) Set
//
//      alpha_p := -min(s), and
//      alpha_d := -min(z).
//
//    Then shift s and z according to the rules:
//
//      s := ( alpha_p > -sqrt(eps)*Max(1,||s||_2) ? s + (1+alpha_p)e : s )
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
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  bool primalInit,
  bool dualInit,
  bool standardShift )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    if( primalInit )
    {
        if( solution.x.Height() != n || solution.x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( solution.s.Height() != k || solution.s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInit )
    {
        if( solution.y.Height() != m || solution.y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( solution.z.Height() != k || solution.z.Width() != 1 )
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
    Ones( ones, k, 1 );
    KKT( problem.A, problem.G, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    Matrix<Real> dSub;
    Permutation p;
    LDL( J, dSub, p, false );

    AffineLPResidual<Matrix<Real>> residual;
    Matrix<Real> u, d;
    Zeros( residual.dualConic, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | 0 A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( residual.dualEquality, n, 1 );
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, solution.x, u, solution.s );
        solution.s *= -1;
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y + G^T z + c = 0 by solving
        //
        //    | 0 A^T G^T | | u |   | -c |
        //    | A  0   0  | | y | = |  0 |,
        //    | G  0  -I  | | z |   |  0 |
        //
        //    where 'u' is an unused dummy variable.
        residual.dualEquality = problem.c;
        Zeros( residual.primalEquality, m, 1 );
        Zeros( residual.primalConic, k, 1 );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, u, solution.y, solution.z );
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real sNorm = Nrm2( solution.s );
    const Real zNorm = Nrm2( solution.z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : s + alpha*e >= 0 }
        // -------------------------------------------
        const auto sMinPair = VectorMinLoc( solution.s );
        const Real alphaPrimal = -sMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( solution.z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( solution.s, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( solution.z, alphaDual+1 );
    }
    else
    {
        LowerClip( solution.s, gammaPrimal );
        LowerClip( solution.z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  bool primalInit,
  bool dualInit,
  bool standardShift )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Grid& grid = problem.A.Grid();
    if( primalInit )
    {
        if( solution.x.Height() != n || solution.x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( solution.s.Height() != k || solution.s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInit )
    {
        if( solution.y.Height() != m || solution.y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( solution.z.Height() != k || solution.z.Width() != 1 )
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
    Ones( ones, k, 1 );
    KKT( problem.A, problem.G, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    DistMatrix<Real> dSub(grid);
    DistPermutation p(grid);
    LDL( J, dSub, p, false );

    AffineLPResidual<DistMatrix<Real>> residual;
    residual.primalEquality.SetGrid(grid);
    residual.primalConic.SetGrid(grid);
    residual.dualEquality.SetGrid(grid);
    residual.dualConic.SetGrid(grid);
    Zeros( residual.dualConic, k, 1 );
    DistMatrix<Real> d(grid), u(grid);
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | 0 A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( residual.dualEquality, n, 1 );
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, solution.x, u, solution.s );
        solution.s *= -1;
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y + G^T z + c = 0 by solving
        //
        //    | 0 A^T G^T | | u |   | -c |
        //    | A  0   0  | | y | = |  0 |,
        //    | G  0  -I  | | z |   |  0 |
        //
        //    where 'u' is an unused dummy variable.
        residual.dualEquality = problem.c;
        Zeros( residual.primalEquality, m, 1 );
        Zeros( residual.primalConic, k, 1 );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, u, solution.y, solution.z );
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real sNorm = Nrm2( solution.s );
    const Real zNorm = Nrm2( solution.z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : s + alpha*e >= 0 }
        // -------------------------------------------
        const auto sMinPair = VectorMinLoc( solution.s );
        const Real alphaPrimal = -sMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( solution.z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( solution.s, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( solution.z, alphaDual+1 );
    }
    else
    {
        LowerClip( solution.s, gammaPrimal );
        LowerClip( solution.z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const SparseMatrix<Real>& JStatic,
  const Matrix<Real>& regTmp,
        SparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    qp::affine::Initialize
    ( JStatic, regTmp,
      problem.b, problem.c, problem.h,
      solution.x, solution.y, solution.z, solution.s,
      sparseLDLFact,
      primalInit, dualInit, standardShift, solveCtrl );
}

template<typename Real>
void Initialize
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const DistSparseMatrix<Real>& JStatic,
  const DistMultiVec<Real>& regTmp,
        DistSparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    qp::affine::Initialize
    ( JStatic, regTmp,
      problem.b, problem.c, problem.h,
      solution.x, solution.y, solution.z, solution.s,
      sparseLDLFact,
      primalInit, dualInit, standardShift, solveCtrl );
}

#define PROTO(Real) \
  template void Initialize \
  ( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    bool primalInit, \
    bool dualInit, \
    bool standardShift ); \
  template void Initialize \
  ( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          AffineLPSolution<DistMatrix<Real>>& solution, \
    bool primalInit, \
    bool dualInit, \
    bool standardShift ); \
  template void Initialize \
  ( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    const SparseMatrix<Real>& JStatic, \
    const Matrix<Real>& regTmp, \
          SparseLDLFactorization<Real>& sparseLDLFact, \
    bool primalInit, \
    bool dualInit, \
    bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl ); \
  template void Initialize \
  ( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          AffineLPSolution<DistMultiVec<Real>>& solution, \
    const DistSparseMatrix<Real>& JStatic, \
    const DistMultiVec<Real>& regTmp, \
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

} // namespace affine
} // namespace lp
} // namespace El
