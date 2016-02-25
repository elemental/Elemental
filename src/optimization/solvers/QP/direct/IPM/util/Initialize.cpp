/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "../util.hpp"

namespace El {
namespace qp {
namespace direct {

//
// Despite the fact that the CVXOPT documentation [1] suggests a single-stage
// procedure for initializing (x,y,z,s), a post-processed two-stage procedure 
// is currently used by the code [2], which, in the case that G = -I and h = 0,
//
// 1) Minimize || x ||^2, s.t. A x = b  by solving
//
//    | Q A^T -I | |  x |   | 0 |
//    | A  0   0 | |  u | = | b |,
//    |-I  0  -I | | -s |   | 0 |
//
//   where 'u' is an unused dummy variable. A Schur-complement manipulation
//   yields
//
//    | Q+I A^T | | x |   | 0 |
//    | A    0  | | u | = | b |.
//
// 2) Minimize || z ||^2, s.t. A^T y - z + c in range(Q) by solving
//
//    | Q A^T -I | | u |   | -c |
//    | A  0   0 | | y | = |  0 |,
//    |-I  0  -I | | z |   |  0 |
//
//    where 'u = -z' is an unused dummy variable. A Schur-complement 
//    manipulation yields
//
//    | Q+I A^T | | -z |   | -c |
//    | A    0  | |  y | = |  0 |.
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
( const Matrix<Real>& Q,
  const Matrix<Real>& A, 
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift )
{
    DEBUG_ONLY(CSE cse("qp::direct::Initialize"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( primalInit )
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != n || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO: Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    Matrix<Real> J, ones;
    Ones( ones, n, 1 );
    AugmentedKKT( Q, A, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    Matrix<Real> dSub;
    Permutation p;
    LDL( J, dSub, p, false );

    Matrix<Real> rc, rb, rmu, u, v, d;
    Zeros( rmu, n, 1 );
    if( !primalInit )
    {
        // Minimize || x ||^2, s.t. A x = b  by solving
        //
        //    | Q+I A^T | | x |   | 0 |
        //    | A    0  | | u | = | b |,
        //
        // where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        Zeros( rmu, n, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution( ones, ones, rmu, d, x, u, v );
    }
    if( !dualInit ) 
    {
        // Minimize || z ||^2, s.t. A^T y - z + c in range(Q) by solving
        //
        //    | Q+I A^T | | -z |   | -c |
        //    | A    0  | |  y | = |  0 |.
        rc = c;
        Zeros( rb, m, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution( ones, ones, rmu, d, z, y, u );
        z *= -1;
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real xNorm = Nrm2( x );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(xNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : x + alpha*e >= 0 }
        // -------------------------------------------
        const auto xMinPair = VectorMinLoc( x );
        const Real alphaPrimal = -xMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized x was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( x, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( x, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& b,
  const ElementalMatrix<Real>& c,
        ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& y,
        ElementalMatrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift )
{
    DEBUG_ONLY(CSE cse("qp::direct::Initialize"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    if( primalInit )
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != n || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO: Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    DistMatrix<Real> J(g), ones(g);
    Ones( ones, n, 1 );
    AugmentedKKT( Q, A, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    DistMatrix<Real> dSub(g);
    DistPermutation p(g);
    LDL( J, dSub, p, false );

    DistMatrix<Real> rc(g), rb(g), rmu(g), u(g), v(g), d(g);
    Zeros( rmu, n, 1 );
    if( !primalInit )
    {
        // Minimize || x ||^2, s.t. A x = b  by solving
        //
        //    | Q+I A^T | | x |   | 0 |
        //    | A    0  | | u | = | b |,
        //
        // where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        Zeros( rmu, n, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution( ones, ones, rmu, d, x, u, v );
    }
    if( !dualInit ) 
    {
        // Minimize || z ||^2, s.t. A^T y - z + c in range(Q) by solving
        //
        //    | Q+I A^T | | -z |   | -c |
        //    | A    0  | |  y | = |  0 |.
        rc = c;
        Zeros( rb, m, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandAugmentedSolution( ones, ones, rmu, d, z, y, u );
        z *= -1;
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real xNorm = Nrm2( x );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(xNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : x + alpha*e >= 0 }
        // -------------------------------------------
        const auto xMinPair = VectorMinLoc( x );
        const Real alphaPrimal = -xMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized x was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( x, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( x, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        vector<Int>& map,
        vector<Int>& invMap, 
        ldl::Separator& rootSep,
        ldl::NodeInfo& info,
  bool primalInit, bool dualInit, bool standardShift, 
  const RegSolveCtrl<Real>& solveCtrl )
{
    DEBUG_ONLY(CSE cse("lp::direct::Initialize"))
    const Int m = A.Height();
    const Int n = A.Width();

    const Real eps = limits::Epsilon<Real>();
    const Real gamma = Pow(eps,Real(0.25));
    const Real delta = Pow(eps,Real(0.25));
    const Real gammaTmp = 0;
    const Real deltaTmp = 0;

    if( primalInit )
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != n || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO: Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    SparseMatrix<Real> J, JOrig;
    Matrix<Real> ones;
    Ones( ones, n, 1 );
    AugmentedKKT( Q, A, gamma, delta, ones, ones, JOrig, false );
    J = JOrig;

    // (Approximately) factor the KKT matrix
    // =====================================
    Matrix<Real> reg;
    reg.Resize( n+m, 1 );
    for( Int i=0; i<n+m; ++i )
    {
        if( i < n )
            reg.Set( i, 0, gammaTmp*gammaTmp );
        else
            reg.Set( i, 0, -deltaTmp*deltaTmp );
    }
    UpdateRealPartOfDiagonal( J, Real(1), reg );

    NestedDissection( J.LockedGraph(), map, rootSep, info );
    InvertMap( map, invMap );

    ldl::Front<Real> JFront;
    JFront.Pull( J, map, info );
    LDL( info, JFront, LDL_2D );

    // Compute the proposed step from the KKT system
    // ---------------------------------------------
    Matrix<Real> rc, rb, rmu, d, u, v;
    Zeros( rmu, n, 1 );
    if( !primalInit )
    {
        // Minimize || x ||^2, s.t. A x = b  by solving
        //
        //    | Q+I A^T | | x |   | 0 |
        //    | A    0  | | u | = | b |,
        //
        // where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        Zeros( rmu, n, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, invMap, info, JFront, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandAugmentedSolution( ones, ones, rmu, d, x, u, v );
    }
    if( !dualInit ) 
    {
        // Minimize || z ||^2, s.t. A^T y - z + c in range(Q) by solving
        //
        //    | Q+I A^T | | -z |   | -c |
        //    | A    0  | |  y | = |  0 |.
        rc = c;
        Zeros( rb, m, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, invMap, info, JFront, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandAugmentedSolution( ones, ones, rmu, d, z, y, u );
        z *= -1;
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real xNorm = Nrm2( x );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(xNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : x + alpha*e >= 0 }
        // -------------------------------------------
        const auto xMinPair = VectorMinLoc( x );
        const Real alphaPrimal = -xMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized x was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( x, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( x, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMap& map,
        DistMap& invMap, 
        ldl::DistSeparator& rootSep,
        ldl::DistNodeInfo& info,
        vector<Int>& mappedSources,
        vector<Int>& mappedTargets,
        vector<Int>& colOffs,
  bool primalInit, bool dualInit, bool standardShift, 
  const RegSolveCtrl<Real>& solveCtrl )
{
    DEBUG_ONLY(CSE cse("lp::direct::Initialize"))
    const Int m = A.Height();
    const Int n = A.Width();

    const Real eps = limits::Epsilon<Real>();
    const Real gamma = Pow(eps,Real(0.25));
    const Real delta = Pow(eps,Real(0.25));
    const Real gammaTmp = 0;
    const Real deltaTmp = 0;

    mpi::Comm comm = A.Comm();
    if( primalInit )
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != n || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO: Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    DistSparseMatrix<Real> J(comm), JOrig(comm);
    DistMultiVec<Real> ones(comm);
    Ones( ones, n, 1 );
    AugmentedKKT( Q, A, gamma, delta, ones, ones, JOrig, false );
    J = JOrig;

    // (Approximately) factor the KKT matrix
    // =====================================
    DistMultiVec<Real> reg(comm);
    reg.Resize( n+m, 1 );
    for( Int iLoc=0; iLoc<reg.LocalHeight(); ++iLoc )
    {
        const Int i = reg.GlobalRow(iLoc);
        if( i < n )
            reg.SetLocal( iLoc, 0, gammaTmp*gammaTmp );
        else
            reg.SetLocal( iLoc, 0, -deltaTmp*deltaTmp );
    }
    UpdateRealPartOfDiagonal( J, Real(1), reg );

    NestedDissection( J.LockedDistGraph(), map, rootSep, info );
    InvertMap( map, invMap );

    ldl::DistFront<Real> JFront;
    JFront.Pull( J, map, rootSep, info, mappedSources, mappedTargets, colOffs );
    LDL( info, JFront, LDL_2D );

    // Compute the proposed step from the KKT system
    // ---------------------------------------------
    DistMultiVec<Real> rc(comm), rb(comm), rmu(comm), d(comm), u(comm), v(comm);
    Zeros( rmu, n, 1 );
    if( !primalInit )
    {
        // Minimize || x ||^2, s.t. A x = b  by solving
        //
        //    | Q+I A^T | | x |   | 0 |
        //    | A    0  | | u | = | b |,
        //
        // where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        Zeros( rmu, n, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, invMap, info, JFront, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandAugmentedSolution( ones, ones, rmu, d, x, u, v );
    }
    if( !dualInit ) 
    {
        // Minimize || z ||^2, s.t. A^T y - z + c in range(Q) by solving
        //
        //    | Q+I A^T | | -z |   | -c |
        //    | A    0  | |  y | = |  0 |.
        rc = c;
        Zeros( rb, m, 1 );
        AugmentedKKTRHS( ones, rc, rb, rmu, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, invMap, info, JFront, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandAugmentedSolution( ones, ones, rmu, d, z, y, u );
        z *= -1;
    }

    const Real epsilon = limits::Epsilon<Real>();
    const Real xNorm = Nrm2( x );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(xNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : x + alpha*e >= 0 }
        // -------------------------------------------
        const auto xMinPair = VectorMinLoc( x );
        const Real alphaPrimal = -xMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized x was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( x, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( x, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

#define PROTO(Real) \
  template void Initialize \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    bool primalInit, bool dualInit, bool standardShift ); \
  template void Initialize \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
    bool primalInit, bool dualInit, bool standardShift ); \
  template void Initialize \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          vector<Int>& map, \
          vector<Int>& invMap, \
          ldl::Separator& rootSep, \
          ldl::NodeInfo& info, \
    bool primalInit, bool dualInit, bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl ); \
  template void Initialize \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMap& map, \
          DistMap& invMap, \
          ldl::DistSeparator& rootSep, \
          ldl::DistNodeInfo& info, \
          vector<Int>& mappedSources, \
          vector<Int>& mappedTargets, \
          vector<Int>& colOffs, \
    bool primalInit, bool dualInit, bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace qp
} // namespace El
