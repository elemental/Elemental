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
( const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift )
{
    DEBUG_ONLY(CSE cse("lp::direct::Initialize"))
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
    AugmentedKKT( A, ones, ones, J );

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
        //    | I A^T | | x |   | 0 |
        //    | A  0  | | u | = | b |,
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
        // Minimize || z ||^2, s.t. A^T y - z + c = 0 by solving
        //
        //    | I A^T | | -z |   | -c |
        //    | A  0  | |  y | = |  0 |.
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
( const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& b, const ElementalMatrix<Real>& c,
        ElementalMatrix<Real>& x,       ElementalMatrix<Real>& y,
        ElementalMatrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift )
{
    DEBUG_ONLY(CSE cse("lp::direct::Initialize"))
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
    AugmentedKKT( A, ones, ones, J );

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
        //    | I A^T | | x |   | 0 |
        //    | A  0  | | u | = | b |,
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
        // Minimize || z ||^2, s.t. A^T y - z + c = 0 by solving
        //
        //    | I A^T | | -z |   | -c |
        //    | A  0  | |  y | = |  0 |.
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
( const SparseMatrix<Real>& A, 
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
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::direct::Initialize
    ( Q, A, b, c, x, y, z, map, invMap, rootSep, info,
      primalInit, dualInit, standardShift, solveCtrl );
}

template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A, 
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
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> Q(comm);
    Q.Resize( n, n );
    qp::direct::Initialize
    ( Q, A, b, c, x, y, z, map, invMap, rootSep, info, 
      mappedSources, mappedTargets, colOffs,
      primalInit, dualInit, standardShift, solveCtrl );
}

#define PROTO(Real) \
  template void Initialize \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    bool primalInit, bool dualInit, bool standardShift ); \
  template void Initialize \
  ( const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
    bool primalInit, bool dualInit, bool standardShift ); \
  template void Initialize \
  ( const SparseMatrix<Real>& A, \
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
  ( const DistSparseMatrix<Real>& A, \
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
} // namespace lp
} // namespace El
