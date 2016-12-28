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
namespace socp {
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
//      alpha_p := -min_i min eig(s_i),
//      alpha_d := -min_i min eig(z_i),
//
//    where min eig(s_i) is the minimum (Jordan) eigenvalue of s restricted to
//    its i'th subcone.
//
//    Then shift s and z according to the rules:
//
//      s := ( alpha_p > -sqrt(eps)*Max(1,||s||_2) ? s + (1+alpha_p)e : s )
//      z := ( alpha_d > -sqrt(eps)*Max(1,||z||_2) ? z + (1+alpha_d)e : z ),
//
//    where 'eps' is the machine precision and 'e' is the identity of the
//    product cone.
//
//    TODO(poulson):
//    Since the post-processing in step (3) has a large discontinuity as the
//    minimum entry approaches sqrt(eps)*Max(1,||q||_2), we can also provide
//    the ability to instead use an per-subcone lower clip.
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
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    if( primalInit )
    {
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( s.Height() != k || s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != k || z.Width() != 1 )
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
    const bool onlyLower = true;
    lp::affine::KKT( A, G, ones, ones, J, onlyLower );

    // Factor the KKT matrix
    // =====================
    Matrix<Real> dSub;
    Permutation p;
    LDL( J, dSub, p, false );

    // w should be equal to the identity element, and so should sqrt(w)
    Matrix<Real> rc, rb, rh, rmu, u, d;
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | 0 A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;
        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        lp::affine::ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
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
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        lp::affine::ExpandCoreSolution( m, n, k, d, u, y, z );
    }

    const Real eps = limits::Epsilon<Real>();
    const Real sNorm = Nrm2( s );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(eps)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(eps)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : s + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaPrimal = -soc::MinEig( s, orders, firstInds );
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaDual = -soc::MinEig( z, orders, firstInds );
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            soc::Shift( s, alphaPrimal+1, orders, firstInds );
        if( alphaDual >= -gammaDual )
            soc::Shift( z, alphaDual+1, orders, firstInds );
    }
    else
    {
        soc::PushInto( s, orders, firstInds, 1+gammaPrimal );
        soc::PushInto( z, orders, firstInds, 1+gammaDual );
    }
}

template<typename Real>
void Initialize
( const DistMatrix<Real>& A,
  const DistMatrix<Real>& G,
  const DistMatrix<Real>& b,
  const DistMatrix<Real>& c,
  const DistMatrix<Real>& h,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        DistMatrix<Real>& x,
        DistMatrix<Real>& y,
        DistMatrix<Real>& z,
        DistMatrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift,
  Int cutoff )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Grid& g = A.Grid();
    if( primalInit )
    {
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( s.Height() != k || s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != k || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO(poulson): Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    DistMatrix<Real> J(g), ones(g);
    Ones( ones, k, 1 );
    const bool onlyLower = true;
    lp::affine::KKT( A, G, ones, ones, J, onlyLower );

    // Factor the KKT matrix
    // =====================
    DistMatrix<Real> dSub(g);
    DistPermutation p(g);
    LDL( J, dSub, p, false );

    DistMatrix<Real> rc(g), rb(g), rh(g), rmu(g), d(g), u(g);
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | 0 A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;
        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        lp::affine::ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
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
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        lp::affine::ExpandCoreSolution( m, n, k, d, u, y, z );
    }

    const Real eps = limits::Epsilon<Real>();
    const Real sNorm = Nrm2( s );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(eps)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(eps)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : s + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaPrimal = -soc::MinEig( s, orders, firstInds, cutoff );
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaDual = -soc::MinEig( z, orders, firstInds, cutoff );
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            soc::Shift( s, alphaPrimal+1, orders, firstInds );
        if( alphaDual >= -gammaDual )
            soc::Shift( z, alphaDual+1, orders, firstInds );
    }
    else
    {
        soc::PushInto( s, orders, firstInds, 1+gammaPrimal, cutoff );
        soc::PushInto( z, orders, firstInds, 1+gammaDual, cutoff );
    }
}

template<typename Real>
void Initialize
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE

    // TODO(poulson): Expose as control parameters
    const Real eps = limits::Epsilon<Real>();
    const Real gamma = Pow(eps,Real(0.25));
    const Real delta = Pow(eps,Real(0.25));
    const Real beta  = Pow(eps,Real(0.25));
    const Real gammaTmp = 0;
    const Real deltaTmp = 0;
    const Real betaTmp  = 0;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    if( primalInit )
    {
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( s.Height() != k || s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != k || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO(poulson): Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    SparseMatrix<Real> JOrig;
    Matrix<Real> ones;
    Ones( ones, k, 1 );
    const bool onlyLower = false;
    lp::affine::StaticKKT( A, G, gamma, delta, beta, JOrig, onlyLower );
    lp::affine::FinishKKT( m, n, ones, ones, JOrig );
    auto J = JOrig;
    J.FreezeSparsity();

    // (Approximately) factor the KKT matrix
    // =====================================
    Matrix<Real> reg;
    reg.Resize( n+m+k, 1 );
    for( Int i=0; i<reg.Height(); ++i )
    {
        if( i < n )        reg.Set( i, 0, gammaTmp*gammaTmp );
        else if( i < n+m ) reg.Set( i, 0, -deltaTmp*deltaTmp );
        else               reg.Set( i, 0, -betaTmp*betaTmp );
    }
    UpdateRealPartOfDiagonal( J, Real(1), reg );

    SparseLDLFactorization<Real> sparseLDLFact;
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
    sparseLDLFact.Factor();

    Matrix<Real> rc, rb, rh, rmu, u, d;
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | 0 A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;
        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );
        lp::affine::ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
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
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );
        lp::affine::ExpandCoreSolution( m, n, k, d, u, y, z );
    }

    const Real sNorm = Nrm2( s );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(eps)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(eps)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : s + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaPrimal = -soc::MinEig( s, orders, firstInds );
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaDual = -soc::MinEig( z, orders, firstInds );
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            soc::Shift( s, alphaPrimal+1, orders, firstInds );
        if( alphaDual >= -gammaDual )
            soc::Shift( z, alphaDual+1, orders, firstInds );
    }
    else
    {
        soc::PushInto( s, orders, firstInds, 1+gammaPrimal );
        soc::PushInto( z, orders, firstInds, 1+gammaDual );
    }
}

template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  bool primalInit, bool dualInit, bool standardShift,
  Int cutoffPar,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE

    // TODO(poulson): Expose as control parameters
    const Real eps = limits::Epsilon<Real>();
    const Real gamma = Pow(eps,Real(0.25));
    const Real delta = Pow(eps,Real(0.25));
    const Real beta  = Pow(eps,Real(0.25));
    const Real gammaTmp = 0;
    const Real deltaTmp = 0;
    const Real betaTmp  = 0;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Grid& grid = A.Grid();
    if( primalInit )
    {
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( s.Height() != k || s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != k || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO(poulson): Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    DistSparseMatrix<Real> JOrig(grid);
    DistMultiVec<Real> ones(grid);
    Ones( ones, k, 1 );
    const bool onlyLower = false;
    lp::affine::StaticKKT( A, G, gamma, delta, beta, JOrig, onlyLower );
    lp::affine::FinishKKT( m, n, ones, ones, JOrig );
    auto J = JOrig;
    J.FreezeSparsity();

    // (Approximately) factor the KKT matrix
    // =====================================
    DistMultiVec<Real> reg(grid);
    reg.Resize( n+m+k, 1 );
    for( Int iLoc=0; iLoc<reg.LocalHeight(); ++iLoc )
    {
        const Int i = reg.FirstLocalRow() + iLoc;
        if( i < n )        reg.SetLocal( iLoc, 0,  gammaTmp*gammaTmp );
        else if( i < n+m ) reg.SetLocal( iLoc, 0, -deltaTmp*deltaTmp );
        else               reg.SetLocal( iLoc, 0, -betaTmp*betaTmp   );
    }
    UpdateRealPartOfDiagonal( J, Real(1), reg );

    DistSparseLDLFactorization<Real> sparseLDLFact;
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
    sparseLDLFact.Factor( LDL_2D );

    DistMultiVec<Real> rc(grid), rb(grid), rh(grid), rmu(grid), u(grid),
                       d(grid);
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | 0 A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;

        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );
        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );
        lp::affine::ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
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
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );

        lp::affine::KKTRHS( rc, rb, rh, rmu, ones, d );
        reg_ldl::RegularizedSolveAfter
        ( JOrig, reg, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );
        lp::affine::ExpandCoreSolution( m, n, k, d, u, y, z );
    }

    const Real sNorm = Nrm2( s );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(eps)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(eps)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : s + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaPrimal =
          -soc::MinEig( s, orders, firstInds, cutoffPar );
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const Real alphaDual = -soc::MinEig( z, orders, firstInds, cutoffPar );
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            soc::Shift( s, alphaPrimal+1, orders, firstInds );
        if( alphaDual >= -gammaDual )
            soc::Shift( z, alphaDual+1, orders, firstInds );
    }
    else
    {
        soc::PushInto( s, orders, firstInds, 1+gammaPrimal, cutoffPar );
        soc::PushInto( z, orders, firstInds, 1+gammaPrimal, cutoffPar );
    }
}

#define PROTO(Real) \
  template void Initialize \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    bool primalInit, bool dualInit, bool standardShift ); \
  template void Initialize \
  ( const DistMatrix<Real>& A, \
    const DistMatrix<Real>& G, \
    const DistMatrix<Real>& b, \
    const DistMatrix<Real>& c, \
    const DistMatrix<Real>& h, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
          DistMatrix<Real>& x, \
          DistMatrix<Real>& y, \
          DistMatrix<Real>& z, \
          DistMatrix<Real>& s, \
    bool primalInit, bool dualInit, bool standardShift, Int cutoff ); \
  template void Initialize \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    bool primalInit, bool dualInit, bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl ); \
  template void Initialize \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMultiVec<Real>& s, \
    bool primalInit, bool dualInit, bool standardShift, Int cutoffPar,\
    const RegSolveCtrl<Real>& solveCtrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace affine
} // namespace socp
} // namespace El
