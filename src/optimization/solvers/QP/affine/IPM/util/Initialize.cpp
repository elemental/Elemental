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
namespace qp {
namespace affine {

//
// Despite the fact that the CVXOPT QP documentation [1] suggests a single-stage
// initialization for conic QPs, we make use of an analogue of their two-stage
// LP initialization scheme.
//
// 1) In order to compute the initial primal variables, (x,s), we solve 
//
//      argmin_{x,s}
//      { (1/2) x^T Q x + (1/2) || s ||^2 | A x = b, G x + s = h }
//
//    via seeking an extreme point of the Lagrangian
//
//      L(x,s;u,z) = (1/2) x^T Q x + (1/2) s^T s +
//                   u^T (A x - b) + z^T (G x + s - h),
//
//    which must satisfy:

//      Grad_x L = Q x + A^T u + G^T z = 0,
//      Grad_u L = A x - b = 0,
//      Grad_z L = G x + s - h = 0,
//      Grad_s L = s + z = 0.
//
//    The latter implies that s = -z, so we end up with the system
//
//      | Q A^T G^T | |  x |   | 0 |
//      | A  0   0  | |  u | = | b |,
//      | G  0  -I  | | -s |   | h |
//
//   where 'u' is an unused dummy variable.
//
//   In the case of sparse matrices {Q,A,G}, we instead seek an extreme point
//   of the regularized Lagrangian
//
//     L(x,s;u,z) = (1/2) x^T Q x + (1/2) s^T s +
//                  u^T (A x - b) + z^T (G x + s - h) +
//                  (1/2) gamma_x || x ||_2^2 -
//                  (1/2) gamma_y || u ||_2^2 -
//                  (1/2) gamma_z || z ||_2^2.
//
//   The gradients then take the form
//
//     Grad_x L = Q x + A^T u + G^T z + gamma_x x = 0,
//     Grad_u L = A x - b - gamma_y u = 0,
//     Grad_z L = G x + s - h - gamma_z z = 0,
//     Grad_s L = s + z = 0.
//
//   Then s = -z again, and we have the regularized (symmetric quasi-definite)
//   system
//
//     | (Q + gamma_x I),     A^T,           G^T       | |  x | = | 0 |
//     |        A,        -gamma_y I,         0,       | |  u | = | b |.
//     |        G,             0,     -(gamma_z + 1) I | | -s | = | h |
//
//   Clearly small amounts of regularization on z should not make much of a
//   difference, but the regularization on x and y could be of importance to the
//   solver stability.
// 
// 2) In order to initialize the dual variables, (y,z), if 'Q = 0' then we could
//    solve
//
//      argmin_{y,z} { (1/2) || z ||_2^2 | A^T y + G^T z + c = 0 }
//
//    by seeking an extremal point of the Lagrangian
//
//      L(y,z;u) = (1/2) || z ||_2^2 + u^T (A^T y + G^T z + c).
//
//    It would follow that
//
//      Grad_u L = A^T y + G^T z + c = 0,
//      Grad_y L = A u = 0,
//      Grad_z L = G u + z = 0,
//      
//    which could be rearranged into the form
//
//      | 0, A^T, G^T | | -u | = | -c |.
//      | A,  0,   0  | |  y |   |  0 |
//      | G,  0,  -I  | |  z |   |  0 |
//
//    If Q is nontrivial, then we can instead investigate
//
//      argmin_{w,y,z} { (1/2) w^T Q w + (1/2) || z ||_2^2 |
//                       Q w + A^T y + G^T z + c = 0 },
//
//    by seeking an extremal point of the Lagrangian
//
//      L(y,z;w,u) = (1/2) w^T Q w + (1/2) || z ||_2^2 +
//                   u^T (Q w + A^T y + G^T z + c).
//
//    It follows that
//
//      Grad_u L = Q w + A^T y + G^T z + c = 0,
//      Grad_y L = A u = 0,
//      Grad_z L = G u + z = 0,
//      Grad_w L = Q w + Q u = 0.
//
//    Thus, Q w = - Q u, so we can form the system of equations
//      
//      | Q, A^T, G^T | | -u | = | -c |.
//      | A,  0,   0  | |  y |   |  0 |
//      | G,  0,  -I  | |  z |   |  0 |
//
//    TODO(poulson): Describe regularization strategy.
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

template<typename Real>
void Initialize
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  bool primalInit,
  bool dualInit,
  bool standardShift )
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
    KKT( Q, A, G, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    Matrix<Real> dSub;
    Permutation p;
    LDL( J, dSub, p, false );

    Matrix<Real> rc, rb, rh, rmu, u, d;
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | Q A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y + G^T z + c in range(Q) by solving
        //
        //    | Q A^T G^T | | u |   | -c |
        //    | A  0   0  | | y | = |  0 |,
        //    | G  0  -I  | | z |   |  0 |
        //
        //    where 'u' is an unused dummy variable.
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, u, y, z );
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
        const auto sMinPair = VectorMinLoc( s );
        const Real alphaPrimal = -sMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( s, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( s, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& G,
  const ElementalMatrix<Real>& b,
  const ElementalMatrix<Real>& c,
  const ElementalMatrix<Real>& h,
        ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& y,
        ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& s,
  bool primalInit,
  bool dualInit,
  bool standardShift )
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
    KKT( Q, A, G, ones, ones, J );

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
        //    | Q A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y + G^T z + c in range(Q) by solving
        //
        //    | Q A^T G^T | | u |   | -c |
        //    | A  0   0  | | y | = |  0 |,
        //    | G  0  -I  | | z |   |  0 |
        //
        //    where 'u' is an unused dummy variable.
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, u, y, z );
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
        const auto sMinPair = VectorMinLoc( s );
        const Real alphaPrimal = -sMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( s, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( s, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const SparseMatrix<Real>& JStatic,
  const Matrix<Real>& regLarge,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
        SparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    const Int m = b.Height();
    const Int n = c.Height();
    const Int k = h.Height();
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
    auto JOrig = JStatic;
    JOrig.FreezeSparsity();
    Matrix<Real> ones;
    Ones( ones, k, 1 );
    FinishKKT( m, n, ones, ones, JOrig );
    auto J = JOrig;
    J.FreezeSparsity();
    UpdateRealPartOfDiagonal( J, Real(1), regLarge );

    // Analyze the sparsity pattern of the KKT system
    // ==============================================
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );

    // (Approximately) factor the KKT matrix
    // =====================================
    sparseLDLFact.Factor();

    // Compute the proposed step from the KKT system
    // ---------------------------------------------
    Matrix<Real> rc, rb, rh, rmu, u, d;
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | Q A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;
        KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, regLarge, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y + G^T z + c in range(Q) by solving
        //
        //    | Q A^T G^T | | u |   | -c |
        //    | A  0   0  | | y | = |  0 |,
        //    | G  0  -I  | | z |   |  0 |
        //
        //    where 'u' is an unused dummy variable.
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, regLarge, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandCoreSolution( m, n, k, d, u, y, z );
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
        const auto sMinPair = VectorMinLoc( s );
        const Real alphaPrimal = -sMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( s, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( s, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& JStatic,
  const DistMultiVec<Real>& regLarge,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
        DistSparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    const Int m = b.Height();
    const Int n = c.Height();
    const Int k = h.Height();
    const Grid& grid = JStatic.Grid();
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
    JOrig = JStatic;
    JOrig.FreezeSparsity();
    JOrig.LockedDistGraph().multMeta = JStatic.LockedDistGraph().multMeta;
    DistMultiVec<Real> ones(grid);
    Ones( ones, k, 1 );
    FinishKKT( m, n, ones, ones, JOrig );
    auto J = JOrig;
    J.FreezeSparsity();
    J.LockedDistGraph().multMeta = JStatic.LockedDistGraph().multMeta;
    UpdateRealPartOfDiagonal( J, Real(1), regLarge );

    // Analyze the nonzero pattern
    // ===========================
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );

    // (Approximately) factor the KKT matrix
    // =====================================
    // TODO(poulson): Consider selective inversion
    sparseLDLFact.Factor( LDL_2D );

    DistMultiVec<Real> rc(grid), rb(grid), rh(grid), rmu(grid), u(grid),
                       d(grid);
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | Q A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        rh *= -1;
        KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, regLarge, sparseLDLFact, d,
          solveCtrl.relTol,
          solveCtrl.maxRefineIts,
          solveCtrl.progress );

        ExpandCoreSolution( m, n, k, d, x, u, s );
        s *= -1;
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y + G^T z + c in range(Q) by solving
        //
        //    | Q A^T G^T | | u |   | -c |
        //    | A  0   0  | | y | = |  0 |,
        //    | G  0  -I  | | z |   |  0 |
        //
        //    where 'u' is an unused dummy variable.
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, regLarge, sparseLDLFact, d,
          solveCtrl.relTol,
          solveCtrl.maxRefineIts,
          solveCtrl.progress );

        ExpandCoreSolution( m, n, k, d, u, y, z );
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
        const auto sMinPair = VectorMinLoc( s );
        const Real alphaPrimal = -sMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( s, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( s, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
}

#define PROTO(Real) \
  template void Initialize \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    bool primalInit, bool dualInit, bool standardShift ); \
  template void Initialize \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& G, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
    const ElementalMatrix<Real>& h, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& s, \
    bool primalInit, bool dualInit, bool standardShift ); \
  template void Initialize \
  ( const SparseMatrix<Real>& JStatic, \
    const Matrix<Real>& regLarge, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
          SparseLDLFactorization<Real>& sparseLDLFact, \
    bool primalInit, bool dualInit, bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl ); \
  template void Initialize \
  ( const DistSparseMatrix<Real>& JStatic, \
    const DistMultiVec<Real>& regLarge, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMultiVec<Real>& s, \
          DistSparseLDLFactorization<Real>& sparseLDLFact, \
    bool primalInit, bool dualInit, bool standardShift, \
    const RegSolveCtrl<Real>& solveCtrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace affine
} // namespace qp
} // namespace El
