/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,       Matrix<Real>& s,
  bool primalInitialized, bool dualInitialized,
  bool standardShift )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::Initialize"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    if( primalInitialized )
    {
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( s.Height() != k || s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInitialized )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != k || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInitialized && dualInitialized )
    {
        // TODO: Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    Matrix<Real> J, ones;
    Ones( ones, k, 1 );
    KKT( A, G, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    Matrix<Real> dSub;
    Matrix<Int> p;
    LDL( J, dSub, p, false );

    Matrix<Real> rc, rb, rh, rmu, u, d;
    Zeros( rmu, k, 1 );
    if( !primalInitialized )
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
        Scale( Real(-1), rb );
        rh = h;
        Scale( Real(-1), rh );
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, x, u, s );
        Scale( Real(-1), s );
    }
    if( !dualInitialized )
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
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, u, y, z );
    }

    // alpha_p := min { alpha : s + alpha*e >= 0 }
    // ===========================================
    const auto sMinPair = VectorMin( s );
    const Real alphaPrimal = -sMinPair.value;
    if( alphaPrimal >= Real(0) && primalInitialized )
        RuntimeError("initialized s was non-positive");

    // alpha_d := min { alpha : z + alpha*e >= 0 }
    // ===========================================
    const auto zMinPair = VectorMin( z );
    const Real alphaDual = -zMinPair.value;
    if( alphaDual >= Real(0) && dualInitialized )
        RuntimeError("initialized z was non-positive");

    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real sNorm = Nrm2( s );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
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
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s,
  bool primalInitialized, bool dualInitialized,
  bool standardShift )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::Initialize"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Grid& g = A.Grid();
    if( primalInitialized ) 
    {
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( s.Height() != k || s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInitialized )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != k || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInitialized && dualInitialized )
    {
        // TODO: Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    DistMatrix<Real> J(g), ones(g);
    Ones( ones, k, 1 );
    KKT( A, G, ones, ones, J );

    // Factor the KKT matrix
    // =====================
    DistMatrix<Real> dSub(g);
    DistMatrix<Int> p(g);
    LDL( J, dSub, p, false );

    DistMatrix<Real> rc(g), rb(g), rh(g), rmu(g), d(g), u(g);
    Zeros( rmu, k, 1 );
    if( !primalInitialized )
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
        Scale( Real(-1), rb );
        rh = h;
        Scale( Real(-1), rh );
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, x, u, s );
        Scale( Real(-1), s );
    }
    if( !dualInitialized )
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
        KKTRHS( rc, rb, rh, rmu, ones, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandCoreSolution( m, n, k, d, u, y, z );
    }

    // alpha_p := min { alpha : s + alpha*e >= 0 }
    // ===========================================
    const auto sMinPair = VectorMin( s );
    const Real alphaPrimal = -sMinPair.value;
    if( alphaPrimal >= Real(0) && primalInitialized )
        RuntimeError("initialized s was non-positive");

    // alpha_d := min { alpha : z + alpha*e >= 0 }
    // ===========================================
    const auto zMinPair = VectorMin( z );
    const Real alphaDual = -zMinPair.value;
    if( alphaDual >= Real(0) && dualInitialized )
        RuntimeError("initialized z was non-positive");

    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real sNorm = Nrm2( s );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(epsilon)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(epsilon)*Max(zNorm,Real(1));
    if( standardShift )
    {
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
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,             Matrix<Real>& s,
        vector<Int>& map,            vector<Int>& invMap, 
        Separator& rootSep,          SymmNodeInfo& info,
  bool primalInitialized, bool dualInitialized,
  bool standardShift,     bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::Initialize"))
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::affine::Initialize
    ( Q, A, G, b, c, h, x, y, z, s, map, invMap, rootSep, info,
      primalInitialized, dualInitialized, standardShift, progress );
}

template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A,  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,      const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,            DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,            DistMultiVec<Real>& s,
        DistMap& map,                     DistMap& invMap, 
        DistSeparator& rootSep,           DistSymmNodeInfo& info,
  bool primalInitialized, bool dualInitialized,
  bool standardShift,     bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::Initialize"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> Q(comm);
    Q.Resize( n, n );
    qp::affine::Initialize
    ( Q, A, G, b, c, h, x, y, z, s, 
      map, invMap, rootSep, info, 
      primalInitialized, dualInitialized, standardShift, progress );
}

#define PROTO(Real) \
  template void Initialize \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z,       Matrix<Real>& s, \
    bool primalInitialized, bool dualInitialized, \
    bool standardShift ); \
  template void Initialize \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s, \
    bool primalInitialized, bool dualInitialized, \
    bool standardShift ); \
  template void Initialize \
  ( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z,             Matrix<Real>& s, \
          vector<Int>& map,            vector<Int>& invMap, \
          Separator& rootSep,          SymmNodeInfo& info, \
    bool primalInitialized, bool dualInitialized, \
    bool standardShift,     bool progress ); \
  template void Initialize \
  ( const DistSparseMatrix<Real>& A,  const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b,      const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x,            DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z,            DistMultiVec<Real>& s, \
          DistMap& map,                     DistMap& invMap, \
          DistSeparator& rootSep,           DistSymmNodeInfo& info, \
    bool primalInitialized, bool dualInitialized, \
    bool standardShift,     bool progress );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace lp
} // namespace El
