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

// The full KKT system is of the form
//
//   | 0 A^T      G     | | x |   |        -c             |
//   | A 0        0     | | y |   |         b             |,
//   | G 0    -(z <> s) | | z | = | -z <> (s o z + tau e) |
//
// and the particular system solved is of the form
//
//   | 0 A^T      G     | | dx |   |     -r_c        |
//   | A 0        0     | | dy |   |     -r_b        |,
//   | G 0    -(z <> s) | | dz | = | -rh + z <> r_mu |
//
// where
//
//   r_c  = A^T y + G^T z + c,
//   r_b  = A x - b,
//   r_h  = G x + s - h,
//   r_mu = s o z - tau e

template<typename Real>
void KKT
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := G
    // ========
    Jzx = G;

    // Jzz := - z <> s
    // ===============
    Matrix<Real> t;
    t = s;
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := G^T
        // ==========
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKT
( const DistMatrix<Real>& A,
  const DistMatrix<Real>& G,
  const DistMatrix<Real>& s,
  const DistMatrix<Real>& z,
        DistMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := G
    // ========
    Jzx = G;

    // Jzz := - z <> s
    // ===============
    DistMatrix<Real,MC,STAR> t(s);
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := G
        // ========
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::affine::KKT( Q, A, G, s, z, J, onlyLower );
}

template<typename Real>
void StaticKKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::affine::StaticKKT( Q, A, G, gamma, delta, beta, J, onlyLower );
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int n = A.Width();
    DistSparseMatrix<Real> Q(A.Grid());
    Q.Resize( n, n );
    qp::affine::KKT( Q, A, G, s, z, J, onlyLower );
}

template<typename Real>
void StaticKKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        DistSparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int n = A.Width();
    DistSparseMatrix<Real> Q(A.Grid());
    Q.Resize( n, n );
    qp::affine::StaticKKT( Q, A, G, gamma, delta, beta, J, onlyLower );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const DistMatrix<Real>& A, \
    const DistMatrix<Real>& G, \
    const DistMatrix<Real>& s, \
    const DistMatrix<Real>& z, \
          DistMatrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void StaticKKT \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          Real beta, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void StaticKKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          Real beta, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower );

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
