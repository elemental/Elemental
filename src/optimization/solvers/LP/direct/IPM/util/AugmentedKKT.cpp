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

// Form 
//
//    | (x <> z)  A^T | | dx | = | -r_c - x <> r_mu |,
//    |    A       0  | | dy |   | -r_b             |
//
// where 
//
//    r_b  = A x - b,
//    r_c  = A^T y - z + c,
//    r_mu = x o z - tau e,
//
// and dz can be computed using
//
//   dz = - x <> (r_mu + z o dx)
//

template<typename Real>
void AugmentedKKT
( const Matrix<Real>& A, 
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), yInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd);
    Matrix<Real> d( z );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d );
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );
}

template<typename Real>
void AugmentedKKT
( const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& JPre, 
  bool onlyLower )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();

    DistMatrixWriteProxy<Real,Real,MC,MR> JProx( JPre );
    auto& J = JProx.Get();

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), yInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd);
    DistMatrix<Real,MC,STAR> d( z );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d );
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );
}

template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& A, 
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    DEBUG_CSE
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Zeros( Q, n, n );
    qp::direct::AugmentedKKT( Q, A, gamma, delta, x, z, J, onlyLower );
}

template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower )
{
    DEBUG_CSE
    const Int n = A.Width();
    DistSparseMatrix<Real> Q(A.Comm());
    Zeros( Q, n, n );
    qp::direct::AugmentedKKT( Q, A, gamma, delta, x, z, J, onlyLower );
}

#define PROTO(Real) \
  template void AugmentedKKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          Matrix<Real>& J, \
    bool onlyLower ); \
  template void AugmentedKKT \
  ( const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& J, \
    bool onlyLower ); \
  template void AugmentedKKT \
  ( const SparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void AugmentedKKT \
  ( const DistSparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower );

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
