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
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::AugmentedKKT"))
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
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& x,    const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::AugmentedKKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);
    auto& J = *JPtr;

    Zeros( J, m+n, m+n );
    const IR xInd(0,n), yInd(n,n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd);
    DistMatrix<Real,STAR,STAR> d( z );
    DiagonalSolve( LEFT, NORMAL, x, d );
    Diagonal( Jxx, d.Matrix() );
    Jyx = A;
    if( !onlyLower )
        Transpose( A, Jxy );
}

template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& x,       const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::AugmentedKKT"))
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::direct::AugmentedKKT( Q, A, x, z, J, onlyLower );
}

template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::AugmentedKKT"))
    const Int n = A.Width();
    DistSparseMatrix<Real> Q(A.Comm());
    Q.Resize( n, n );
    qp::direct::AugmentedKKT( Q, A, x, z, J, onlyLower );
}

#define PROTO(Real) \
  template void AugmentedKKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
    Matrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z, \
    AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& x, const Matrix<Real>& z, \
    SparseMatrix<Real>& J, bool onlyLower ); \
  template void AugmentedKKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& x, const DistMultiVec<Real>& z, \
    DistSparseMatrix<Real>& J, bool onlyLower );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace lp
} // namespace El
