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

// The full (regularized) KKT system is of the form
//
//   | gamma^2*I     A^T           G       | | dx |   |     -r_c             |
//   |      A    -delta^2*I        0       | | dy |   |     -r_b             |,
//   |      G         0      -(inv(z) o s) | | dz | = | -r_h + inv(z) o r_mu |
//
// where, in the case of a naive update with barrier parameter tau, has
//
//   r_c  = A^T y + G^T z + c,
//   r_b  = A x - b,
//   r_h  = G x + s - h,
//   r_mu = s o z - tau e

template<typename Real>
void KKT
( const Matrix<Real>& A,
  const Matrix<Real>& G,
        Real gamma,
        Real delta,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 

    ShiftDiagonal( Jxx, gamma*gamma );
    ShiftDiagonal( Jyy, -delta*delta ); 

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
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
        Real gamma,
        Real delta,
  const AbstractDistMatrix<Real>& s,
  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& JPre,
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::affine::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);    auto& J = *JPtr;

    Zeros( J, n+m+k, n+m+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    ShiftDiagonal( Jxx, gamma*gamma );
    ShiftDiagonal( Jyy, -delta*delta );

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
        Real gamma,
        Real delta,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::affine::KKT"))
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::affine::KKT( Q, A, G, gamma, delta, s, z, J, onlyLower );
}

template<typename Real>
void StaticKKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::affine::StaticKKT"))
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::affine::StaticKKT( Q, A, G, gamma, delta, J, onlyLower );
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::affine::KKT"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> Q(comm);
    Q.Resize( n, n );
    qp::affine::KKT( Q, A, G, gamma, delta, s, z, J, onlyLower );
}

template<typename Real>
void StaticKKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        DistSparseMatrix<Real>& J,
  bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::affine::StaticKKT"))
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    DistSparseMatrix<Real> Q(comm);
    Q.Resize( n, n );
    qp::affine::StaticKKT( Q, A, G, gamma, delta, J, onlyLower );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
    const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void StaticKKT \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          SparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
    const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower ); \
  template void StaticKKT \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
          Real gamma, \
          Real delta, \
          DistSparseMatrix<Real>& J, \
    bool onlyLower );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace lp
} // namespace El
