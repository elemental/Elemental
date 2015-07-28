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

// The full KKT system is of the form
//
//   |  gamma^2*I     A^T          -I     | | dx |   |   -r_c        |
//   |       A    -delta^2*I        0     | | dy |   |   -r_b        |,
//   |      -I         0      -inv(z) o x | | dz | = | inv(z) o r_mu |
//
// where, in the case of a naive update with barrier parameter tau,
//
//   r_c = A^T y - z + c,
//   r_b = A x - b,
//   r_mu = x o z - tau e.

template<typename Real>
void KKT
( const Matrix<Real>& A, 
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, 2*n+m, 2*n+m );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 

    // Jxx := gamma^2*I
    // ================
    ShiftDiagonal( Jxx, gamma*gamma );

    // Jyy := -delta^2*I
    // =================
    ShiftDiagonal( Jyy, -delta*delta );

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := -I
    // =========
    ShiftDiagonal( Jzx, Real(-1) );

    // Jzz := - z <> x
    // ===============
    Matrix<Real> t;
    t = x;
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy ); 

        // Jxz := -I
        // =========
        ShiftDiagonal( Jxz, Real(-1) );
    }
}

template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A, 
        Real gamma,
        Real delta,
  const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);    auto& J = *JPtr;

    Zeros( J, 2*n+m, 2*n+m );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,2*n+m);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd);
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd);
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd);

    ShiftDiagonal( Jxx, gamma*gamma );
    ShiftDiagonal( Jyy, -delta*delta );

    // Jyx := A
    // ========
    Jyx = A;

    // Jzx := -I
    // =========
    ShiftDiagonal( Jzx, Real(-1) );

    // Jzz := - z <> x
    // ===============
    DistMatrix<Real,MC,STAR> t(x);
    DiagonalSolve( LEFT, NORMAL, z, t );
    t *= -1;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        // Jxy := A^T
        // ==========
        Transpose( A, Jxy );

        // Jxz := -I
        // =========
        ShiftDiagonal( Jxz, Real(-1) );
    }
}

template<typename Real>
void KKT
( const SparseMatrix<Real>& A, 
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::KKT"))
    const Int n = A.Width();
    SparseMatrix<Real> Q;
    Q.Resize( n, n );
    qp::direct::KKT( Q, A, gamma, delta, x, z, J, onlyLower );
}

template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A, 
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CSE cse("lp::direct::KKT"))
    const Int n = A.Width();
    DistSparseMatrix<Real> Q(A.Comm());
    Q.Resize( n, n );
    qp::direct::KKT( Q, A, gamma, delta, x, z, J, onlyLower );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const SparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const Matrix<Real>& x, \
    const Matrix<Real>& z, \
          SparseMatrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const DistSparseMatrix<Real>& A, \
          Real gamma, \
          Real delta, \
    const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& z, \
          DistSparseMatrix<Real>& J, bool onlyLower );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace lp
} // namespace El
