/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace lp {
namespace dual {

//     | 0 A^T   G^T  |
// J = | A  0     0   |, with the variable ordering (x,y,z)
//     | G  0  -W^T W |

// TODO: Add support for non-Nesterov/Todd other scalings
template<typename Real>
void KKT
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& w,
  Matrix<Real>& J, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    // t := - diag(l)^T diag(l)
    Matrix<Real> t;
    t = w;
    DiagonalScale( LEFT, NORMAL, w, t );
    Scale( Real(-1), t );

    Zeros( J, m+n+k, m+n+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m+k,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 
    Jyx = A;
    Jzx = G;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        Transpose( A, Jxy );
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& w, 
  AbstractDistMatrix<Real>& JPre, bool onlyLower )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre);    auto& J = *JPtr;

    // t := - diag(w)^T diag(w)
    DistMatrix<Real,MC,STAR> t(A.Grid());
    t = w;
    DiagonalScale( LEFT, NORMAL, w, t );
    Scale( Real(-1), t );

    Zeros( J, m+n+k, m+n+k );
    const IR xInd(0,n), yInd(n,n+m), zInd(n+m+k,n+m+k);
    auto Jxx = J(xInd,xInd); auto Jxy = J(xInd,yInd); auto Jxz = J(xInd,zInd); 
    auto Jyx = J(yInd,xInd); auto Jyy = J(yInd,yInd); auto Jyz = J(yInd,zInd); 
    auto Jzx = J(zInd,xInd); auto Jzy = J(zInd,yInd); auto Jzz = J(zInd,zInd); 
    Jyx = A;
    Jzx = G;
    Diagonal( Jzz, t );

    if( !onlyLower )
    {
        Transpose( A, Jxy );
        Transpose( G, Jxz );
    }
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc, const Matrix<Real>& rb,
  const Matrix<Real>& rh, const Matrix<Real>& rmu,
  const Matrix<Real>& w, const Matrix<Real>& l,
  Matrix<Real>& d )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::KKTRHS"))
    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Zeros( d, n+m+k, 1 );

    auto dx = d(xInd,IR(0,1));
    dx = rc;
    Scale( Real(-1), dx );

    auto dy = d(yInd,IR(0,1));
    dy = rb;
    Scale( Real(-1), dy );

    // dz = -r_h + W^T (l <> r_mu)
    // ===========================
    auto dz = d(zInd,IR(0,1));
    dz = rmu;
    DiagonalSolve( LEFT, NORMAL, l, dz );
    DiagonalScale( LEFT, NORMAL, w, dz );
    Axpy( Real(-1), rh, dz );
}

template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rb, 
  const AbstractDistMatrix<Real>& rh, const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& w,  const AbstractDistMatrix<Real>& l,
        AbstractDistMatrix<Real>& dPre )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::KKTRHS"))
    const Int m = rb.Height();
    const Int n = rc.Height();
    const Int k = rh.Height();

    auto dPtr = WriteProxy<Real,MC,MR>(&dPre); 
    auto& d = *dPtr;

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Zeros( d, n+m+k, 1 );

    auto dx = d(xInd,IR(0,1));
    Copy( rc, dx );
    Scale( Real(-1), dx );

    auto dy = d(yInd,IR(0,1));
    Copy( rb, dy );
    Scale( Real(-1), dy );

    // dz = -r_h + W^T (l <> r_mu)
    // ===========================
    auto dz = d(zInd,IR(0,1));
    Copy( rmu, dz );
    DiagonalSolve( LEFT, NORMAL, l, dz );
    DiagonalScale( LEFT, NORMAL, w, dz );
    Axpy( Real(-1), rh, dz );
}

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, Int k, const Matrix<Real>& d, 
  const Matrix<Real>& rmu,
  const Matrix<Real>& w, const Matrix<Real>& l,
  Matrix<Real>& dx, Matrix<Real>& dy, 
  Matrix<Real>& dz, Matrix<Real>& ds )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::ExpandKKTSolution"))
    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    dx = d(xInd,IR(0,1));
    dy = d(yInd,IR(0,1));
    dz = d(zInd,IR(0,1));

    // ds = -W^T (l <> r_mu + W dz)
    // ============================
    ds = dz;
    DiagonalScale( LEFT, NORMAL, w, ds );

    Matrix<Real> t;
    t = rmu;
    DiagonalSolve( LEFT, NORMAL, l, t );
    Axpy( Real(1), t, ds );

    DiagonalScale( LEFT, NORMAL, w, ds );
    Scale( Real(-1), ds );
}

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, Int k, const AbstractDistMatrix<Real>& dPre, 
  const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& w,   const AbstractDistMatrix<Real>& l,
  AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
  AbstractDistMatrix<Real>& dz, AbstractDistMatrix<Real>& ds )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::ExpandKKTSolution"))
    
    auto dPtr = ReadProxy<Real,MC,MR>(&dPre);    
    auto& d = *dPtr;

    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR xInd(0,n), yInd(n,n+m), zInd(n+m,n+m+k);
    Copy( d(xInd,IR(0,1)), dx );
    Copy( d(yInd,IR(0,1)), dy );
    Copy( d(zInd,IR(0,1)), dz );

    // ds = -W^T (l <> r_mu + W dz)
    // ============================
    Copy( dz, ds );
    DiagonalScale( LEFT, NORMAL, w, ds );

    DistMatrix<Real> t(d.Grid());
    t = rmu;
    DiagonalSolve( LEFT, NORMAL, l, t );
    Axpy( Real(1), t, ds );

    DiagonalScale( LEFT, NORMAL, w, ds );
    Scale( Real(-1), ds );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& w, \
    Matrix<Real>& J, bool onlyLower ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& w, \
    AbstractDistMatrix<Real>& J, bool onlyLower ); \
  template void KKTRHS \
  ( const Matrix<Real>& rc, const Matrix<Real>& rb, \
    const Matrix<Real>& rh, const Matrix<Real>& rmu, \
    const Matrix<Real>& w,  const Matrix<Real>& l, \
          Matrix<Real>& d ); \
  template void KKTRHS \
  ( const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rb, \
    const AbstractDistMatrix<Real>& rh, const AbstractDistMatrix<Real>& rmu, \
    const AbstractDistMatrix<Real>& w,  const AbstractDistMatrix<Real>& l, \
          AbstractDistMatrix<Real>& d ); \
  template void ExpandKKTSolution \
  ( Int m, Int n, Int k, const Matrix<Real>& d, \
    const Matrix<Real>& rmu, \
    const Matrix<Real>& w, const Matrix<Real>& l, \
    Matrix<Real>& dx, Matrix<Real>& dy, \
    Matrix<Real>& dz, Matrix<Real>& ds ); \
  template void ExpandKKTSolution \
  ( Int m, Int n, Int k, const AbstractDistMatrix<Real>& d, \
    const AbstractDistMatrix<Real>& rmu, \
    const AbstractDistMatrix<Real>& w, const AbstractDistMatrix<Real>& l, \
    AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, \
    AbstractDistMatrix<Real>& dz, AbstractDistMatrix<Real>& ds );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace dual
} // namespace lp
} // namespace El
