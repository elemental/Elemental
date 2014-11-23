/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace lin_prog {

template<typename Real>
void FormSystem
( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l,
  Real tau, Matrix<Real>& J, Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormSystem"))
    const Int m = A.Height();
    const Int n = A.Width();

    // Form the Jacobian, J
    // ====================
    Zeros( J, 2*n+m, 2*n+m );
    IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);
    auto Jss = J(sInd,sInd); auto Jsx = J(sInd,xInd); auto Jsl = J(sInd,lInd);
    auto Jxs = J(xInd,sInd); auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jls = J(lInd,sInd); auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    Diagonal( Jss, x );
    Diagonal( Jsx, s );
    Identity( Jxs, n, n );
    Transpose( A, Jxl ); 
    Jlx = A;

    // Form the right-hand side, y
    // ===========================
    Zeros( y, 2*n+m, 1 );
    auto ys = y(sInd,IR(0,1));
    auto yx = y(xInd,IR(0,1));
    auto yl = y(lInd,IR(0,1));
    for( Int i=0; i<n; ++i )
        ys.Set( i, 0, -x.Get(i,0)*s.Get(i,0)+tau );
    yx = c;
    Gemv( TRANSPOSE, Real(-1), A, l, Real(1), yx );
    Axpy( Real(-1), s, yx );
    yl = b;
    Gemv( NORMAL, Real(-1), A, x, Real(1), yl );
}

template<typename Real>
void FormSystem
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& sPre, const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Real>& l,
  Real tau, AbstractDistMatrix<Real>& JPre, AbstractDistMatrix<Real>& yPre )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::FormSystem"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto JPtr = WriteProxy<Real,MC,MR>(&JPre); auto& J = *JPtr;
    auto yPtr = WriteProxy<Real,MC,MR>(&yPre); auto& y = *yPtr;

    auto xPtr = ReadProxy<Real,STAR,STAR>(&xPre); auto& x = *xPtr;
    auto sPtr = ReadProxy<Real,STAR,STAR>(&sPre); auto& s = *sPtr;

    // Form the Jacobian, J
    // ====================
    Zeros( J, 2*n+m, 2*n+m );
    IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);
    auto Jss = J(sInd,sInd); auto Jsx = J(sInd,xInd); auto Jsl = J(sInd,lInd);
    auto Jxs = J(xInd,sInd); auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jls = J(lInd,sInd); auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    Diagonal( Jss, x.LockedMatrix() );
    Diagonal( Jsx, s.LockedMatrix() );
    Identity( Jxs, n, n );
    Transpose( A, Jxl ); 
    Jlx = A;

    // Form the right-hand side, y
    // ===========================
    Zeros( y, 2*n+m, 1 );
    auto ys = y(sInd,IR(0,1));
    auto yx = y(xInd,IR(0,1));
    auto yl = y(lInd,IR(0,1));
    for( Int iLoc=0; iLoc<ys.LocalHeight(); ++iLoc )
    {
        const Int i = ys.GlobalRow(iLoc);
        ys.SetLocal( iLoc, 0, -x.Get(i,0)*s.Get(i,0)+tau );
    }
    yx = c;
    Gemv( TRANSPOSE, Real(-1), A, l, Real(1), yx );
    Axpy( Real(-1), s, yx );
    yl = b;
    Gemv( NORMAL, Real(-1), A, x, Real(1), yl );
}

template<typename Real>
void SolveSystem
( Int m, Int n, Matrix<Real>& J, Matrix<Real>& y, 
  Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::SolveSystem"))
    if( J.Height() != 2*n+m || J.Width() != 2*n+m )
        LogicError("Jacobian was the wrong size");
    if( y.Height() != 2*n+m || y.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    GaussianElimination( J, y );
    IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);    
    ds = y(sInd,IR(0,1));
    dx = y(xInd,IR(0,1));
    dl = y(lInd,IR(0,1));
}

template<typename Real>
void SolveSystem
( Int m, Int n, 
  AbstractDistMatrix<Real>& J, 
  AbstractDistMatrix<Real>& yPre, 
  AbstractDistMatrix<Real>& ds, 
  AbstractDistMatrix<Real>& dx, 
  AbstractDistMatrix<Real>& dl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::SolveSystem"))
    auto yPtr = ReadWriteProxy<Real,MC,MR>(&yPre);    
    auto& y = *yPtr;

    if( J.Height() != 2*n+m || J.Width() != 2*n+m )
        LogicError("Jacobian was the wrong size");
    if( y.Height() != 2*n+m || y.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    GaussianElimination( J, y );
    IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);    
    ds = y(sInd,IR(0,1));
    dx = y(xInd,IR(0,1));
    dl = y(lInd,IR(0,1));
}

#define PROTO(Real) \
  template void FormSystem \
  ( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& s, const Matrix<Real>& x, const Matrix<Real>& l, \
    Real tau, Matrix<Real>& J, Matrix<Real>& y ); \
  template void FormSystem \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& l, \
    Real tau, AbstractDistMatrix<Real>& J, AbstractDistMatrix<Real>& y ); \
  template void SolveSystem \
  ( Int m, Int n, Matrix<Real>& J, Matrix<Real>& y, \
    Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl ); \
  template void SolveSystem \
  ( Int m, Int n, AbstractDistMatrix<Real>& J, AbstractDistMatrix<Real>& y, \
    AbstractDistMatrix<Real>& ds, AbstractDistMatrix<Real>& dx, \
    AbstractDistMatrix<Real>& dl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
