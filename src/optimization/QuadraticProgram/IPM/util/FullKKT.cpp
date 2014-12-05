/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace quad_prog {

//     | X  S 0   |
// J = | I -Q A^T |, with the variable ordering (s,x,l)
//     | 0  A 0   |

template<typename Real>
void KKT
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& s, const Matrix<Real>& x,
        Matrix<Real>& J )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    Zeros( J, 2*n+m, 2*n+m );
    const IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);
    auto Jss = J(sInd,sInd); auto Jsx = J(sInd,xInd); auto Jsl = J(sInd,lInd);
    auto Jxs = J(xInd,sInd); auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jls = J(lInd,sInd); auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    Diagonal( Jss, x );
    Diagonal( Jsx, s );
    Identity( Jxs, n, n );
    Jxx = Q;
    Scale( Real(-1), Jxx );
    Transpose( A, Jxl ); 
    Jlx = A;
}

template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& Q,    const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& sPre, const AbstractDistMatrix<Real>& xPre,
  AbstractDistMatrix<Real>& JPre )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::KKT"))
    const Int m = A.Height();
    const Int n = A.Width();

    auto xPtr = ReadProxy<Real,STAR,STAR>(&xPre); auto& x = *xPtr;
    auto sPtr = ReadProxy<Real,STAR,STAR>(&sPre); auto& s = *sPtr;
    auto JPtr = WriteProxy<Real,MC,MR>(&JPre); auto& J = *JPtr;

    Zeros( J, 2*n+m, 2*n+m );
    IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);
    auto Jss = J(sInd,sInd); auto Jsx = J(sInd,xInd); auto Jsl = J(sInd,lInd);
    auto Jxs = J(xInd,sInd); auto Jxx = J(xInd,xInd); auto Jxl = J(xInd,lInd);
    auto Jls = J(lInd,sInd); auto Jlx = J(lInd,xInd); auto Jll = J(lInd,lInd);
    Diagonal( Jss, x.LockedMatrix() );
    Diagonal( Jsx, s.LockedMatrix() );
    Identity( Jxs, n, n );
    Jxx = Q;
    Scale( Real(-1), Jxx );
    Transpose( A, Jxl );
    Jlx = A;
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
  Matrix<Real>& y )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::KKTRHS"))
    const Int m = rb.Height();
    const Int n = rc.Height();
    const IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);
    Zeros( y, 2*n+m, 1 );

    auto ys = y(sInd,IR(0,1));
    ys = rmu;
    Scale( Real(-1), ys );

    auto yx = y(xInd,IR(0,1));
    yx = rc;
    Scale( Real(-1), yx );

    auto yl = y(lInd,IR(0,1));
    yl = rb;
    Scale( Real(-1), yl );
}

template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, 
  const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& yPre )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::KKTRHS"))

    auto yPtr = WriteProxy<Real,MC,MR>(&yPre); 
    auto& y = *yPtr;

    const Int m = rb.Height();
    const Int n = rc.Height();
    const IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);
    Zeros( y, 2*n+m, 1 );

    auto ys = y(sInd,IR(0,1));
    Copy( rmu, ys );
    Scale( Real(-1), ys );

    auto yx = y(xInd,IR(0,1));
    Copy( rc, yx );
    Scale( Real(-1), yx );

    auto yl = y(lInd,IR(0,1));
    Copy( rb, yl );
    Scale( Real(-1), yl );
}

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const Matrix<Real>& y, 
  Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::ExpandKKTSolution"))
    if( y.Height() != 2*n+m || y.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);    
    ds = y(sInd,IR(0,1));
    dx = y(xInd,IR(0,1));
    dl = y(lInd,IR(0,1));
}

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const AbstractDistMatrix<Real>& yPre, 
  AbstractDistMatrix<Real>& ds, AbstractDistMatrix<Real>& dx, 
  AbstractDistMatrix<Real>& dl )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::ExpandKKTSolution"))
    
    auto yPtr = ReadProxy<Real,MC,MR>(&yPre);    
    auto& y = *yPtr;

    if( y.Height() != 2*n+m || y.Width() != 1 )
        LogicError("Right-hand side was the wrong size");

    const IR sInd(0,n), xInd(n,2*n), lInd(2*n,2*n+m);    
    Copy( y(sInd,IR(0,1)), ds );
    Copy( y(xInd,IR(0,1)), dx );
    Copy( y(lInd,IR(0,1)), dl );
}

#define PROTO(Real) \
  template void KKT \
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& s, const Matrix<Real>& x, \
    Matrix<Real>& J ); \
  template void KKT \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x, \
    AbstractDistMatrix<Real>& J ); \
  template void KKTRHS \
  ( const Matrix<Real>& rmu, const Matrix<Real>& rc, \
    const Matrix<Real>& rb, Matrix<Real>& y ); \
  template void KKTRHS \
  ( const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, \
    const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& y ); \
  template void ExpandKKTSolution \
  ( Int m, Int n, const Matrix<Real>& yPre, \
    Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl ); \
  template void ExpandKKTSolution \
  ( Int m, Int n, const AbstractDistMatrix<Real>& yPre, \
    AbstractDistMatrix<Real>& ds, AbstractDistMatrix<Real>& dx, \
    AbstractDistMatrix<Real>& dl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace quad_prog
} // namespace El
