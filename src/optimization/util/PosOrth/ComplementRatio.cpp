/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace pos_orth {

// Compute Max( s o z ) / Min( s o z ) to determine if we need to recenter

template<typename Real,typename>
Real ComplementRatio
( const Matrix<Real>& s,
  const Matrix<Real>& z )
{
    DEBUG_CSE
    const Int k = s.Height();
    const Real* sBuf = s.LockedBuffer();
    const Real* zBuf = z.LockedBuffer();

    Real maxProd = 0;
    for( Int i=0; i<k; ++i )
        maxProd = Max( sBuf[i]*zBuf[i], maxProd );

    Real minProd = maxProd;
    for( Int i=0; i<k; ++i )
        minProd = Min( sBuf[i]*zBuf[i], minProd );

    return maxProd/minProd;
}

template<typename Real,typename>
Real ComplementRatio
( const ElementalMatrix<Real>& sPre,
  const ElementalMatrix<Real>& zPre )
{
    DEBUG_CSE

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      sProx( sPre, ctrl ),
      zProx( zPre, ctrl );
    auto& s = sProx.GetLocked();
    auto& z = zProx.GetLocked();

    const Int localHeight = s.LocalHeight();
    const Real* sBuf = s.LockedBuffer();
    const Real* zBuf = z.LockedBuffer();

    Real maxLocProd = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        maxLocProd = Max( sBuf[iLoc]*zBuf[iLoc], maxLocProd );
    const Real maxProd = mpi::AllReduce( maxLocProd, mpi::MAX, s.DistComm() );
    
    Real minLocProd = maxProd;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        minLocProd = Min( sBuf[iLoc]*zBuf[iLoc], minLocProd );
    const Real minProd = mpi::AllReduce( minLocProd, mpi::MIN, s.DistComm() );

    return maxProd/minProd;
}

template<typename Real,typename>
Real ComplementRatio
( const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z )
{
    DEBUG_CSE
    const Int localHeight = s.LocalHeight();
    const Real* sBuf = s.LockedMatrix().LockedBuffer();
    const Real* zBuf = z.LockedMatrix().LockedBuffer();

    Real maxLocProd = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        maxLocProd = Max( sBuf[iLoc]*zBuf[iLoc], maxLocProd );
    const Real maxProd = mpi::AllReduce( maxLocProd, mpi::MAX, s.Comm() );
    
    Real minLocProd = maxProd;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        minLocProd = Min( sBuf[iLoc]*zBuf[iLoc], minLocProd );
    const Real minProd = mpi::AllReduce( minLocProd, mpi::MIN, s.Comm() );

    return maxProd/minProd;
}

#define PROTO(Real) \
  template Real ComplementRatio \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& z ); \
  template Real ComplementRatio \
  ( const ElementalMatrix<Real>& s, \
    const ElementalMatrix<Real>& z ); \
  template Real ComplementRatio \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace pos_orth
} // namespace El
