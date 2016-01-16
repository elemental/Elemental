/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace soc {

template<typename Real,typename>
void Reflect
(       Matrix<Real>& x, 
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::Reflect"))

    const Int height = x.Height();

    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    const Int* firstIndBuf = firstInds.LockedBuffer();
    Real* xBuf = x.Buffer(); 
    
    for( Int i=0; i<height; ++i )
        if( i != firstIndBuf[i] )
            xBuf[i] = -xBuf[i];
}

template<typename Real,typename>
void Reflect
(       ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre )
{
    DEBUG_ONLY(CSE cse("soc::Reflect"))
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadWriteProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl ); 
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    DEBUG_ONLY(
      const Int height = x.Height();
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    const Int localHeight = x.LocalHeight();
    Real* xBuf = x.Buffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( x.GlobalRow(iLoc) != firstIndBuf[iLoc] )
            xBuf[iLoc] = -xBuf[iLoc];
}

template<typename Real,typename>
void Reflect
(       DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::Reflect"))

    DEBUG_ONLY(
      const Int height = x.Height();
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    const Int firstLocalRow = x.FirstLocalRow();
    const Int localHeight = x.LocalHeight();
    Real* xBuf = x.Matrix().Buffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( iLoc+firstLocalRow != firstIndBuf[iLoc] )
            xBuf[iLoc] = -xBuf[iLoc];
}

#define PROTO(Real) \
  template void Reflect \
  (       Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Reflect \
  (       ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds ); \
  template void Reflect \
  (       DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace soc
} // namespace El
