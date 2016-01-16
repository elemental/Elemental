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
void Shift
(       Matrix<Real>& x,
        Real shift,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::Shift"))
    const Int height = x.Height();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    Real* xBuf = x.Buffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();
    for( Int i=0; i<height; ++i )
        if( i == firstIndBuf[i] )
            xBuf[i] += shift;
}

template<typename Real,typename>
void Shift
(       ElementalMatrix<Real>& xPre,
        Real shift,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre )
{
    DEBUG_ONLY(CSE cse("soc::Shift"))
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

    const Int height = x.Height();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    Real* xBuf = x.Buffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( x.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            xBuf[iLoc] += shift;
}

template<typename Real,typename>
void Shift
(       DistMultiVec<Real>& x,
        Real shift,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::Shift"))

    const Int height = x.Height();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 ) 
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();
    Real* xBuf = x.Matrix().Buffer();

    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( x.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            xBuf[iLoc] += shift;
}

#define PROTO(Real) \
  template void Shift \
  (       Matrix<Real>& x, Real shift, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Shift \
  (       ElementalMatrix<Real>& x, Real shift, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds ); \
  template void Shift \
  (       DistMultiVec<Real>& x, Real shift, \
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
