/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace soc {

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void Reflect
(       Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    EL_DEBUG_CSE
    const Int height = x.Height();
    EL_DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )
    for( Int i=0; i<height; ++i )
        if( i != firstInds(i) )
            x(i) = -x(i);
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void Reflect
(       AbstractDistMatrix<Real>& xPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre )
{
    EL_DEBUG_CSE
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

    EL_DEBUG_ONLY(
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

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void Reflect
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds )
{
    EL_DEBUG_CSE

    EL_DEBUG_ONLY(
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
  (       AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds ); \
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
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
