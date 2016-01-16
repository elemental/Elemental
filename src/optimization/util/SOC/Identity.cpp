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
void Identity
(       Matrix<Real>& x, 
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::Identity"))
    const Int height = orders.Height();
    DEBUG_ONLY(
      if( firstInds.Height() != height || 
          firstInds.Width() != 1 || orders.Width() != 1 )
          LogicError("orders and firstInds should vectors of the same height");
    )

    const Int* firstIndBuf = firstInds.LockedBuffer();

    Zeros( x, height, 1 );
    Real* xBuf = x.Buffer();
    for( Int i=0; i<height; ++i )
        if( i == firstIndBuf[i] )
            xBuf[i] = 1;
}

template<typename Real,typename>
void Identity
(       ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre )
{
    DEBUG_ONLY(CSE cse("soc::Identity"))
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixWriteProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    const Int height = orders.Height();
    DEBUG_ONLY(
      if( firstInds.Height() != height || 
          firstInds.Width() != 1 || orders.Width() != 1 )
          LogicError("orders and firstInds should vectors of the same height");
    )

    const Int* firstIndBuf = firstInds.LockedBuffer();

    Zeros( x, height, 1 );
    Real* xBuf = x.Buffer();
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstIndBuf[iLoc] )
            xBuf[iLoc] = 1;
    }
}

template<typename Real,typename>
void Identity
(       DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::Identity"))

    const Int height = orders.Height();
    DEBUG_ONLY(
      if( firstInds.Height() != height ||
          firstInds.Width() != 1 || orders.Width() != 1 )
          LogicError("orders and firstInds should vectors of the same height");
    )

    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    x.SetComm( orders.Comm() );
    Zeros( x, height, 1 );
    Real* xBuf = x.Matrix().Buffer();
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstIndBuf[iLoc] )
            xBuf[iLoc] = 1;
    }
}

#define PROTO(Real) \
  template void Identity \
  (       Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Identity \
  (       ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds ); \
  template void Identity \
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
