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
void LowerNorms
( const Matrix<Real>& x,
        Matrix<Real>& lowerNorms,
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

    auto xLower = x;
    for( Int i=0; i<height; ++i )
        if( i == firstInds(i) )
            xLower(i) = 0;

    soc::Dots( xLower, xLower, lowerNorms, orders, firstInds );
    for( Int i=0; i<height; ++i )
        if( i == firstInds(i) )
            lowerNorms(i) = Sqrt(lowerNorms(i));
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void LowerNorms
( const AbstractDistMatrix<Real>& xPre,
        AbstractDistMatrix<Real>& lowerNormsPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    EL_DEBUG_CSE
    AssertSameGrids( xPre, lowerNormsPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixWriteProxy<Real,Real,VC,STAR>
      lowerNormsProx( lowerNormsPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& lowerNorms = lowerNormsProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    const Int localHeight = x.LocalHeight();
    EL_DEBUG_ONLY(
      const Int height = x.Height();
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    const Int* firstIndBuf = firstInds.LockedBuffer();

    auto xLower = x;
    Real* xLowerBuf = xLower.Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( xLower.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            xLowerBuf[iLoc] = 0;

    soc::Dots( xLower, xLower, lowerNorms, orders, firstInds, cutoff );
    Real* lowerNormBuf = lowerNorms.Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( lowerNorms.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            lowerNormBuf[iLoc] = Sqrt(lowerNormBuf[iLoc]);
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void LowerNorms
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& lowerNorms,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    EL_DEBUG_CSE
    const Int localHeight = x.LocalHeight();
    EL_DEBUG_ONLY(
      const Int height = x.Height();
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    auto xLower = x;
    Real* xLowerBuf = xLower.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( xLower.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            xLowerBuf[iLoc] = 0;

    soc::Dots( xLower, xLower, lowerNorms, orders, firstInds, cutoff );
    Real* lowerNormBuf = lowerNorms.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( lowerNorms.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            lowerNormBuf[iLoc] = Sqrt(lowerNormBuf[iLoc]);
}

#define PROTO(Real) \
  template void LowerNorms \
  ( const Matrix<Real>& x, \
          Matrix<Real>& lowerNorms, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void LowerNorms \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& lowerNorms, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void LowerNorms \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& lowerNorms, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
