/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void SOCLowerNorms
( const Matrix<Real>& x, 
        Matrix<Real>& lowerNorms,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCLowerNorms"))
    const Int height = x.Height();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    auto xLower = x;
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) )
            xLower.Set( i, 0, Real(0) );

    SOCDots( xLower, xLower, lowerNorms, orders, firstInds );
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) )
            lowerNorms.Set( i, 0, Sqrt(lowerNorms.Get(i,0)) );
}

template<typename Real>
void SOCLowerNorms
( const AbstractDistMatrix<Real>& xPre, 
        AbstractDistMatrix<Real>& lowerNormsPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCLowerNorms"))
    AssertSameGrids( xPre, lowerNormsPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr          = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto lowerNormsPtr = WriteProxy<Real,VC,STAR>(&lowerNormsPre,ctrl);
    auto ordersPtr     = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr  = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& lowerNorms = *lowerNormsPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    auto xLower = x;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( xLower.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            xLower.SetLocal( iLoc, 0, Real(0) );

    SOCDots( xLower, xLower, lowerNorms, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( lowerNorms.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )    
            lowerNorms.SetLocal( iLoc, 0, Sqrt(lowerNorms.GetLocal(iLoc,0)) );
}

template<typename Real>
void SOCLowerNorms
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& lowerNorms,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCLowerNorms"))
    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    auto xLower = x;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( xLower.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            xLower.SetLocal( iLoc, 0, Real(0) );

    SOCDots( xLower, xLower, lowerNorms, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( lowerNorms.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            lowerNorms.SetLocal( iLoc, 0, Sqrt(lowerNorms.GetLocal(iLoc,0)) );
}

#define PROTO(Real) \
  template void SOCLowerNorms \
  ( const Matrix<Real>& x, \
          Matrix<Real>& lowerNorms, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCLowerNorms \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& lowerNorms, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCLowerNorms \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& lowerNorms, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
