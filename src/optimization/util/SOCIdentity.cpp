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
void SOCIdentity
(       Matrix<Real>& x, 
  const Matrix<Int>& orders, const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCIdentity"))
    const Int height = orders.Height();
    if( firstInds.Height() != height || 
        firstInds.Width() != 1 || orders.Width() != 1 )
        LogicError("orders and firstInds should vectors of the same height");

    Zeros( x, height, 1 );
    for( Int i=0; i<height; ++i )
    {
        if( i == firstInds.Get(i,0) )
            x.Set( i, 0, Real(1) );
    }
}

template<typename Real>
void SOCIdentity
(       AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre )
{
    DEBUG_ONLY(CSE cse("SOCIdentity"))
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = WriteProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Int height = orders.Height();
    if( firstInds.Height() != height || 
        firstInds.Width() != 1 || orders.Width() != 1 )
        LogicError("orders and firstInds should vectors of the same height");

    Zeros( x, height, 1 );
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstInds.GetLocal(iLoc,0) )
            x.SetLocal( iLoc, 0, Real(1) );
    }
}

template<typename Real>
void SOCIdentity
(       DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCIdentity"))

    const Int height = orders.Height();
    if( firstInds.Height() != height ||
        firstInds.Width() != 1 || orders.Width() != 1 )
        LogicError("orders and firstInds should vectors of the same height");

    x.SetComm( orders.Comm() );
    Zeros( x, height, 1 );
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstInds.GetLocal(iLoc,0) )
            x.SetLocal( iLoc, 0, Real(1) );
    }
}

#define PROTO(Real) \
  template void SOCIdentity \
  (       Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCIdentity \
  (       AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds ); \
  template void SOCIdentity \
  (       DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
