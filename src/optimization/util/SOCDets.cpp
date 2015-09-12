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
void SOCDets
( const Matrix<Real>& x, 
        Matrix<Real>& d,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCDets"))
    auto Rx = x;
    SOCReflect( Rx, orders, firstInds );
    SOCDots( x, Rx, d, orders, firstInds );
}

template<typename Real>
void SOCDets
( const ElementalMatrix<Real>& xPre, 
        ElementalMatrix<Real>& dPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCDets"))
    AssertSameGrids( xPre, dPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto dPtr = WriteProxy<Real,VC,STAR>(&dPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& d = *dPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    auto Rx = x;
    SOCReflect( Rx, orders, firstInds );
    SOCDots( x, Rx, d, orders, firstInds, cutoff );
}

template<typename Real>
void SOCDets
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& d,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCDets"))
    auto Rx = x;
    SOCReflect( Rx, orders, firstInds );
    SOCDots( x, Rx, d, orders, firstInds, cutoff );
}

#define PROTO(Real) \
  template void SOCDets \
  ( const Matrix<Real>& x, \
          Matrix<Real>& d, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCDets \
  ( const ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& d, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCDets \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& d, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
