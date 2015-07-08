/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// inv(x) = (R x) / det(x)

template<typename Real>
void SOCInverse
( const Matrix<Real>& x, 
        Matrix<Real>& xInv,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCInverse"))

    Matrix<Real> dInv;
    SOCDets( x, dInv, orders, firstInds );
    ConeBroadcast( dInv, orders, firstInds );
    auto entryInv = [=]( Real alpha ) { return Real(1)/alpha; };
    EntrywiseMap( dInv, function<Real(Real)>(entryInv) );

    auto Rx = x;
    SOCReflect( Rx, orders, firstInds );

    Hadamard( dInv, Rx, xInv ); 
}

template<typename Real>
void SOCInverse
( const AbstractDistMatrix<Real>& xPre, 
        AbstractDistMatrix<Real>& xInvPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCInverse"))
    AssertSameGrids( xPre, xInvPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto xInvPtr = WriteProxy<Real,VC,STAR>(&xInvPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& xInv = *xInvPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> dInv(x.Grid()); 
    SOCDets( x, dInv, orders, firstInds, cutoff );
    ConeBroadcast( dInv, orders, firstInds );
    auto entryInv = [=]( Real alpha ) { return Real(1)/alpha; };
    EntrywiseMap( dInv, function<Real(Real)>(entryInv) );

    auto Rx = x;
    SOCReflect( Rx, orders, firstInds );

    Hadamard( dInv, Rx, xInv );
}

template<typename Real>
void SOCInverse
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& xInv,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCInverse"))

    DistMultiVec<Real> dInv(x.Comm());
    SOCDets( x, dInv, orders, firstInds, cutoff );
    ConeBroadcast( dInv, orders, firstInds );
    auto entryInv = [=]( Real alpha ) { return Real(1)/alpha; };
    EntrywiseMap( dInv, function<Real(Real)>(entryInv) );

    auto Rx = x;
    SOCReflect( Rx, orders, firstInds );

    Hadamard( dInv, Rx, xInv );
}

#define PROTO(Real) \
  template void SOCInverse \
  ( const Matrix<Real>& x, \
          Matrix<Real>& xInv, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCInverse \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& xInv, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCInverse \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& xInv, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
