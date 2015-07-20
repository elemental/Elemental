/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// x o y = [ x^T y; x0 y1 + y0 x1 ]

template<typename Real>
void SOCApply
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCApply"))
    SOCDots( x, y, z, orders, firstInds );
    auto xRoots = x;
    auto yRoots = y;
    ConeBroadcast( xRoots, orders, firstInds );
    ConeBroadcast( yRoots, orders, firstInds );
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i != firstInds.Get(i,0) )
            z.Update( i, 0, xRoots.Get(i,0)*y.Get(i,0) +
                            yRoots.Get(i,0)*x.Get(i,0) );
}

template<typename Real>
void SOCApply
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Real>& yPre,
        AbstractDistMatrix<Real>& zPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApply"))
    AssertSameGrids( xPre, yPre, zPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto yPtr = ReadProxy<Real,VC,STAR>(&yPre,ctrl);
    auto zPtr = WriteProxy<Real,VC,STAR>(&zPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& y = *yPtr;
    auto& z = *zPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    SOCDots( x, y, z, orders, firstInds );
    auto xRoots = x;
    auto yRoots = y;
    ConeBroadcast( xRoots, orders, firstInds );
    ConeBroadcast( yRoots, orders, firstInds );
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( i != firstInd )
            z.UpdateLocal
            ( iLoc, 0, xRoots.GetLocal(iLoc,0)*y.GetLocal(iLoc,0) +
                       yRoots.GetLocal(iLoc,0)*x.GetLocal(iLoc,0) );
    }
}

template<typename Real>
void SOCApply
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApply"))
    SOCDots( x, y, z, orders, firstInds );
    auto xRoots = x;
    auto yRoots = y;
    ConeBroadcast( xRoots, orders, firstInds );
    ConeBroadcast( yRoots, orders, firstInds );
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int firstInd = firstInds.GetLocal(iLoc,0);
        if( i != firstInd )
            z.UpdateLocal
            ( iLoc, 0, xRoots.GetLocal(iLoc,0)*y.GetLocal(iLoc,0) +
                       yRoots.GetLocal(iLoc,0)*x.GetLocal(iLoc,0) );
    }
}

template<typename Real>
void SOCApply
( const Matrix<Real>& x, 
        Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCApply"))
    // TODO?: Optimize
    Matrix<Real> z;
    SOCApply( x, y, z, orders, firstInds );
    y = z;
}

template<typename Real>
void SOCApply
( const AbstractDistMatrix<Real>& x, 
        AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApply"))
    // TODO?: Optimize
    DistMatrix<Real,VC,STAR> z(x.Grid());
    SOCApply( x, y, z, orders, firstInds, cutoff );
    y = z;
}

template<typename Real>
void SOCApply
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApply"))
    // TODO?: Optimize
    DistMultiVec<Real> z(x.Comm());
    SOCApply( x, y, z, orders, firstInds, cutoff );
    y = z;
}

#define PROTO(Real) \
  template void SOCApply \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& y, \
          Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCApply \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCApply \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff ); \
  template void SOCApply \
  ( const Matrix<Real>& x, \
          Matrix<Real>& y, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCApply \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCApply \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
