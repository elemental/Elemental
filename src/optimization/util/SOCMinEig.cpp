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
void SOCMinEig
( const Matrix<Real>& x, 
        Matrix<Real>& minEigs,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCMinEig"))
    SOCLowerNorms( x, minEigs, orders, firstInds );
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            minEigs.Set( i, 0, x.Get(i,0)-minEigs.Get(i,0) );
}

template<typename Real>
void SOCMinEig
( const AbstractDistMatrix<Real>& xPre, 
        AbstractDistMatrix<Real>& minEigsPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMinEig"))
    AssertSameGrids( xPre, minEigsPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr         = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto minEigsPtr   = WriteProxy<Real,VC,STAR>(&minEigsPre,ctrl);
    auto ordersPtr    = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& minEigs = *minEigsPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    SOCLowerNorms( x, minEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)-minEigs.GetLocal(iLoc,0) );
}

template<typename Real>
void SOCMinEig
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& minEigs,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMinEig"))
    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    SOCLowerNorms( x, minEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)-minEigs.GetLocal(iLoc,0) );
}

template<typename Real>
Real SOCMinEig
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCMinEig"))
    Matrix<Real> minEigs;
    SOCMinEig( x, minEigs, orders, firstInds );

    Real minEig = std::numeric_limits<Real>::max();
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            minEig = Min(minEigs.Get(i,0),minEig);
    return minEig;
}

template<typename Real>
Real SOCMinEig
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMinEig"))
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr         = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto ordersPtr    = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> minEigs(x.Grid());
    SOCMinEig( x, minEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real minEigLocal = std::numeric_limits<Real>::max();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigLocal = Min(minEigLocal,minEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( minEigLocal, mpi::MIN, x.DistComm() );
}

template<typename Real>
Real SOCMinEig
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMinEig"))
    DistMultiVec<Real> minEigs(x.Comm());
    SOCMinEig( x, minEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real minEigLocal = std::numeric_limits<Real>::max();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigLocal = Min(minEigLocal,minEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( minEigLocal, mpi::MIN, x.Comm() );
}

#define PROTO(Real) \
  template void SOCMinEig \
  ( const Matrix<Real>& x, \
          Matrix<Real>& minEigs, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCMinEig \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& minEigs, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCMinEig \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& minEigs, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff ); \
  template Real SOCMinEig \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Real SOCMinEig \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template Real SOCMinEig \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
