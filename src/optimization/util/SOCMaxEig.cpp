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
void SOCMaxEig
( const Matrix<Real>& x, 
        Matrix<Real>& maxEigs,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCMaxEig"))
    SOCLowerNorms( x, maxEigs, orders, firstInds );
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            maxEigs.Set( i, 0, x.Get(i,0)+maxEigs.Get(i,0) );
}

template<typename Real>
void SOCMaxEig
( const AbstractDistMatrix<Real>& xPre, 
        AbstractDistMatrix<Real>& maxEigsPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMaxEig"))
    AssertSameGrids( xPre, maxEigsPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr         = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto maxEigsPtr   = WriteProxy<Real,VC,STAR>(&maxEigsPre,ctrl);
    auto ordersPtr    = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& maxEigs = *maxEigsPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    SOCLowerNorms( x, maxEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)+maxEigs.GetLocal(iLoc,0) );
}

template<typename Real>
void SOCMaxEig
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& maxEigs,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMaxEig"))
    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    SOCLowerNorms( x, maxEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)+maxEigs.GetLocal(iLoc,0) );
}

template<typename Real>
Real SOCMaxEig
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCMaxEig"))
    Matrix<Real> maxEigs;
    SOCMaxEig( x, maxEigs, orders, firstInds );

    Real maxEig = std::numeric_limits<Real>::min();
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            maxEig = Max(maxEigs.Get(i,0),maxEig);
    return maxEig;
}

template<typename Real>
Real SOCMaxEig
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMaxEig"))
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

    DistMatrix<Real,VC,STAR> maxEigs(x.Grid());
    SOCMaxEig( x, maxEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real maxEigLocal = std::numeric_limits<Real>::min();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigLocal = Max(maxEigLocal,maxEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( maxEigLocal, mpi::MAX, x.DistComm() );
}

template<typename Real>
Real SOCMaxEig
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCMaxEig"))
    DistMultiVec<Real> maxEigs(x.Comm());
    SOCMaxEig( x, maxEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real maxEigLocal = std::numeric_limits<Real>::min();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigLocal = Max(maxEigLocal,maxEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( maxEigLocal, mpi::MAX, x.Comm() );
}

#define PROTO(Real) \
  template void SOCMaxEig \
  ( const Matrix<Real>& x, \
          Matrix<Real>& maxEigs, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCMaxEig \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& maxEigs, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCMaxEig \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& maxEigs, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff ); \
  template Real SOCMaxEig \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Real SOCMaxEig \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template Real SOCMaxEig \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
