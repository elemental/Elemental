/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace soc {

// TODO: Lower-level access

template<typename Real>
void MaxEig
( const Matrix<Real>& x, 
        Matrix<Real>& maxEigs,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    soc::LowerNorms( x, maxEigs, orders, firstInds );
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            maxEigs.Set( i, 0, x.Get(i,0)+maxEigs.Get(i,0) );
}

template<typename Real>
void MaxEig
( const ElementalMatrix<Real>& xPre, 
        ElementalMatrix<Real>& maxEigsPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
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

    soc::LowerNorms( x, maxEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)+maxEigs.GetLocal(iLoc,0) );
}

template<typename Real>
void MaxEig
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& maxEigs,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    soc::LowerNorms( x, maxEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)+maxEigs.GetLocal(iLoc,0) );
}

template<typename Real>
Real MaxEig
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    Matrix<Real> maxEigs;
    soc::MaxEig( x, maxEigs, orders, firstInds );

    Real maxEig = std::numeric_limits<Real>::min();
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            maxEig = Max(maxEigs.Get(i,0),maxEig);
    return maxEig;
}

template<typename Real>
Real MaxEig
( const ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
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
    soc::MaxEig( x, maxEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real maxEigLocal = std::numeric_limits<Real>::min();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigLocal = Max(maxEigLocal,maxEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( maxEigLocal, mpi::MAX, x.DistComm() );
}

template<typename Real>
Real MaxEig
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    DistMultiVec<Real> maxEigs(x.Comm());
    soc::MaxEig( x, maxEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real maxEigLocal = std::numeric_limits<Real>::min();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            maxEigLocal = Max(maxEigLocal,maxEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( maxEigLocal, mpi::MAX, x.Comm() );
}

#define PROTO(Real) \
  template void MaxEig \
  ( const Matrix<Real>& x, \
          Matrix<Real>& maxEigs, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void MaxEig \
  ( const ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& maxEigs, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, Int cutoff ); \
  template void MaxEig \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& maxEigs, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff ); \
  template Real MaxEig \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Real MaxEig \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, Int cutoff ); \
  template Real MaxEig \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace soc
} // namespace El
