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
void MinEig
( const Matrix<Real>& x, 
        Matrix<Real>& minEigs,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    soc::LowerNorms( x, minEigs, orders, firstInds );
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            minEigs.Set( i, 0, x.Get(i,0)-minEigs.Get(i,0) );
}

template<typename Real>
void MinEig
( const ElementalMatrix<Real>& xPre, 
        ElementalMatrix<Real>& minEigsPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
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

    soc::LowerNorms( x, minEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)-minEigs.GetLocal(iLoc,0) );
}

template<typename Real>
void MinEig
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& minEigs,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
        LogicError("x, orders, and firstInds should be column vectors");
    if( orders.Height() != height || firstInds.Height() != height )
        LogicError("orders and firstInds should be of the same height as x");

    soc::LowerNorms( x, minEigs, orders, firstInds, cutoff );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigs.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)-minEigs.GetLocal(iLoc,0) );
}

template<typename Real>
Real MinEig
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    Matrix<Real> minEigs;
    soc::MinEig( x, minEigs, orders, firstInds );

    Real minEig = std::numeric_limits<Real>::max();
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstInds.Get(i,0) ) 
            minEig = Min(minEigs.Get(i,0),minEig);
    return minEig;
}

template<typename Real>
Real MinEig
( const ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
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
    soc::MinEig( x, minEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real minEigLocal = std::numeric_limits<Real>::max();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigLocal = Min(minEigLocal,minEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( minEigLocal, mpi::MIN, x.DistComm() );
}

template<typename Real>
Real MinEig
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    DistMultiVec<Real> minEigs(x.Comm());
    soc::MinEig( x, minEigs, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    Real minEigLocal = std::numeric_limits<Real>::max();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstInds.GetLocal(iLoc,0) )
            minEigLocal = Min(minEigLocal,minEigs.GetLocal(iLoc,0));
    return mpi::AllReduce( minEigLocal, mpi::MIN, x.Comm() );
}

#define PROTO(Real) \
  template void MinEig \
  ( const Matrix<Real>& x, \
          Matrix<Real>& minEigs, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void MinEig \
  ( const ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& minEigs, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, Int cutoff ); \
  template void MinEig \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& minEigs, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff ); \
  template Real MinEig \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Real MinEig \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, Int cutoff ); \
  template Real MinEig \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace soc
} // namespace El
