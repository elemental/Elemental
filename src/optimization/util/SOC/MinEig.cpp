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

template<typename Real,typename>
void MinEig
( const Matrix<Real>& x, 
        Matrix<Real>& minEigs,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    soc::LowerNorms( x, minEigs, orders, firstInds );

    const Int height = x.Height();

          Real* minEigBuf = minEigs.Buffer();
    const Real* xBuf = x.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    for( Int i=0; i<height; ++i )
        if( i == firstIndBuf[i] ) 
            minEigBuf[i] = xBuf[i]-minEigBuf[i];
}

template<typename Real,typename>
void MinEig
( const ElementalMatrix<Real>& xPre, 
        ElementalMatrix<Real>& minEigsPre,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    AssertSameGrids( xPre, minEigsPre, orders, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr         = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto minEigsPtr   = WriteProxy<Real,VC,STAR>(&minEigsPre,ctrl);
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& minEigs = *minEigsPtr;
    auto& firstInds = *firstIndsPtr;

    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    soc::LowerNorms( x, minEigs, orders, firstInds, cutoff );

          Real* minEigBuf = minEigs.Buffer();
    const Real* xBuf = x.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            minEigBuf[iLoc] = xBuf[iLoc] - minEigBuf[iLoc];
}

template<typename Real,typename>
void MinEig
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& minEigs,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    soc::LowerNorms( x, minEigs, orders, firstInds, cutoff );

          Real* minEigBuf = minEigs.Matrix().Buffer();
    const Real* xBuf = x.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            minEigBuf[iLoc] = xBuf[iLoc] - minEigBuf[iLoc];
}

template<typename Real,typename>
Real MinEig
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    Matrix<Real> minEigs;
    soc::MinEig( x, minEigs, orders, firstInds );

    const Real* minEigBuf = minEigs.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    Real minEig = std::numeric_limits<Real>::max();
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstIndBuf[i] ) 
            minEig = Min(minEigBuf[i],minEig);
    return minEig;
}

template<typename Real,typename>
Real MinEig
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    AssertSameGrids( x, orders, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> minEigs(x.Grid());
    soc::MinEig( x, minEigs, orders, firstInds, cutoff );

    const Real* minEigBuf = minEigs.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    const Int localHeight = x.LocalHeight();
    Real minEigLocal = std::numeric_limits<Real>::max();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            minEigLocal = Min(minEigLocal,minEigBuf[iLoc]);
    return mpi::AllReduce( minEigLocal, mpi::MIN, x.DistComm() );
}

template<typename Real,typename>
Real MinEig
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MinEig"))
    DistMultiVec<Real> minEigs(x.Comm());
    soc::MinEig( x, minEigs, orders, firstInds, cutoff );

    const Real* minEigBuf = minEigs.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    const Int localHeight = x.LocalHeight();
    Real minEigLocal = std::numeric_limits<Real>::max();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( minEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            minEigLocal = Min(minEigLocal,minEigBuf[iLoc]);
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
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void MinEig \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& minEigs, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff ); \
  template Real MinEig \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Real MinEig \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template Real MinEig \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace soc
} // namespace El
