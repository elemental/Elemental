/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace soc {

template<typename Real,typename>
void MaxEig
( const Matrix<Real>& x, 
        Matrix<Real>& maxEigs,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    soc::LowerNorms( x, maxEigs, orders, firstInds );

          Real* maxEigBuf = maxEigs.Buffer();
    const Real* xBuf = x.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstIndBuf[i] ) 
            maxEigBuf[i] = xBuf[i]+maxEigBuf[i];
}

template<typename Real,typename>
void MaxEig
( const ElementalMatrix<Real>& xPre, 
        ElementalMatrix<Real>& maxEigsPre,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    AssertSameGrids( xPre, maxEigsPre, orders, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixWriteProxy<Real,Real,VC,STAR>
      maxEigsProx( maxEigsPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& maxEigs = maxEigsProx.Get();
    auto& firstInds = firstIndsProx.GetLocked();

    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    soc::LowerNorms( x, maxEigs, orders, firstInds, cutoff );

          Real* maxEigBuf = maxEigs.Buffer();
    const Real* xBuf = x.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            maxEigBuf[iLoc] = xBuf[iLoc] + maxEigBuf[iLoc];
}

template<typename Real,typename>
void MaxEig
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& maxEigs,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    const Int height = x.Height();
    const Int localHeight = x.LocalHeight();
    DEBUG_ONLY(
      if( x.Width() != 1 || orders.Width() != 1 || firstInds.Width() != 1 )
          LogicError("x, orders, and firstInds should be column vectors");
      if( orders.Height() != height || firstInds.Height() != height )
          LogicError("orders and firstInds should be of the same height as x");
    )

    soc::LowerNorms( x, maxEigs, orders, firstInds, cutoff );

          Real* maxEigBuf = maxEigs.Matrix().Buffer();
    const Real* xBuf = x.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            maxEigBuf[iLoc] = xBuf[iLoc] + maxEigBuf[iLoc];
}

template<typename Real,typename>
Real MaxEig
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    Matrix<Real> maxEigs;
    soc::MaxEig( x, maxEigs, orders, firstInds );

    const Real* maxEigBuf = maxEigs.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    Real maxEig = limits::Lowest<Real>();
    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i == firstIndBuf[i] ) 
            maxEig = Max(maxEigBuf[i],maxEig);
    return maxEig;
}

template<typename Real,typename>
Real MaxEig
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    AssertSameGrids( x, orders, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Int,Int,VC,STAR> firstIndsProx( firstIndsPre, ctrl );
    auto& firstInds = firstIndsProx.GetLocked();

    DistMatrix<Real,VC,STAR> maxEigs(x.Grid());
    soc::MaxEig( x, maxEigs, orders, firstInds, cutoff );

    const Real* maxEigBuf = maxEigs.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    const Int localHeight = x.LocalHeight();
    Real maxEigLocal = limits::Lowest<Real>();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            maxEigLocal = Max(maxEigLocal,maxEigBuf[iLoc]);
    return mpi::AllReduce( maxEigLocal, mpi::MAX, x.DistComm() );
}

template<typename Real,typename>
Real MaxEig
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::MaxEig"))
    DistMultiVec<Real> maxEigs(x.Comm());
    soc::MaxEig( x, maxEigs, orders, firstInds, cutoff );

    const Real* maxEigBuf = maxEigs.LockedMatrix().LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    const Int localHeight = x.LocalHeight();
    Real maxEigLocal = limits::Lowest<Real>();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        if( maxEigs.GlobalRow(iLoc) == firstIndBuf[iLoc] )
            maxEigLocal = Max(maxEigLocal,maxEigBuf[iLoc]);
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
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void MaxEig \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& maxEigs, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff ); \
  template Real MaxEig \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Real MaxEig \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template Real MaxEig \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace soc
} // namespace El
