/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace soc {

// x o y = [ x^T y; x0 y1 + y0 x1 ]

template<typename Real,typename>
void Apply
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_CSE
    soc::Dots( x, y, z, orders, firstInds );
    auto xRoots = x;
    auto yRoots = y;
    cone::Broadcast( xRoots, orders, firstInds );
    cone::Broadcast( yRoots, orders, firstInds );

    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
        if( i != firstInds(i) )
            z(i) += xRoots(i)*y(i) + yRoots(i)*x(i);
}

template<typename Real,typename>
void Apply
( const ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Real>& yPre,
        ElementalMatrix<Real>& zPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_CSE
    AssertSameGrids( xPre, yPre, zPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl ),
      yProx( yPre, ctrl );
    DistMatrixWriteProxy<Real,Real,VC,STAR>
      zProx( zPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& y = yProx.GetLocked();
    auto& z = zProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    soc::Dots( x, y, z, orders, firstInds );
    auto xRoots = x;
    auto yRoots = y;
    cone::Broadcast( xRoots, orders, firstInds );
    cone::Broadcast( yRoots, orders, firstInds );

    const Int localHeight = x.LocalHeight();
    const Real* xBuf     = x.LockedBuffer();
    const Real* xRootBuf = xRoots.LockedBuffer();
    const Real* yBuf     = y.LockedBuffer();
    const Real* yRootBuf = yRoots.LockedBuffer();
          Real* zBuf     = z.Buffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Int firstInd = firstIndBuf[iLoc];
        if( i != firstInd )
            zBuf[iLoc] += xRootBuf[iLoc]*yBuf[iLoc] + yRootBuf[iLoc]*xBuf[iLoc];
    }
}

template<typename Real,typename>
void Apply
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_CSE
    soc::Dots( x, y, z, orders, firstInds );
    auto xRoots = x;
    auto yRoots = y;
    cone::Broadcast( xRoots, orders, firstInds );
    cone::Broadcast( yRoots, orders, firstInds );

    const Int firstLocalRow = x.FirstLocalRow();
    const Int localHeight = x.LocalHeight();
    const Real* xBuf     = x.LockedMatrix().LockedBuffer();
    const Real* xRootBuf = xRoots.LockedMatrix().LockedBuffer();
    const Real* yBuf     = y.LockedMatrix().LockedBuffer();
    const Real* yRootBuf = yRoots.LockedMatrix().LockedBuffer();
          Real* zBuf     = z.Matrix().Buffer();
    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = iLoc + firstLocalRow;
        const Int firstInd = firstIndBuf[iLoc];
        if( i != firstInd )
            zBuf[iLoc] += xRootBuf[iLoc]*yBuf[iLoc] + yRootBuf[iLoc]*xBuf[iLoc];
    }
}

template<typename Real,typename>
void Apply
( const Matrix<Real>& x, 
        Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_CSE
    // TODO?: Optimize
    Matrix<Real> z;
    soc::Apply( x, y, z, orders, firstInds );
    y = z;
}

template<typename Real,typename>
void Apply
( const ElementalMatrix<Real>& x, 
        ElementalMatrix<Real>& y,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Int cutoff )
{
    DEBUG_CSE
    // TODO?: Optimize
    DistMatrix<Real,VC,STAR> z(x.Grid());
    soc::Apply( x, y, z, orders, firstInds, cutoff );
    y = z;
}

template<typename Real,typename>
void Apply
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_CSE
    // TODO?: Optimize
    DistMultiVec<Real> z(x.Comm());
    soc::Apply( x, y, z, orders, firstInds, cutoff );
    y = z;
}

#define PROTO(Real) \
  template void Apply \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& y, \
          Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Apply \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void Apply \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff ); \
  template void Apply \
  ( const Matrix<Real>& x, \
          Matrix<Real>& y, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Apply \
  ( const ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void Apply \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
