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

template<typename Real,typename>
void PushInto
(       Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  Real minDist )
{
    DEBUG_CSE

    Matrix<Real> d;
    soc::LowerNorms( x, d, orders, firstInds );

    const Int height = x.Height();
    for( Int i=0; i<height; ++i )
    {
        Real& x0 = x(i);
        const Real lowerNorm = d(i);
        const Int firstInd = firstInds(i);
        if( i == firstInd && x0-lowerNorm < minDist )
            x0 = minDist + lowerNorm;
    }
}

template<typename Real,typename>
void PushInto
(       ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Real minDist, Int cutoff )
{
    DEBUG_CSE
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadWriteProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    DistMatrix<Real,VC,STAR> d(x.Grid());
    soc::LowerNorms( x, d, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    auto& xLoc = x.Matrix();
    auto& dLoc = d.LockedMatrix();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        Real& x0 = xLoc(iLoc);
        const Real lowerNorm = dLoc(iLoc);
        if( i == firstIndsLoc(iLoc) && x0-lowerNorm < minDist )
            x0 = minDist + lowerNorm;
    }
}

template<typename Real,typename>
void PushInto
(       DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, 
  Real minDist, Int cutoff )
{
    DEBUG_CSE

    DistMultiVec<Real> d(x.Comm());
    soc::LowerNorms( x, d, orders, firstInds, cutoff );

    const int localHeight = x.LocalHeight();
    auto& xLoc = x.Matrix();
    auto& dLoc = d.LockedMatrix();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        Real& x0 = xLoc(iLoc);
        const Real lowerNorm = dLoc(iLoc);
        if( i == firstIndsLoc(iLoc) && x0-lowerNorm < minDist )
            x0 = minDist + lowerNorm;
    }
}

#define PROTO(Real) \
  template void PushInto \
  (       Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    Real minDist ); \
  template void PushInto \
  (       ElementalMatrix<Real>& x, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Real minDist, Int cutoff ); \
  template void PushInto \
  (       DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Real minDist, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
