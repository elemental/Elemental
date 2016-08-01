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
void Dets
( const Matrix<Real>& x, 
        Matrix<Real>& d,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_CSE
    auto Rx = x;
    soc::Reflect( Rx, orders, firstInds );
    soc::Dots( x, Rx, d, orders, firstInds );
}

template<typename Real,typename>
void Dets
( const ElementalMatrix<Real>& xPre, 
        ElementalMatrix<Real>& dPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_CSE
    AssertSameGrids( xPre, dPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixWriteProxy<Real,Real,VC,STAR>
      dProx( dPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& d = dProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    auto Rx = x;
    soc::Reflect( Rx, orders, firstInds );
    soc::Dots( x, Rx, d, orders, firstInds, cutoff );
}

template<typename Real,typename>
void Dets
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& d,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_CSE
    auto Rx = x;
    soc::Reflect( Rx, orders, firstInds );
    soc::Dots( x, Rx, d, orders, firstInds, cutoff );
}

#define PROTO(Real) \
  template void Dets \
  ( const Matrix<Real>& x, \
          Matrix<Real>& d, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Dets \
  ( const ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& d, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, Int cutoff ); \
  template void Dets \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& d, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
