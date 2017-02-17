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

// inv(x) = (R x) / det(x)

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void Inverse
( const Matrix<Real>& x,
        Matrix<Real>& xInv,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    EL_DEBUG_CSE

    Matrix<Real> dInv;
    soc::Dets( x, dInv, orders, firstInds );
    cone::Broadcast( dInv, orders, firstInds );
    auto entryInv = []( const Real& alpha ) { return Real(1)/alpha; };
    EntrywiseMap( dInv, MakeFunction(entryInv) );

    auto Rx = x;
    soc::Reflect( Rx, orders, firstInds );

    Hadamard( dInv, Rx, xInv );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void Inverse
( const AbstractDistMatrix<Real>& xPre,
        AbstractDistMatrix<Real>& xInvPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    EL_DEBUG_CSE
    AssertSameGrids( xPre, xInvPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixWriteProxy<Real,Real,VC,STAR>
      xInvProx( xInvPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& xInv = xInvProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    DistMatrix<Real,VC,STAR> dInv(x.Grid());
    soc::Dets( x, dInv, orders, firstInds, cutoff );
    cone::Broadcast( dInv, orders, firstInds );
    auto entryInv = []( const Real& alpha ) { return Real(1)/alpha; };
    EntrywiseMap( dInv, MakeFunction(entryInv) );

    auto Rx = x;
    soc::Reflect( Rx, orders, firstInds );

    Hadamard( dInv, Rx, xInv );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void Inverse
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& xInv,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    EL_DEBUG_CSE

    DistMultiVec<Real> dInv(x.Grid());
    soc::Dets( x, dInv, orders, firstInds, cutoff );
    cone::Broadcast( dInv, orders, firstInds );
    auto entryInv = []( const Real& alpha ) { return Real(1)/alpha; };
    EntrywiseMap( dInv, MakeFunction(entryInv) );

    auto Rx = x;
    soc::Reflect( Rx, orders, firstInds );

    Hadamard( dInv, Rx, xInv );
}

#define PROTO(Real) \
  template void Inverse \
  ( const Matrix<Real>& x, \
          Matrix<Real>& xInv, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void Inverse \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& xInv, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void Inverse \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& xInv, \
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
