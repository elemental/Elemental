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

// Members of second-order cones are stored contiguously within the column
// vector x, with the corresponding order of the cone each member belongs to
// stored in the same index of 'order', and the first index of the cone
// being listed in the same index of 'firstInd'.
template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Int NumOutside
( const Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    EL_DEBUG_CSE
    Matrix<Real> d;
    soc::Dets( x, d, orders, firstInds );

    Int numNonSO = 0;
    const Int height = x.Height();
    for( Int i=0; i<height; )
    {
        const Int order = orders(i);
        const Int firstInd = firstInds(i);
        if( i != firstInd )
            LogicError("Inconsistency in orders and firstInds");
        const Real det = d(i);
        if( det < Real(0) )
            ++numNonSO;
        i += order;
    }
    return numNonSO;
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Int NumOutside
( const AbstractDistMatrix<Real>& xPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    EL_DEBUG_CSE
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    DistMatrix<Real,VC,STAR> d(x.Grid());
    soc::Dets( x, d, orders, firstInds, cutoff );

    Int numLocalNonSOC = 0;
    const Int localHeight = x.LocalHeight();
    auto& dLoc = d.LockedMatrix();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstIndsLoc(iLoc) && dLoc(iLoc) < Real(0) )
            ++numLocalNonSOC;
    }
    return mpi::AllReduce( numLocalNonSOC, x.DistComm() );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Int NumOutside
( const DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    EL_DEBUG_CSE
    const Grid& grid = x.Grid();

    DistMultiVec<Real> d(grid);
    soc::Dets( x, d, orders, firstInds, cutoff );

    Int numLocalNonSOC = 0;
    const int localHeight = x.LocalHeight();
    auto& dLoc = d.LockedMatrix();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstIndsLoc(iLoc) && dLoc(iLoc) < Real(0) )
            ++numLocalNonSOC;
    }
    return mpi::AllReduce( numLocalNonSOC, grid.Comm() );
}

#define PROTO(Real) \
  template Int NumOutside \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Int NumOutside \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Int cutoff ); \
  template Int NumOutside \
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
#include <El/macros/Instantiate.h>
} // namespace soc
} // namespace El
