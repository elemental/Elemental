/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Members of second-order cones are stored contiguously within the column
// vector x, with the corresponding order of the cone each member belongs to
// stored in the same index of 'order', and the first index of the cone 
// being listed in the same index of 'firstInd'.
template<typename Real>
Int NumNonSOC
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("NumNonSOC"))
    Matrix<Real> d;
    SOCDets( x, d, orders, firstInds );

    Int numNonSO = 0;
    const Int height = x.Height();
    for( Int i=0; i<height; )
    {
        const Int order = orders.Get(i,0);
        const Int firstInd = firstInds.Get(i,0);
        if( i != firstInd )
            LogicError("Inconsistency in orders and firstInds");
        const Int det = d.Get(i,0);
        if( det < Real(0) )
            ++numNonSO;
        i += order;
    }
    return numNonSO;
}

template<typename Real>
Int NumNonSOC
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("NumNonSOC"))
    AssertSameGrids( xPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> d(x.Grid());
    SOCDets( x, d, orders, firstInds, cutoff );

    Int numLocalNonSOC = 0;
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstInds.GetLocal(iLoc,0) && d.GetLocal(iLoc,0) < Real(0) )
            ++numLocalNonSOC;
    }
    return mpi::AllReduce( numLocalNonSOC, x.DistComm() );
}

template<typename Real>
Int NumNonSOC
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("NumNonSOC"))

    DistMultiVec<Real> d(x.Comm());
    SOCDets( x, d, orders, firstInds, cutoff );

    Int numLocalNonSOC = 0;
    const int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i == firstInds.GetLocal(iLoc,0) && d.GetLocal(iLoc,0) < Real(0) )
            ++numLocalNonSOC;
    }
    return mpi::AllReduce( numLocalNonSOC, x.Comm() );
}

#define PROTO(Real) \
  template Int NumNonSOC \
  ( const Matrix<Real>& x, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template Int NumNonSOC \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template Int NumNonSOC \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
