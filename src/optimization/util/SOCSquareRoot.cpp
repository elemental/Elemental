/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// sqrt(x) = [ eta_0; x_1/(2 eta_0) ],
// where eta_0 = sqrt(x_0 + sqrt(det(x))) / sqrt(2).

template<typename Real>
void SOCSquareRoot
( const Matrix<Real>& x, 
        Matrix<Real>& xRoot,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCSquareRoot"))

    Matrix<Real> d;
    SOCDets( x, d, orders, firstInds );
    ConeBroadcast( d, orders, firstInds );

    const Int height = x.Height();
    Zeros( xRoot, height, 1 );
    for( Int i=0; i<height; )
    {
        const Int order = orders.Get(i,0);
        const Int firstInd = firstInds.Get(i,0);
        if( i != firstInd )       
            LogicError("Inconsistency in orders and firstInds");

        const Real eta0 = Sqrt(x.Get(i,0)+Sqrt(d.Get(i,0)))/Sqrt(Real(2));
        xRoot.Set( i, 0, eta0 );
        for( Int k=1; k<order; ++k )
            xRoot.Set( i+k, 0, x.Get(i+k,0)/(2*eta0) );

        i += order;
    }
}

template<typename Real>
void SOCSquareRoot
( const AbstractDistMatrix<Real>& xPre, 
        AbstractDistMatrix<Real>& xRootPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCSquareRoot"))
    AssertSameGrids( xPre, xRootPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto xRootPtr = WriteProxy<Real,VC,STAR>(&xRootPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& xRoot = *xRootPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> d(x.Grid());
    SOCDets( x, d, orders, firstInds );
    ConeBroadcast( d, orders, firstInds );

    auto roots = x;
    ConeBroadcast( roots, orders, firstInds );

    const Int localHeight = x.LocalHeight();
    xRoot.SetGrid( x.Grid() );
    Zeros( xRoot, x.Height(), 1 );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Real x0 = roots.GetLocal(iLoc,0);
        const Real det = d.GetLocal(iLoc,0);
        const Real eta0 = Sqrt(x0+Sqrt(det))/Sqrt(Real(2));
        if( i == firstInds.GetLocal(iLoc,0) )
            xRoot.SetLocal( iLoc, 0, eta0 );
        else
            xRoot.SetLocal( iLoc, 0, x.GetLocal(iLoc,0)/(2*eta0) );
    }
}

template<typename Real>
void SOCSquareRoot
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& xRoot,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCSquareRoot"))

    DistMultiVec<Real> d(x.Comm());
    SOCDets( x, d, orders, firstInds );
    ConeBroadcast( d, orders, firstInds );

    auto roots = x;
    ConeBroadcast( roots, orders, firstInds );

    const Int localHeight = x.LocalHeight();
    xRoot.SetComm( x.Comm() );
    Zeros( xRoot, x.Height(), 1 );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        const Real x0 = roots.GetLocal(iLoc,0);
        const Real det = d.GetLocal(iLoc,0);
        const Real eta0 = Sqrt(x0+Sqrt(det))/Sqrt(Real(2));
        if( i == firstInds.GetLocal(iLoc,0) )
            xRoot.SetLocal( iLoc, 0, eta0 );
        else
            xRoot.SetLocal( iLoc, 0, x.GetLocal(iLoc,0)/(2*eta0) );
    }
}

#define PROTO(Real) \
  template void SOCSquareRoot \
  ( const Matrix<Real>& x, \
          Matrix<Real>& xRoot, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCSquareRoot \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& xRoot, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCSquareRoot \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& xRoot, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
