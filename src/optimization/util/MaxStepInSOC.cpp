/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// NOTE: It is assumed that x is in the SOC, but some minimal branching is 
//       needed based upon whether or not y is also in the cone 
//       (which is determined from the combination of its determinant and the 
//       sign of its first entry).

template<typename Real>
Real MaxStepInSOC
( const Matrix<Real>& x, const Matrix<Real>& y, 
  const Matrix<Int>& orders, const Matrix<Int>& firstInds,
  Real upperBound )
{
    DEBUG_ONLY(CSE cse("MaxStepInSOC"))
    LogicError("This routine is not yet finished");

    Matrix<Real> xDets, yDets, xTRys, maxSteps;
    SOCDets( x, xDets, orders, firstInds );
    SOCDets( y, yDets, orders, firstInds );

    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    SOCDots( x, Ry, xTRys, orders, firstInds );

    Real alpha = upperBound;
    Int i = 0;
    const Int height = x.Height();
    while( i < height )
    {
        if( i != firstInds.Get(i,0) )
            LogicError("Inconsistency in orders and firstInds");

        const Real xDet = xDets.Get(i,0);
        const Real yDet = yDets.Get(i,0);
        const Real xTRy = xTRys.Get(i,0);
 
        Real maxStep;
        if( y.Get(i,0) >= Real(0) && yDet >= Real(0) )
        {
            // y is in the positive cone, so we could step arbitrarily far
            maxStep = upperBound;
        }
        else if( yDet != Real(0) )
        {
            // TODO: Use a tolerance rather than checking strictly for zero.
            //       Perhaps lapack::SafeMin<Real>() would be appropriate.
            // TODO: Compute the square-root more carefully
            Real root = (xTRy - Sqrt(xTRy*xTRy-xDet*yDet))/yDet;
            if( root >= Real(0) )
                maxStep = root;
            else
                maxStep = upperBound;
        }
        else
        {
        
        }
        alpha = Min(alpha,maxStep);

        i += orders.Get(i,0);
    }
    return alpha;
}

template<typename Real>
Real MaxStepInSOC
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Real>& yPre, 
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  Real upperBound, Int cutoff )
{
    DEBUG_ONLY(CSE cse("MaxStepInSOC"))
    LogicError("This routine is not yet finished");

    ProxyCtrl control;
    control.colConstrain = true;
    control.colAlign = 0;
    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,control);
    auto yPtr = ReadProxy<Real,VC,STAR>(&yPre,control);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,control);
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,control);
    auto& x = *xPtr;
    auto& y = *yPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Grid& g = x.Grid();
    DistMatrix<Real,VC,STAR> xDets(g), yDets(g), xTRys(g);
    SOCDets( x, xDets, orders, firstInds, cutoff );
    SOCDets( y, yDets, orders, firstInds, cutoff );

    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    SOCDots( x, Ry, xTRys, orders, firstInds, cutoff );

    Real alpha = upperBound;
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i != firstInds.GetLocal(iLoc,0) )
            continue;

        const Real xDet = xDets.GetLocal(iLoc,0);
        const Real yDet = yDets.GetLocal(iLoc,0);
        const Real xTRy = xTRys.GetLocal(iLoc,0);
        const Real maxStep = (xTRy - Sqrt(xTRy*xTRy-xDet*yDet))/yDet;
        alpha = Min(alpha,maxStep);
    }

    return mpi::AllReduce( alpha, mpi::MIN, x.DistComm() );
}

template<typename Real>
Real MaxStepInSOC
( const DistMultiVec<Real>& x,     
  const DistMultiVec<Real>& y, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Real upperBound, Int cutoff )
{
    DEBUG_ONLY(CSE cse("MaxStepInSOC"))
    LogicError("This routine is not yet finished");

    DistMultiVec<Real> xDets, yDets, xTRys, maxSteps;
    SOCDets( x, xDets, orders, firstInds, cutoff );
    SOCDets( y, yDets, orders, firstInds, cutoff );

    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    SOCDots( x, Ry, xTRys, orders, firstInds, cutoff );

    Real alpha = upperBound;
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i != firstInds.GetLocal(iLoc,0) )
            continue;
        
        const Real xDet = xDets.GetLocal(iLoc,0);
        const Real yDet = yDets.GetLocal(iLoc,0);
        const Real xTRy = xTRys.GetLocal(iLoc,0);
        const Real maxStep = (xTRy - Sqrt(xTRy*xTRy-xDet*yDet))/yDet;
        alpha = Min(alpha,maxStep);
    }
    return mpi::AllReduce( alpha, mpi::MIN, x.Comm() );
}

#define PROTO(Real) \
  template Real MaxStepInSOC \
  ( const Matrix<Real>& s,     const Matrix<Real>& ds, \
    const Matrix<Int>& orders, const Matrix<Int>& firstInds, \
    Real upperBound ); \
  template Real MaxStepInSOC \
  ( const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& ds, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Real upperBound, Int cutoff ); \
  template Real MaxStepInSOC \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& ds, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Real upperBound, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
