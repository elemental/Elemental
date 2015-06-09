/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Given x within a Second-Order Cone and an arbitrary search direction y in 
// R^n, find the maximum value of alpha* in [0,upperBound] such that 
// x + t y lies within the cone for all values of t in [0,alpha*].
//
// Thus, we seek the minimum positive solution to the quadratic equation
//   
//   det(x + alpha y) = det(y) alpha^2 + 2 (x^T R y) alpha + det(x) = 0,
//
// which is degenerate in the case where y lies on the boundary of either the
// second-order cone or its negative reflection, in which case we may simply
// set
//
//   alpha = -det(x) / (2 x^T R y).
//
// It also happens that, when y lies within the reflection of the orthogonal
// complement of x, that x^T R y = 0, and the only positive solution to the 
// original formula simplifies to
//
//   alpha = sqrt(-det(x)/det(y)).
//
// In the simple case where y lies within the second-order cone, we are free
// to make alpha as large as possible.
//
// In all cases, after finding the maximum unconstrained step-length, alpha,
// we set 
//
//   alpha* = min(alpha,upperBound).
//
// Our algorithm for each subcone, given a predefined tolerance for inversion,
// delta, can be summarized via the steps:
//
//   Compute det(x), det(y), and x^T R y.
//
//   If y_0 >= 0 and det(y) >= 0, 
//
//     Set alpha* = upperBound.
//
//   Else if |det(y)| <= delta,
//
//     Set alpha* = min(-det(x)/(2 x^T R y),upperBound).
//
//   Else
//
//     plusRoot := (-(x^T R y) + sqrt((x^T R y)^2 - det(x) det(y)))/det(y),
//     minusRoot := (-(x^T R y) - sqrt((x^T R y)^2 - det(x) det(y)))/det(y),
//     minRoot := min(plusRoot,minusRoot)
//     maxRoot := max(plusRoot,minusRoot)
//     
//     If minRoot >= 0:
//
//       alpha* = min(minRoot,upperBound)
//     
//     Else:
//
//       alpha* = min(maxRoot,upperBound)
//
//     End
//
//   End
//
// In order to compute the step-length for a product cone, the minimum over 
// the set of subcone step-lengths should be taken.
// 

namespace {

template<typename Real>
Real ChooseStepLength
( Real y0, Real xDet, Real yDet, Real xTRy, Real upperBound, 
  Real delta=lapack::MachineSafeMin<Real>() )
{
    DEBUG_ONLY(CSE cse("ChooseStepLength"))
    if( y0 >= Real(0) && yDet >= Real(0) ) 
    {
        return upperBound;
    }
    else if( Abs(yDet) <= delta )
    {
        return Min(-xDet/(2*xTRy),upperBound);
    }
    else
    {
        Real sqrtDiscrim = Sqrt((xTRy)*(xTRy)-xDet*yDet);
        Real plusRoot = (-xTRy+sqrtDiscrim)/yDet;
        Real minusRoot = (-xTRy-sqrtDiscrim)/yDet;
        Real minRoot = Min(plusRoot,minusRoot);
        Real maxRoot = Max(plusRoot,minusRoot);
        if( minRoot >= Real(0) )
            return Min(minRoot,upperBound);
        else
            return Min(maxRoot,upperBound);
    }
}

} // anonymous namespace

template<typename Real>
Real MaxStepInSOC
( const Matrix<Real>& x, const Matrix<Real>& y, 
  const Matrix<Int>& orders, const Matrix<Int>& firstInds,
  Real upperBound )
{
    DEBUG_ONLY(CSE cse("MaxStepInSOC"))

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

        const Real y0 = y.Get(i,0);
        const Real xDet = xDets.Get(i,0);
        const Real yDet = yDets.Get(i,0);
        const Real xTRy = xTRys.Get(i,0);

        alpha = ChooseStepLength(y0,xDet,yDet,xTRy,alpha);

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

        const Real y0 = y.GetLocal(iLoc,0);
        const Real xDet = xDets.GetLocal(iLoc,0);
        const Real yDet = yDets.GetLocal(iLoc,0);
        const Real xTRy = xTRys.GetLocal(iLoc,0);
        
        alpha = ChooseStepLength(y0,xDet,yDet,xTRy,alpha);
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
        
        const Real y0 = y.GetLocal(iLoc,0);
        const Real xDet = xDets.GetLocal(iLoc,0);
        const Real yDet = yDets.GetLocal(iLoc,0);
        const Real xTRy = xTRys.GetLocal(iLoc,0);

        alpha = ChooseStepLength(y0,xDet,yDet,xTRy,alpha);
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
