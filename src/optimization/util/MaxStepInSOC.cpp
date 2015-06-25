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
// Alternatively [1, pg. 23], since x is assumed to be a member of the SOC, one 
// can compute an automorphism of the cone which maps x to e, as
// 
//   Q_{x^{-1/2}) x = e,
//
// where Q_z is the quadratic representation of the Jordan algebra member z,
//
//   Q_z = 2 z z^T - det(z) R.
//
// Then, defining G = Q_{x^{-1/2}},
//
//   x + alpha y >= 0 iff G (x + alpha y) = e + alpha (G y) >= 0.
//
// We can then recognize that [2]
//
//   max { t >= 0 | e + t (G y) >= 0 } = max { t >= 0 | e / t + G y >= 0 },
//                                     = min { t >= 0 | e * t + G y >= 0 },
//
// where the result is simply equal to max(|| z_1 ||_2 - z_0,0) if we define 
// z = G y. Then we need only consider two cases:
//
//   If max(|| z_1 ||_2 - z_0,0) == 0, return upperBound
//
//   Otherwise, return min(1 / max(|| z_1 ||_2 - z_0,0),upperBound).
//
// [1] L. Vandenberghe, "The CVXOPT linear and quadratic cone program solvers",
//     2010. Last accessed from 
//     http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
//
// [2] M. Andersen, J. Dahl, and L. Vandenberghe, CVXOPT function misc.max_step,
//     2014. Last accessed from 
//     https://github.com/cvxopt/cvxopt/blob/f3ca94fb997979a54b913f95b816132f7fd44820/src/python/misc.py#L1018
// 

namespace {

template<typename Real>
Real ChooseStepLength
( Real x0, Real y0, Real xDet, Real yDet, Real xTRy, Real upperBound, 
  Real delta=Epsilon<Real>() )
{
    DEBUG_ONLY(CSE cse("ChooseStepLength"))
    Real step;
    if( y0 >= Real(0) && yDet >= Real(0) ) 
    {
        step = upperBound;
    }
    else if( Abs(yDet) <= delta )
    {
        // Fall back to a backstepping line search rather than using the 
        // alpha^2 = 0 approximation alpha = - 2 det(x) / (x^T R y),
        // which has been observed to, in some cases, return 0 instead of the
        // upper bound.
        Real stepRatio = 0.99;
        step = upperBound;
        while( step*step*yDet + 2*step*xTRy + xDet <= 0 || x0+step*y0 <= 0 )
            step *= stepRatio;
    }
    else
    {
        Real discrim = Max(xTRy*xTRy-xDet*yDet,Real(0));
        Real sqrtDiscrim = Sqrt(discrim);
        Real plusRoot = (-xTRy+sqrtDiscrim)/yDet;
        Real minusRoot = (-xTRy-sqrtDiscrim)/yDet;
        Real minRoot = Min(plusRoot,minusRoot);
        Real maxRoot = Max(plusRoot,minusRoot);
        if( minRoot >= Real(0) )
            step = minRoot;
        else
            step = maxRoot;
    }
    step = Max(step,Real(0));
    step = Min(step,upperBound);
    return step;
}

} // anonymous namespace

template<typename Real>
Real MaxStepInSOC
( const Matrix<Real>& x, const Matrix<Real>& y, 
  const Matrix<Int>& orders, const Matrix<Int>& firstInds,
  Real upperBound )
{
    DEBUG_ONLY(CSE cse("MaxStepInSOC"))
    typedef Promote<Real> PReal;

    Matrix<PReal> xProm, yProm;
    Copy( x, xProm );
    Copy( y, yProm );

    Matrix<PReal> xDets, yDets, xTRys, maxSteps;
    SOCDets( xProm, xDets, orders, firstInds );
    SOCDets( yProm, yDets, orders, firstInds );

    auto Ry = yProm;
    SOCReflect( Ry, orders, firstInds );
    SOCDots( xProm, Ry, xTRys, orders, firstInds );

    PReal alpha = upperBound;
    const Int height = x.Height();
    for( Int i=0; i<height; )
    {
        const Int order = orders.Get(i,0);
        const Int firstInd = firstInds.Get(i,0);
        if( i != firstInd )
            LogicError("Inconsistency in orders and firstInds");

        const PReal x0 = xProm.Get(i,0);
        const PReal y0 = yProm.Get(i,0);
        const PReal xDet = xDets.Get(i,0);
        const PReal yDet = yDets.Get(i,0);
        const PReal xTRy = xTRys.Get(i,0);

        alpha = ChooseStepLength(x0,y0,xDet,yDet,xTRy,alpha);

        i += order;
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
    typedef Promote<Real> PReal;

    ProxyCtrl control;
    control.colConstrain = true;
    control.colAlign = 0;
    auto xPtr = ReadProxy<PReal,VC,STAR>(&xPre,control);
    auto yPtr = ReadProxy<PReal,VC,STAR>(&yPre,control);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,control);
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,control);
    auto& x = *xPtr;
    auto& y = *yPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    const Grid& g = x.Grid();
    DistMatrix<PReal,VC,STAR> xDets(g), yDets(g), xTRys(g);
    SOCDets( x, xDets, orders, firstInds, cutoff );
    SOCDets( y, yDets, orders, firstInds, cutoff );

    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    SOCDots( x, Ry, xTRys, orders, firstInds, cutoff );

    PReal alpha = upperBound;
    const Int localHeight = x.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i != firstInds.GetLocal(iLoc,0) )
            continue;

        const PReal x0 = x.GetLocal(iLoc,0);
        const PReal y0 = y.GetLocal(iLoc,0);
        const PReal xDet = xDets.GetLocal(iLoc,0);
        const PReal yDet = yDets.GetLocal(iLoc,0);
        const PReal xTRy = xTRys.GetLocal(iLoc,0);
        
        alpha = ChooseStepLength(x0,y0,xDet,yDet,xTRy,alpha);
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
    typedef Promote<Real> PReal;
    mpi::Comm comm = x.Comm();

    DistMultiVec<PReal> xProm(comm), yProm(comm),
                        xDets(comm), yDets(comm), xTRys(comm), maxSteps(comm);
    Copy( x, xProm );
    Copy( y, yProm );
    SOCDets( xProm, xDets, orders, firstInds, cutoff );
    SOCDets( yProm, yDets, orders, firstInds, cutoff );

    auto Ry = yProm;
    SOCReflect( Ry, orders, firstInds );
    SOCDots( xProm, Ry, xTRys, orders, firstInds, cutoff );

    PReal alpha = upperBound;
    const Int localHeight = xProm.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = xProm.GlobalRow(iLoc);
        if( i != firstInds.GetLocal(iLoc,0) )
            continue;
        
        const PReal x0 = xProm.GetLocal(iLoc,0);
        const PReal y0 = yProm.GetLocal(iLoc,0);
        const PReal xDet = xDets.GetLocal(iLoc,0);
        const PReal yDet = yDets.GetLocal(iLoc,0);
        const PReal xTRy = xTRys.GetLocal(iLoc,0);

        alpha = ChooseStepLength(x0,y0,xDet,yDet,xTRy,alpha);
    }
    return mpi::AllReduce( alpha, mpi::MIN, comm );
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
