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

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real ChooseStepLength
( const Real& x0,
  const Real& y0,
  const Real& xDet,
  const Real& yDet,
  const Real& xTRy,
  const Real& upperBound, 
  const Real& delta=limits::Epsilon<Real>() )
{
    DEBUG_CSE
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

template<typename Real,typename>
Real MaxStep
( const Matrix<Real>& x,
  const Matrix<Real>& y, 
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  Real upperBound )
{
    DEBUG_CSE
    typedef Promote<Real> PReal;
    const Int height = x.Height();

    Matrix<PReal> xProm, yProm;
    Copy( x, xProm );
    Copy( y, yProm );

    Matrix<PReal> xDets, yDets, xTRys, maxSteps;
    soc::Dets( xProm, xDets, orders, firstInds );
    soc::Dets( yProm, yDets, orders, firstInds );

    auto Ry = yProm;
    soc::Reflect( Ry, orders, firstInds );
    soc::Dots( xProm, Ry, xTRys, orders, firstInds );

    const Int* orderBuf = orders.LockedBuffer();
    const Int* firstIndBuf = firstInds.LockedBuffer();
    const PReal* xBuf = xProm.LockedBuffer();
    const PReal* yBuf = yProm.LockedBuffer();
    const PReal* xDetBuf = xDets.LockedBuffer();
    const PReal* yDetBuf = yDets.LockedBuffer();
    const PReal* xTRyBuf = xTRys.LockedBuffer();

    PReal alpha = upperBound;
    for( Int i=0; i<height; )
    {
        const Int order = orderBuf[i];
        const Int firstInd = firstIndBuf[i];
        DEBUG_ONLY(
          if( i != firstInd )
              LogicError("Inconsistency in orders and firstInds");
        )

        const PReal x0 = xBuf[i];
        const PReal y0 = yBuf[i];
        const PReal xDet = xDetBuf[i];
        const PReal yDet = yDetBuf[i];
        const PReal xTRy = xTRyBuf[i];

        alpha = ChooseStepLength(x0,y0,xDet,yDet,xTRy,alpha);

        i += order;
    }
    return Real(alpha);
}

template<typename Real,typename>
Real MaxStep
( const ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Real>& yPre, 
  const ElementalMatrix<Int>& ordersPre,
  const ElementalMatrix<Int>& firstIndsPre,
  Real upperBound, Int cutoff )
{
    DEBUG_CSE
    typedef Promote<Real> PReal;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.colAlign = 0;

    DistMatrixReadProxy<Real,PReal,VC,STAR>
      xProx( xPre, control ),
      yProx( yPre, control );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, control ),
      firstIndsProx( firstIndsPre, control );
    auto& x = xProx.GetLocked();
    auto& y = yProx.GetLocked();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    const Grid& g = x.Grid();
    DistMatrix<PReal,VC,STAR> xDets(g), yDets(g), xTRys(g);
    soc::Dets( x, xDets, orders, firstInds, cutoff );
    soc::Dets( y, yDets, orders, firstInds, cutoff );

    auto Ry = y;
    soc::Reflect( Ry, orders, firstInds );
    soc::Dots( x, Ry, xTRys, orders, firstInds, cutoff );

    const Int localHeight = x.LocalHeight();
    const Int* firstIndBuf = firstInds.LockedBuffer();
    const PReal* xBuf = x.LockedBuffer();
    const PReal* yBuf = y.LockedBuffer();
    const PReal* xDetBuf = xDets.LockedBuffer();
    const PReal* yDetBuf = yDets.LockedBuffer();
    const PReal* xTRyBuf = xTRys.LockedBuffer();

    PReal alpha = upperBound;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = x.GlobalRow(iLoc);
        if( i != firstIndBuf[iLoc] )
            continue;

        const PReal x0 = xBuf[iLoc];
        const PReal y0 = yBuf[iLoc];
        const PReal xDet = xDetBuf[iLoc];
        const PReal yDet = yDetBuf[iLoc];
        const PReal xTRy = xTRyBuf[iLoc];
        
        alpha = ChooseStepLength(x0,y0,xDet,yDet,xTRy,alpha);
    }

    return Real(mpi::AllReduce( alpha, mpi::MIN, x.DistComm() ));
}

template<typename Real,typename>
Real MaxStep
( const DistMultiVec<Real>& x,     
  const DistMultiVec<Real>& y, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Real upperBound, Int cutoff )
{
    DEBUG_CSE
    typedef Promote<Real> PReal;
    mpi::Comm comm = x.Comm();

    DistMultiVec<PReal> xProm(comm), yProm(comm),
                        xDets(comm), yDets(comm), xTRys(comm), maxSteps(comm);
    Copy( x, xProm );
    Copy( y, yProm );
    soc::Dets( xProm, xDets, orders, firstInds, cutoff );
    soc::Dets( yProm, yDets, orders, firstInds, cutoff );

    auto Ry = yProm;
    soc::Reflect( Ry, orders, firstInds );
    soc::Dots( xProm, Ry, xTRys, orders, firstInds, cutoff );

    const Int localHeight = xProm.LocalHeight();
    const Int firstLocalRow = xProm.FirstLocalRow();

    const Int* firstIndBuf = firstInds.LockedMatrix().LockedBuffer();
    const PReal* xBuf = xProm.LockedMatrix().LockedBuffer();
    const PReal* yBuf = yProm.LockedMatrix().LockedBuffer();
    const PReal* xDetBuf = xDets.LockedMatrix().LockedBuffer();
    const PReal* yDetBuf = yDets.LockedMatrix().LockedBuffer();
    const PReal* xTRyBuf = xTRys.LockedMatrix().LockedBuffer();

    PReal alpha = upperBound;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = iLoc + firstLocalRow;
        if( i != firstIndBuf[iLoc] )
            continue;

        const PReal x0 = xBuf[iLoc];
        const PReal y0 = yBuf[iLoc];
        const PReal xDet = xDetBuf[iLoc];
        const PReal yDet = yDetBuf[iLoc];
        const PReal xTRy = xTRyBuf[iLoc];

        alpha = ChooseStepLength(x0,y0,xDet,yDet,xTRy,alpha);
    }
    return Real(mpi::AllReduce( alpha, mpi::MIN, comm ));
}

#define PROTO(Real) \
  template Real MaxStep \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& ds, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    Real upperBound ); \
  template Real MaxStep \
  ( const ElementalMatrix<Real>& s, \
    const ElementalMatrix<Real>& ds, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Real upperBound, Int cutoff ); \
  template Real MaxStep \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& ds, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Real upperBound, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
