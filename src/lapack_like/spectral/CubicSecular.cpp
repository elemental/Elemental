/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// Solve for an inner root of the secular equation
//
//   f(x) = rho + z(0) / (d(0)-x) + z(1) / (d(1)-x) + z(2) / (d(2)-x),
//
// where each numerator is positive and d(0) < d(1) < d(2).
//
// Just as in LAPACK's {s,d}laed6 [CITATION], we require that the user pass in
// an accurate evaluation of f(0).
//
template<typename Real,typename>
CubicSecularInfo
CubicSecular
( bool initialize,
  bool rightRoot,
  const Real& rho,
  const Matrix<Real>& z,
  const Matrix<Real>& d,
  const Real& originEval,
        Real& root,
  const CubicSecularCtrl& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( z.Height() != 3 || z.Width() != 1 )
          LogicError("z should be a column vector of length 3");
      if( d.Height() != 3 || d.Width() != 1 )
          LogicError("d should be a column vector of length 3");
    )
    const Real zero(0), one(1);
    const Real eps = limits::Epsilon<Real>();
    const Real safeMinToCube = limits::SafeMinToCube<Real>();
    const Real safeMinToCubeInv = one / safeMinToCube;
    const Real safeMinToRootCube = safeMinToCube*safeMinToCube;
    const Real safeMinToRootCubeInv = safeMinToCubeInv*safeMinToCubeInv;
    CubicSecularInfo info;

    Real rootLowerBound = ( rightRoot ? d(1) : d(0) );
    Real rootUpperBound = ( rightRoot ? d(2) : d(1) );

    if( originEval < zero )
    {
        // The root is to the right of zero
        rootLowerBound = zero;
    }
    else
    {
        // The root is to the left of zero
        rootUpperBound = zero;
    }

    Real rootEst = zero;
    if( initialize )
    {
        // Form the relevant quadratic equation
        // TODO(poulson): Document with the derivation
        Real a, bNeg, c;
        if( rightRoot )
        {
            a = rho + z(0) / ((d(0)-d(1))-(d(2)-d(1))/2);
            bNeg = a*(d(1)+d(2)) + z(1) + z(2);
            c = a*d(1)*d(2) + z(1)*d(2) + z(2)*d(1);
        }
        else
        {
            a = rho + z(2) / ((d(2)-d(1))-(d(0)-d(1))/2);
            bNeg = a*(d(0)+d(1)) + z(0) + z(1);
            c = a*d(0)*d(1) + z(0)*d(1) + z(1)*d(0);
        }

        // Normalize the coefficients of the quadratic equation
        const Real maxAbs = Max( Abs(a), Max( Abs(bNeg), Abs(c) ) );
        a /= maxAbs;
        bNeg /= maxAbs;
        c /= maxAbs;

        rootEst = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix );

        if( rootEst < rootLowerBound || rootEst > rootUpperBound )
        {
            // We have a nonsensical result, so default to the average
            rootEst = (rootLowerBound+rootUpperBound) / 2;
        }
        if( d(0) == rootEst || d(1) == rootEst || d(2) == rootEst )
        {
            rootEst = zero;
        }
        else
        {
            // Carefully evaluate f(rootEst) by perturbing originEval=f(0)
            const Real secular = originEval +
              rootEst*z(0)/(d(0)*(d(0)-rootEst)) +
              rootEst*z(1)/(d(1)*(d(1)-rootEst)) +
              rootEst*z(2)/(d(2)*(d(2)-rootEst));

            if( secular <= zero )
            {
                // The root is right of rootEst
                rootLowerBound = rootEst; 
            }
            else
            {
                // The root is left of rootEst
                rootUpperBound = rootEst;
            }

            if( Abs(originEval) <= Abs(secular) )
            {
                // The origin is better than the current estimate
                rootEst = zero;
            }
        }
    }

    bool rescale = false;
    Real scaleInv(1);
    Matrix<Real> zScaled(3,1), dScaled(3,1);
    Real maxDenomAbs;
    if( rightRoot )
    {
        maxDenomAbs = Min( Abs(d(1)-rootEst), Abs(d(2)-rootEst) );
    }
    else
    {
        maxDenomAbs = Min( Abs(d(0)-rootEst), Abs(d(1)-rootEst) );
    }
    if( maxDenomAbs <= safeMinToCube ) 
    {
        rescale = true; 
        Real scale;
        if( maxDenomAbs <= safeMinToRootCube )
        {
            scale = safeMinToRootCubeInv;
            scaleInv = safeMinToRootCube;
        }
        else
        {
            scale = safeMinToCubeInv;
            scaleInv = safeMinToCube;
        }

        for( Int i=0; i<3; ++i )
        {
            zScaled(i) = z(i)*scale;
            dScaled(i) = d(i)*scale;
        }
        rootEst *= scale;
        rootLowerBound *= scale;
        rootUpperBound *= scale;
    }
    else
    {
        zScaled = z;
        dScaled = d;
    }

    // Compute a relative correction to the f(0) to perturb to f(rootEst).
    // Also compute f'(rootEst) and f''(rootEst)/2.
    Real secularRelCorrection = zero;
    Real secularDeriv = zero;
    Real secularSecondDerivHalf = zero;
    for( Int i=0; i<3; ++i )
    {
        // The i'th term is
        //
        //   f_i(x) = z(i) / (d(i)-x),
        //
        // so its first derivative is
        //
        //   f'_i(x) = z(i) / (d(i)-x)^2,
        //
        // and its second derivative is
        //
        //   f''_i(x) = 2 z(i) / (d(i)-x)^3.
        //
        const Real temp = one / (dScaled(i)-rootEst);
        const Real temp1 = zScaled(i)*temp;
        const Real temp2 = temp1*temp;
        const Real temp3 = temp2*temp;
        secularRelCorrection += temp1 / dScaled(i);
        secularDeriv += temp2;
        secularSecondDerivHalf += temp3;
    }
    Real secular = originEval + rootEst*secularRelCorrection;
    ++info.numIterations;
    if( Abs(secular) == zero )
    {
        if( rescale )
            rootEst *= scaleInv;
        root = rootEst;
        return info;
    }
    if( secular <= zero )
    {
        // The root is right of our current estimate
        rootLowerBound = rootEst;
    }
    else
    {
        // The root is left of our current estimate
        rootUpperBound = rootEst;
    }

    // Begin Borges/Gragg/Thornton/Warner scheme
    while( true )
    {
        if( info.numIterations >= ctrl.maxIterations )
        {
            info.converged = false;
            break;
        }

        const Real leftDenom =
          ( rightRoot ? dScaled(1) : dScaled(0) ) - rootEst;
        const Real rightDenom =
          ( rightRoot ? dScaled(2) : dScaled(1) ) - rootEst;

        Real a = secular - (leftDenom+rightDenom)*secularDeriv +
          leftDenom*rightDenom*secularSecondDerivHalf;

        Real bNeg = (leftDenom+rightDenom)*secular -
          leftDenom*rightDenom*secularDeriv;

        Real c = leftDenom*rightDenom*secular;

        // Normalize the coefficients of the quadratic equation
        const Real maxAbs = Max( Abs(a), Max( Abs(bNeg), Abs(c) ) );
        a /= maxAbs; 
        bNeg /= maxAbs;
        c /= maxAbs;

        Real eta = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix );
        if( secular*eta >= zero )
        {
            // The current update does not move in the right direction, so fall
            // back to a small Newton step (as the derivative is likely large).
            eta = -secular / secularDeriv;
        }

        rootEst += eta;
        if( rootEst < rootLowerBound || rootEst > rootUpperBound )
        {
            // We have a nonsensical answer, so restart in the center
            rootEst = (rootLowerBound+rootUpperBound) / 2;
        }
        ++info.numIterations;

        bool converged = false;
        for( Int i=0; i<3; ++i )
        {
            if( dScaled(i)-rootEst == zero )
            {
                converged = true;
                break;
            }
        }
        if( converged )
            break;

        secularRelCorrection = zero;
        secularDeriv = zero;
        secularSecondDerivHalf = zero;
        Real relErrorBound = zero;
        for( Int i=0; i<3; ++i )
        {
            const Real temp = one / (dScaled(i)-rootEst);
            const Real temp1 = zScaled(i)*temp;
            const Real temp2 = temp1*temp;
            const Real temp3 = temp2*temp;
            const Real temp4 = temp1 / dScaled(i);
            secularRelCorrection += temp4;
            relErrorBound += Abs(temp4);
            secularDeriv += temp2;
            secularSecondDerivHalf += temp3;
        }
        secular = originEval + rootEst*secularRelCorrection;
        relErrorBound = 8*(Abs(originEval)+Abs(rootEst)*relErrorBound) +
          Abs(rootEst)*secularDeriv;

        if( Abs(secular) <= eps*relErrorBound )
        {
            // We have converged
            break;
        }
        if( secular <= zero ) 
            rootLowerBound = rootEst;
        else
            rootUpperBound = rootEst;
    }

    if( rescale )
        rootEst *= scaleInv;
    root = rootEst;
    return info;
}

#define PROTO(Real) \
  template CubicSecularInfo CubicSecular \
  ( bool initialize, \
    bool rightRoot, \
    const Real& rho, \
    const Matrix<Real>& z, \
    const Matrix<Real>& d, \
    const Real& originEval, \
          Real& root, \
    const CubicSecularCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT

#include <El/macros/Instantiate.h>

} // namespace El
