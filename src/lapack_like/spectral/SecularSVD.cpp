/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./SecularSVD/TwoByTwo.hpp"

namespace El {

template<typename Real,typename=EnableIf<IsReal<Real>>>
struct CubicSecularSingularValueInfo
{
    Real root;
    Int numIterations = 0;
    bool converged = true;
};

namespace secular_svd {

template<typename Real>
struct State
{
    Int origin;
    bool originOnLeft;

    Matrix<Real> dPlusShift;
    Matrix<Real> dMinusShift;

    Real diagSqDiff;

    Real sigmaEst; 
    Real sigmaRelEst;
    Real sigmaRelLowerBound;
    Real sigmaRelUpperBound;

    // The secular equation and its derivative evaluated at the current
    // eigenvalue estimate
    Real secular, secularOld;
    Real secularDeriv;

    // The secular equation (and its derivative) with the origin term removed
    Real secularMinus;
    Real secularMinusDeriv;

    // The sum of the secular equation terms left of the origin (not inclusive)
    Real psiMinus;
    Real psiMinusDeriv;

    // The sum of the secular equation terms right of the origin (not inclusive)
    Real phi;
    Real phiDeriv;

    Real relErrorBound;

    bool useThreePoles;
    bool geometricAverage;
    bool alternateStrategy;
};

template<typename Real>
struct LastState
{
    Matrix<Real> dPlusShift;
    Matrix<Real> dMinusShift;

    Real diagSqDiff;

    Real sigmaEst; 
    Real sigmaRelEst;

    // The secular equation and its derivative evaluated at the current
    // eigenvalue estimate
    Real secular, secularOld;
    Real secularDeriv;

    // The secular equation (and its derivative) with the origin term removed
    Real secularMinus;
    Real secularMinusDeriv;

    // The sum of the secular equation terms left of the origin (not inclusive)
    Real psiMinus;
    Real psiMinusDeriv;

    // The origin term of psi and its derivative
    Real psiOrigin;
    Real psiOriginDeriv;

    Real relErrorBound;
};

// Solve for an inner root of the secular equation
//
//   f(x) = rho + z(0) / (d(0)-x) + z(1) / (d(1)-x) + z(2) / (d(2)-x),
//
// where each numerator is positive and d(0) < d(1) < d(2).
//
// Just as in LAPACK's {s,d}laed6 [CITATION], we require that the user pass in
// an accurate evaluation of f(0).
//
template<typename Real,typename=EnableIf<IsReal<Real>>>
CubicSecularSingularValueInfo<Real>
CubicSecular
( bool rightRoot,
  const Real& rho,
  const Matrix<Real>& z,
  const Matrix<Real>& d,
  const Real& originEval,
  bool initialize,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    CubicSecularSingularValueInfo<Real> info;

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
    Real scaleInv;
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
        Real scale, scaleInv;
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
        info.root = rootEst;
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
        if( info.numIterations >= ctrl.maxCubicIterations )
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
    info.root = rootEst;
    return info;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void EvaluateSecular
( const Real& rho,
  const Matrix<Real>& u,
        State<Real>& state )
{
    DEBUG_CSE
    const Real one(1);
    const Int n = u.Height();
    const Int origin = state.origin;
    const Real rhoInv = one / rho;
    Real temp;

    // Compute psi_m and psi'_m, where m is the origin index, and an
    // approximation of the error
    // (see Ren-Cang Li, "Solving Secular Equations Stably and Efficiently",
    // LAPACK Working Note 89, 1993 [CITATION], as well as LAPACK's 
    // {s,d}lasd4 [CITATION]). The loop direction is chosen to heuristically
    // sum from small to large components.
    state.psiMinus = 0;
    state.psiMinusDeriv = 0;
    state.relErrorBound = 0;
    for( Int j=0; j<origin; ++j )
    {
        // Compute the j'th term divided by u(j)
        temp = u(j) / (state.dPlusShift(j)*state.dMinusShift(j));
        state.psiMinus += u(j)*temp; // This should be a negative contribution
        state.psiMinusDeriv += temp*temp;
        state.relErrorBound += state.psiMinus;
    }
    state.relErrorBound = Abs(state.relErrorBound); // This should be negation

    // Compute phi, its derivative, and accumulate an approximation of the
    // error in both psi_{m-1} and phi_m, where m is the origin index. The loop
    // direction is chosen to heuristically sum from small to large components.
    state.phi = 0;
    state.phiDeriv = 0;
    for( Int j=n-1; j>origin; --j )
    {
        // Compute the j'th term divided by u(j)
        temp = u(j) / (state.dPlusShift(j)*state.dMinusShift(j));
        state.phi += u(j)*temp;
        state.phiDeriv += temp*temp;
        state.relErrorBound += state.phi;
    }

    // Compute the secular function with the origin term removed
    // (and its derivative)
    state.secularMinus = rhoInv + state.psiMinus + state.phi;
    state.secularMinusDeriv = state.psiMinusDeriv + state.phiDeriv;

    // Cf. LAPACK's {s,d}lasd4 [CITATION] for this computational strategy
    temp = u(origin) / (state.dPlusShift(origin)*state.dMinusShift(origin));
    state.secularDeriv = (state.psiMinusDeriv + temp*temp) + state.phiDeriv;
    temp *= u(origin);
    state.secular = state.secularMinus + temp;
    state.relErrorBound +=
      8*(state.phi-state.psiMinus) + 2*rhoInv + 3*Abs(temp);
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void EvaluateSecularLast
( const Real& rho,
  const Matrix<Real>& u,
        LastState<Real>& state )
{
    DEBUG_CSE
    const Real one(1);
    const Int n = u.Height();
    const Int origin = n-1;
    const Real rhoInv = one / rho;
    Real temp;

    state.psiMinus = 0;
    state.psiMinusDeriv = 0;
    state.relErrorBound = 0;
    for( Int j=0; j<origin; ++j )
    {
        // Compute the j'th term divided by u(j)
        temp = u(j) / (state.dPlusShift(j)*state.dMinusShift(j));
        state.psiMinus += u(j)*temp;
        state.psiMinusDeriv += temp*temp;
        state.relErrorBound += state.psiMinus; // This should be negative
    }
    state.relErrorBound = Abs(state.relErrorBound); // This should be negation

    // Compute the origin term divided by u(origin)
    temp = u(origin) / (state.dPlusShift(origin)*state.dMinusShift(origin));
    state.psiOrigin = u(origin)*temp;
    state.psiOriginDeriv = temp*temp;
    state.relErrorBound -= state.psiOrigin;

    // Compute the secular function with the origin term removed
    // (and its derivative)
    const Real psi = state.psiMinus + state.psiOrigin;
    state.secular = rhoInv + psi;

    state.relErrorBound += rhoInv - 8*psi;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularInitialGuess
( const Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        State<Real>& state,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Real zero(0), one(1);
    const Int k = whichSingularValue;
    const Int n = d.Height();
    const Real diagSqDiffHalf = state.diagSqDiff / 2;
    if( k >= n-1 )
        LogicError("Assumption of inner singular value broken");

    // We can test whether the k'th eigenvalue is in the first or second
    // half of (d(k)^2, d(k+1)^2) by testing the sign of the secular
    // equation evaluated at the center point,
    //
    //   center = (d(k)^2 + d(k+1)^2)/2,
    //
    // and its square-root,
    //
    //   centerRoot = sqrt(center).
    //
    const Real center = (d(k)*d(k) + d(k+1)*d(k+1)) / 2; 
    const Real centerRoot = Sqrt(center);

    // Follow LAPACK's {s,d}lasd4's [CITATION] lead in carefully computing
    //
    //   shiftedCenterRoot = centerRoot - d(k).
    //
    const Real shiftedCenterRoot = diagSqDiffHalf / (d(k) + centerRoot);

    // Again use LAPACK's {s,d}lasd4's [CITATION] suggestion that the
    // diagonal entries plus-or-minus the center roots be computed in a
    // safe (but somewhat obscured) manner.
    for( Int j=0; j<n; ++j ) 
    {
        state.dPlusShift(j) = (d(j) + d(k)) + shiftedCenterRoot;
        state.dMinusShift(j) = (d(j) - d(k)) - shiftedCenterRoot;
    }

    // Given the partition of the secular equation as
    //
    //   f(x) = (1/rho) + sum_{j=0}^{n-1} u(j)^2 / (d(j)^2 - x)
    //        = (1/rho) + psi_k(x) + phi_k(x),
    //
    // where
    //
    //   psi_k(x) = sum_{j=0}^k       u(j)^2 / (d(j)^2 - x),
    //   phi_k(x) = sum_{j=k+1}^{n-1} u(j)^2 / (d(j)^2 - x),
    //
    // we follow the suggestion surrounding Eq's (40) and (41) of
    //
    //   Ren-Cang Li, "Solving Secular Equations Stably and Efficiently",
    //   LAPACK Working Note 89, 1993 [CITATION]
    //
    // and split off the last term of psi_k(x) and first term of phi_k(x)
    // so that the initial condition can be subsequently quickly computed.
    // We denote psi_k and phi_k evaluated at the center point with said 
    // terms removed as 'psiMinus' and 'phiMinus', respectively.
    state.psiMinus = zero;
    for( Int j=0; j<k; ++j )
        state.psiMinus +=
          u(j)*u(j) / (state.dPlusShift(j)*state.dMinusShift(j));
    Real phiMinus = zero; 
    for( Int j=k+2; j<n; ++j )
        phiMinus +=
          u(j)*u(j) / (state.dPlusShift(j)*state.dMinusShift(j));

    // The secular equation can now be expressed as
    //
    //   f(center) = (1/rho) + psiMinus + phiMinus +
    //       u(k)^2 / (d(k)^2 - center) + u(k+1)^2 / (d(k+1)^2 - center),
    //
    // where the top row will turn out to be the 'a' coefficient of a 
    // quadratic equation a x^2 + b x + c = 0 that we will solve to compute
    // the initial guess for the sought-after eigenvalue. But recall that
    // we should carefully compute (d(j)^2 - center).
    const Real a = one/rho + state.psiMinus + phiMinus;
    const Real secularCenter = a +
      u(k)*u(k) / (state.dPlusShift(k)*state.dMinusShift(k)) +
      u(k+1)*u(k+1) / (state.dPlusShift(k+1)*state.dMinusShift(k+1));

    state.geometricAverage = false; // TODO(poulson): Documentation?
    if( secularCenter >= zero )
    {
        // The eigenvalue lives in (d(k)^2, center), so we solve the 
        // quadratic equation
        //
        //   a + u(k)^2 / (d(k)^2 - x) + u(k+1)^2 / (d(k+1)^2 - x) = 0,
        //
        // directly for tau = x - d(k)^2, which yields the quadratic
        // equation a tau^2 + b tau + c = 0, with
        //
        //   b = -a gap - (u(k)^2 + u(k+1)^2),
        //   c = u(k)^2 gap,
        //
        // with gap = d(k+1)^2 - d(k)^2.
        state.originOnLeft = true;
        state.origin = k;
        state.sigmaRelLowerBound = zero; 
        state.sigmaRelUpperBound = shiftedCenterRoot;
        const Real bNeg = a*state.diagSqDiff + u(k)*u(k) + u(k+1)*u(k+1);
        const Real c = u(k)*u(k)*state.diagSqDiff;

        // Compute tau = sigmaEst^2 - d(k)^2
        Real tau = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix ); 
 
        state.sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue( tau, d(state.origin) );

        // TODO(poulson): Document this LAPACK heuristic
        const Real eps = limits::Epsilon<Real>();
        const Real sqrtEps = Sqrt(eps);
        if( (d(k) <= sqrtEps*d(k+1)) &&
            (Abs(u(k)) <= sqrtEps) &&
            (d(k) > Real(0)) )
        {
            state.geometricAverage = true;
            state.sigmaRelEst = Min( 10*d(k), state.sigmaRelUpperBound );
        }
    }
    else
    {
        // The eigenvalue lives in [center, d(k+1)^2), so we solve the
        // quadratic equation
        //
        //   a + u(k)^2 / (d(k)^2 - x) + u(k+1)^2 / (d(k+1)^2 - x) = 0,
        //
        // directly for tau = x - d(k+1)^2, which yields the quadratic
        // equation a tau^2 + b tau + c = 0, with
        //
        //   b = a gap - (u(k)^2 + u(k+1)^2),
        //   c = -u(k)^2 gap,
        //
        // with gap = d(k+1)^2 - d(k)^2.
        state.originOnLeft = false;
        state.origin = k+1;

        // Again follow LAPACK's {s,d}lasd4's [CITATION] lead in carefully
        // computing centerRoot - d(k+1).
        state.sigmaRelLowerBound =
          -diagSqDiffHalf / (d(state.origin) + centerRoot);
        state.sigmaRelUpperBound = 0; 
        const Real bNeg = -a*state.diagSqDiff + u(k)*u(k) + u(k+1)*u(k+1);
        const Real c = -u(k+1)*u(k+1)*state.diagSqDiff;

        const Real tau = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix );

        state.sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue( tau, d(state.origin) );
    }

    state.sigmaEst = state.sigmaRelEst + d(state.origin);
    for( Int j=0; j<n; ++j ) 
    {
        state.dPlusShift(j) = (d(j) + d(state.origin)) + state.sigmaRelEst;
        state.dMinusShift(j) = (d(j) - d(state.origin)) - state.sigmaRelEst;
    }

    EvaluateSecular( rho, u, state );
}

// For seeking the last root of the secular equation in
//
//   (d(n-1)^2, d(n-1)^2+rho).
//
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularInitialGuessLast
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        LastState<Real>& state,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int n = d.Height();
    const Real zero(0), one(1);
    const Real rhoInv = one / rho;

    const Int origin = n-1;

    // Since the largest eigenvalue must live within
    //
    //   d(n-1)^2 + (0, rho ||u||_2^2) = d(n-1)^2 + (0, rho),
    //
    // we use the center as our initial guess. Therefore, put
    // 
    //   tau = (rho ||u||_2^2) / 2 = rho / 2 = sigmaEst^2 - d(n-1)^2.
    //
    Real tau = rho / 2;

    state.sigmaRelEst =
      RelativeEigenvalueToRelativeSingularValue( tau, d(origin) );

    state.sigmaEst = state.sigmaRelEst + d(origin);
    for( Int j=0; j<n; ++j ) 
    {
        state.dPlusShift(j) = (d(j) + d(origin)) + state.sigmaRelEst;
        state.dMinusShift(j) = (d(j) - d(origin)) - state.sigmaRelEst;
    }

    // Form psi_{origin-2}
    Real psiDoubleMinus = zero;
    psiDoubleMinus = zero;
    for( Int j=0; j<origin-1; ++j )
        psiDoubleMinus +=
          u(j)*u(j) / (state.dPlusShift(j)*state.dMinusShift(j));

    state.psiMinus = psiDoubleMinus +
      u(origin-1)*u(origin-1) /
      (state.dPlusShift(origin-1)*state.dMinusShift(origin-1));

    state.psiOrigin = u(origin)*u(origin) /
      (state.dPlusShift(origin)*state.dMinusShift(origin));

    const Real a = rhoInv + psiDoubleMinus;
    state.secularMinus = rhoInv + state.psiMinus;
    state.secular = state.secularMinus + state.psiOrigin;

    const Real diagSum = d(origin) + d(origin-1);
    const Real diagDiff = d(origin) - d(origin-1);
    state.diagSqDiff = diagSum*diagDiff;
    if( state.secular <= zero )
    {
        const Real singValUpperBound = Sqrt(d(origin)*d(origin) + rho);

        // TODO(poulson): Document this branch from LAPACK's {s,d}lasd4
        // [CITATION]
        const Real temp =
          u(origin-1)*u(origin-1) /
           ((d(origin-1) + singValUpperBound)*
            (d(origin) - d(origin-1) + rho/(d(origin) + singValUpperBound))) +
          u(origin)*u(origin) / rho;

        if( a <= temp )
        {
            tau = rho;
        }
        else
        {
            const Real bNeg = -a*state.diagSqDiff + u(origin-1)*u(origin-1) +
              u(origin)*u(origin);
            const Real c = -u(origin)*u(origin)*state.diagSqDiff;
            tau = SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
            state.sigmaRelEst =
              RelativeEigenvalueToRelativeSingularValue( tau, d(origin) );
        }
    }
    else
    {
        // We will approximate the secular equation with
        //
        //  f_{n-1,n-2}(y) + u(n-2)^2 / (d(n-2)^2 - x) +
        //                   u(n-1)^2 / (d(n-1)^2 - x),
        //
        // where f_{n-1,n-2}(y) is the secular equation, with the (n-1)'th and 
        // (n-2)'th terms removed, evaluated at the current estimate. We solve
        // for a root of this equation in terms of tau = x - d(n-1)^2.
        // In particular, we pick the '-' branch of the quadratic equation.

        const Real bNeg = -a*state.diagSqDiff + u(origin-1)*u(origin-1) +
          u(origin)*u(origin);
        const Real c = -u(origin)*u(origin)*state.diagSqDiff;
        tau = SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
        state.sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue( tau, d(origin) );
    }

    state.sigmaEst = state.sigmaRelEst + d(origin);
    for( Int j=0; j<n; ++j ) 
    {
        state.dPlusShift(j) = (d(j) + d(origin)) + state.sigmaRelEst;
        state.dMinusShift(j) = (d(j) - d(origin)) - state.sigmaRelEst;
    }

    EvaluateSecularLast( rho, u, state );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularUpdate
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        State<Real>& state,
  bool initialize,
  const SecularSingularValueCtrl<Real>& ctrl,
  SecularSingularValueInfo<Real>& info )
{
    DEBUG_CSE
    const Real zero(0);
    const Int n = d.Height();
    const Int k = whichSingularValue;
    const Int origin = state.origin;
    Real eta;

    if( state.useThreePoles )
    {
        // Use the "Hybrid Scheme" described in subsection 3.4 of LAWN 89 and
        // implemented within {s,d}lasd4.
 
        // Carefully compute
        //
        //   leftGap = d(origin-1)^2 - sigmaEst^2, and
        //   rightGap = d(origin+1)^2 - sigmaEst^2.
        //
        const Real leftGap =
          state.dPlusShift(origin-1)*state.dMinusShift(origin-1);
        const Real rightGap =
          state.dPlusShift(origin+1)*state.dMinusShift(origin+1);

        // Carefully compute doubleGap = d(origin+1)^2 - d(origin-1)^2
        const Real doubleGap =
          (d(origin+1)-d(origin-1))*(d(origin+1)+d(origin-1));

        Real a;
        Matrix<Real> zCubic(3,1);
        if( state.alternateStrategy )
        {
            a = state.secularMinus - leftGap*state.psiMinusDeriv -
              rightGap*state.phiDeriv;
            zCubic(0) = leftGap*leftGap*state.psiMinusDeriv;
            zCubic(1) = u(origin)*u(origin);
            zCubic(2) = rightGap*rightGap*state.phiDeriv;
        }
        else
        {
            if( state.originOnLeft )
            {
                // Since the shift origin, m, is k, We will interpolate the
                // secular equation as
                //
                //  Q(x; a, s, S) = a + u(m-1)^2 / (d(m-1)^2 - x) + 
                //                      u(m  )^2 / (d(m  )^2 - x) + 
                //                      S        / (d(m+1)^2 - x),
                //
                // which can be effected by applying the Middle Way, the Fixed
                // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme
                //  to
                //
                //   f_m(y) = f(y) - u(m)^2 / (d(m)^2 - x).
                //
                // We balance between the Fixed Weight Method and the Middle
                // Way; their formulae for computing a are identical. In both
                // cases,
                // 
                //  a = f_m - rightGap*f'_m +
                //      (u(m-1)/leftGap)^2*(d(m+1)^2-d(m-1)^2).
                //
                // The computation of S primarily follows the Fixed Weight
                // Method but enforces the non-negativity of psi_{k-2}' when
                // computing
                //
                //   f'_{m-1,m}(y) = psi'_{m-2}(y) + phi'_m(y).
                //
                // The two subscripts on f denote the m-1 and m terms of the
                // secular equation being left out.
                //
                const Real leftRatio = u(origin-1) / leftGap;
                const Real leftDerivTerm = leftRatio*leftRatio;
                const Real a =
                  (state.secularMinus-rightGap*state.secularMinusDeriv) +
                  leftDerivTerm*doubleGap;

                // We could either flip or clip here, but LAPACK tends to flip
                // and clips here, so we will follow suit so that requesting
                // flips leads to mirroring LAPACK.
                const Real psiDoubleMinusDeriv =
                  Max( state.psiMinusDeriv-leftDerivTerm, zero );
    
                zCubic(0) = u(origin-1)*u(origin-1);
                zCubic(1) = u(origin)*u(origin);
                zCubic(2) =
                  rightGap*rightGap*(psiDoubleMinusDeriv+state.phiDeriv);
            }
            else
            {
                // Since the shift origin, m, is k+1, we will interpolate the
                // secular equation as
                //
                //  Q(x; a, s, S) = a + s        / (d(m-1)^2 - x) +
                //                      u(m)^2   / (d(m  )^2 - x) +
                //                      u(m+1)^2 / (d(m+1)^2 - x),
                //
                // which can be effected by applying the Middle Way, the Fixed
                // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme
                //  to
                //
                //   f_m(y) = f(y) - u(m)^2 / (d(m)^2 - x).
                //
                // We balance between the Fixed Weight Method and the Middle
                // Way; their formulae for computing a are identical. In both
                // cases,
                // 
                //  a = f_m - leftGap*f'_m -
                //      (u(m+1)/rightGap)^2*(d(m+1)^2-d(m-1)^2),
                //
                // The computation of s primarily follows the Fixed Weight
                // Method but enforces the non-negativity of phi_{m+1}' when
                // computing
                //
                //   f'_{m,m+1}(y) = psi'_{m-1}(y) + phi'_{m+1}(y).
                //
                // The two subscripts on f denote the m and m+1 terms of the
                // secular equation being left out.
                //
                const Real rightRatio = u(origin+1) / rightGap;
                const Real mPlusDerivTerm = rightRatio*rightRatio;
                const Real a =
                  (state.secularMinus-leftGap*state.secularMinusDeriv) -
                  mPlusDerivTerm*doubleGap;

                // We could either flip or clip here, but LAPACK tends to flip
                // and clips here, so we will follow suit so that requesting
                // flips leads to mirroring LAPACK.
                const Real phiMinusDeriv =
                  Max( state.phiDeriv-mPlusDerivTerm, zero );
            
                zCubic(0) = leftGap*leftGap*(state.psiMinusDeriv+phiMinusDeriv);
                zCubic(1) = u(origin)*u(origin);
                zCubic(2) = u(origin+1)*u(origin+1);
            }
        }
        Matrix<Real> dCubic(3,1);
        dCubic(0) = leftGap;
        dCubic(1) = state.dPlusShift(origin)*state.dMinusShift(origin);
        dCubic(2) = rightGap;

        auto cubicInfo =
          CubicSecular
          ( state.originOnLeft, a, zCubic, dCubic, state.secular,
            initialize, ctrl );
        info.numCubicIterations += cubicInfo.numIterations;
        if( cubicInfo.converged )
        {
            eta = cubicInfo.root;
            if( ctrl.progress )
                Output
                ("Cubic converged in ",cubicInfo.numIterations,
                 " iter's with eta=",eta);
        }
        else
        {
            if( ctrl.progress )
                Output("Cubic did *not* converge");
            ++info.numCubicFailures;
            state.useThreePoles = false;
            const Real kGap = state.dPlusShift(k)*state.dMinusShift(k);
            const Real kp1Gap = state.dPlusShift(k+1)*state.dMinusShift(k+1);
            Real a;
            if( state.originOnLeft )
            {
                const Real temp = u(k) / kGap;
                a = state.secular - kp1Gap*state.secularDeriv +
                  state.diagSqDiff*(temp*temp);
            }
            else
            {
                const Real temp = u(k+1) / kp1Gap;
                a = state.secular - kGap*state.secularDeriv +
                  state.diagSqDiff*(temp*temp);
            }
            Real bNeg = (kGap+kp1Gap)*state.secular -
              kGap*kp1Gap*state.secularDeriv;
            const Real c = kGap*kp1Gap*state.secular;
            if( a == zero && bNeg == zero )
            {
                // a x^2 + b x + c = 0 has collapsed to c = 0, which is
                // nonsense. Follow LAPACK's {s,d}lasd4 [CITATION] in handling
                // this breakdown.
                if( state.originOnLeft )
                {
                    bNeg = u(k)*u(k) + kp1Gap*kp1Gap*state.secularMinusDeriv;
                }
                else
                {
                    bNeg = u(k+1)*u(k+1) + kGap*kGap*state.secularMinusDeriv; 
                }
            }
            eta = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix );
        }
    }
    else
    {
        // Use one step of the Fixed Weight Method, as described in subsection
        // 3.2 of LAPACK Working Note 89 [CITATION]. We solve the quadratic
        // equation implied by the osculatory interpolation of the secular
        // equation as
        //
        //   Q(x; a, s, S) = a + s / (d(k)^2 - x) + S / (d(k+1)^2 - x).
        //
        // As noted after Proposition 3 of LAWN 89, the iteration formula does
        // not depend upon either s or S!

        // Safely compute the distance from our current eigenvalue estimate to
        // the surrounding grid points, d(k)^2 and d(k+1)^2.
        const Real kGap = state.dPlusShift(k)*state.dMinusShift(k);
        const Real kp1Gap = state.dPlusShift(k+1)*state.dMinusShift(k+1);

        Real a;
        if( state.alternateStrategy )
        {
            // TODO(poulson): Document this strategy
            if( state.originOnLeft )
            {
                const Real temp = u(k) / kGap;
                a = state.secular - kGap*(state.psiMinusDeriv+temp*temp) -
                  kp1Gap*state.phiDeriv;
            }
            else
            {
                const Real temp = u(k+1) / kp1Gap;
                a = state.secular - kGap*state.psiMinusDeriv -
                  kp1Gap*(state.phiDeriv+temp*temp);
            }
        }
        else
        {
            if( state.originOnLeft )
            {
                // Apply Eq'ns (32) through (34) from LAWN 89 to set
                //
                //  a = f - kp1Gap f' + (u(k)^2/kGap^2) (d(k+1)^2 - d(k)^2),
                //
                // and implicitly set
                //
                //  s = u(k)^2,
                //  S = kp1Gap^2 (f' - u(k)^2/kGap^2),
                //
                // where y is the current estimate of the k'th eigenvalue and 
                // f and f' are the secular equation and its derivative
                // evaluated at y.
                const Real temp = u(k) / kGap;
                a = state.secular - kp1Gap*state.secularDeriv +
                  (temp*temp)*state.diagSqDiff;
            }
            else
            {
                // Apply Eq'ns (35) through (37) from LAWN 89 to set
                //
                //  a = f - kGap f' - (u(k+1)^2/kp1Gap^2) (d(k+1)^2 - d(k)^2),
                //
                // and implicitly set
                //
                //  s = kGap^2 (f' - u(k+1)^2/kp1Gap^2),
                //  S = u(k+1)^2.
                //
                const Real temp = u(k+1) / kp1Gap;
                a = state.secular - kGap*state.secularDeriv -
                  (temp*temp)*state.diagSqDiff;
            }
        }

        // See Proposition 3 from LAWN 89 [CITATION] for these formulae for the
        // coefficients of the quadratic equation a x^2 - bNeg x + c = 0.
        Real bNeg = (kGap+kp1Gap)*state.secular -
          kGap*kp1Gap*state.secularDeriv;
        Real c = kGap*kp1Gap*state.secular;

        if( a == zero && bNeg == zero )
        {
            // a eta^2 + b eta + c = 0 has collapsed to c = 0, which is
            // nonsensical to solve for eta. We therefore follow the lead of
            // LAPACK's {s,d}lasd4 [CITATION] and use a mysterious patch-up.
            //
            // TODO(poulson): Provide motivation for these formulae.
            if( state.alternateStrategy )
            {
                if( state.originOnLeft )
                {
                    const Real temp = u(k) / kGap;
                    bNeg = kGap*kGap*(state.psiMinusDeriv+temp*temp) +
                           kp1Gap*kp1Gap*state.phiDeriv;
                }
                else
                {
                    const Real temp = u(k+1) / kp1Gap;
                    bNeg = kGap*kGap*state.psiMinusDeriv +
                           kp1Gap*kp1Gap*(state.phiDeriv+temp*temp);
                }
            }
            else
            {
                if( state.originOnLeft )
                    bNeg = u(k)*u(k) + kp1Gap*kp1Gap*state.secularMinusDeriv;
                else
                    bNeg = u(k+1)*u(k+1) + kGap*kGap*state.secularMinusDeriv;
            }
        }
        eta = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix );
    }

    if( state.secular*eta >= zero )
    {
        if( ctrl.progress )
            Output("Falling back to Newton step");
        // The current update does not move in the right direction, so fall
        // back to a small Newton step (as the derivative is likely large).
        eta = -state.secular / state.secularDeriv;
    }

    // Convert eta from an eigenvalue update to singular value update
    eta = RelativeEigenvalueToRelativeSingularValue( eta, state.sigmaEst );

    // Since sigmaRelEst = sigmaEst - d(origin), the following sets
    // sigmaRelEstProp = sigmaEstProp - d(origin).
    const Real sigmaRelEstProp = state.sigmaRelEst + eta;
    if( sigmaRelEstProp > state.sigmaRelUpperBound ||
        sigmaRelEstProp < state.sigmaRelLowerBound )
    {
        if( ctrl.progress )
            Output("Stepped out of bounds");
        // We again follow LAPACK's {s,d}lasd4 for how to handle this breakdown

        if( state.secular < zero )
        {
            // Since the secular equation implies that the root is to the right
            // of the un-updated estimate, but our update pushed us out of
            // bounds, we default to the update pushing us to the center of the
            // feasible domain, [sigmaRelEst, sigmaRelUpperBound).
            eta = (state.sigmaRelUpperBound - state.sigmaRelEst) / 2;
        }
        else
        {
            // Since the secular equation implies that the root is to the left
            // of the un-updated estimate, but our update pushed us out of
            // bounds, we default to the update being in the center of the
            // feasible domain, (sigmaRelLowerBound,sigmaRelEst].
            eta = (state.sigmaRelLowerBound - state.sigmaRelEst) / 2;
        }

        if( state.geometricAverage )
        {
            if( state.secular < zero )
            {
                if( state.sigmaRelEst > zero )
                    eta = Sqrt(state.sigmaRelUpperBound*state.sigmaRelEst) -
                      state.sigmaRelEst;
            }
            else
            {
                if( state.sigmaRelLowerBound > zero )
                    eta = Sqrt(state.sigmaRelLowerBound*state.sigmaRelEst) -
                      state.sigmaRelEst;
            }
        }
    }

    state.sigmaEst += eta;
    state.sigmaRelEst += eta;
    state.secularOld = state.secular;
    for( Int j=0; j<n; ++j )
    {
        state.dPlusShift(j) += eta;
        state.dMinusShift(j) -= eta;
    }

    EvaluateSecular( rho, u, state );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularUpdateLast
( const Real& rho,
  const Matrix<Real>& u,
        LastState<Real>& state,
  bool initialize,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Real zero(0);
    const Int n = state.dPlusShift.Height();
    const Int origin = n-1;

    const Real kGap = state.dPlusShift(origin)*state.dMinusShift(origin);
    const Real km1Gap = state.dPlusShift(origin-1)*state.dMinusShift(origin-1);

    Real a = state.secular - km1Gap*state.psiMinusDeriv -
      kGap*state.psiOriginDeriv;

    Real eta;
    if( initialize )
    {
        // Force flipping always?
        if( ctrl.negativeFix == CLIP_NEGATIVES )
            a = Max( a, zero );
        else
            a = Abs( a );

        if( a == zero )
        {
            // TODO(poulson): Explain this special case
            eta = rho - state.sigmaEst*state.sigmaEst;
        }
        else
        {
            // We Approach From the Right with psi_{n-2}(x) osculatorally
            // interpolated as G(x; d(n-2)^2, r, s) and phi_{n-2}(x) exactly
            // represented as F(x; p, q). See the discussion surrounding Eq'n
            // (23) of LAWN 89 [CITATION], but keep in mind that we use the
            // definitions (a,b,c) corresponding to the standard quadratic 
            // equation a x^2 + b x + c = 0 rather than the notation of LAWN 89.
            const Real bNeg =
              (kGap+km1Gap)*state.secular -
              kGap*km1Gap*(state.psiMinusDeriv+state.psiOriginDeriv);
            const Real c = kGap*km1Gap*state.secular;
            eta = SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
        }
    }
    else
    {
        const Real bNeg =
          (kGap+km1Gap)*state.secular -
          kGap*km1Gap*(state.psiMinusDeriv+state.psiOriginDeriv);
        const Real c = kGap*km1Gap*state.secular;
        eta = SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
    }

    if( state.secular*eta >= zero )
    {
        // The current update does not move in the right direction, so fall back
        // to a small Newton step (as the derivative is likely large).
        if( ctrl.progress )
            Output("Falling back to Newton step");
        eta = -state.secular / (state.psiMinusDeriv+state.psiOriginDeriv);
    }

    if( initialize )
    {
        // eta = Min( eta, rho+kGap )
        if( eta-kGap > rho )
        {
            if( ctrl.progress )
                Output("Stepped out of bounds");
            eta = rho + kGap;
        }
    }
    else
    {
        if( eta-kGap < zero )
        {
            if( ctrl.progress )
                Output("Stepped out of bounds");
            eta /= 2;
        }
    }

    // Convert eta from an eigenvalue update to singular value update
    eta = RelativeEigenvalueToRelativeSingularValue( eta, state.sigmaEst );

    state.sigmaEst += eta;
    state.sigmaRelEst += eta;
    state.secularOld = state.secular;
    for( Int j=0; j<n; ++j )
    {
        state.dPlusShift(j) += eta;
        state.dMinusShift(j) -= eta;
    }
    EvaluateSecularLast( rho, u, state );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularSingularValueInfo<Real>
SecularInner
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        State<Real>& state,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    Real temp; // for frequent usage as a temporary product
    const Real zero(0);
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichSingularValue;
    const Int n = d.Height();
    DEBUG_ONLY(
      if( k == n-1 )
          LogicError("SecularInner meant for inner singular values");
      if( n <= 2 )
          LogicError("SecularInner meant for n > 2");
    )

    SecularSingularValueInfo<Real> info;
    state.dMinusShift.Resize(n,1);
    state.dPlusShift.Resize(n,1);

    const Real diagSum = d(k+1) + d(k);
    const Real diagDiff = d(k+1) - d(k);
    state.diagSqDiff = diagSum*diagDiff; // d(k+1)^2 - d(k)^2

    SecularInitialGuess( k, d, rho, u, state, ctrl );
    ++info.numIterations;

    // See the discussion in LAWN 89 [CITATION] on the usage of three poles
    state.useThreePoles =
      (state.originOnLeft ? state.secularMinus<zero : state.secularMinus>zero);
    if( state.origin == 0 || state.origin == n-1 )
        state.useThreePoles = false;

    // Check if we have already converged
    if( Abs(state.secular) <= eps*state.relErrorBound )
    {
        info.singularValue = state.sigmaEst;
        return info;
    }

    if( state.secular <= zero )
    {
        // Mark the sought eigenvalue as right of our current estimate
        state.sigmaRelLowerBound =
          Max( state.sigmaRelLowerBound, state.sigmaRelEst );
    }
    else
    {
        // Mark the sought eigenvalue as left of our current estimate
        state.sigmaRelUpperBound =
          Min( state.sigmaRelUpperBound, state.sigmaRelEst );
    }

    // Compute the first update to our estimate
    bool initialize = true;
    state.alternateStrategy = false;
    SecularUpdate( k, d, rho, u, state, initialize, ctrl, info );
    ++info.numIterations;
   
    // This strategy was described in Ren-Cang Li's LAWN 89
    if( state.originOnLeft )
    {
        if( -state.secular > Abs(state.secularOld)*ctrl.sufficientDecay )
        {
            state.alternateStrategy = true;        
            ++info.numAlternations;
        }
    }
    else
    {
        if( state.secular > Abs(state.secularOld)*ctrl.sufficientDecay )
        {
            state.alternateStrategy = true;
            ++info.numAlternations;
        }
    }

    initialize = false;
    while( true )
    {
        DEBUG_ONLY(
          if( !limits::IsFinite(state.sigmaEst) )
          {
              RuntimeError("Produced non-finite sigmaEst=",state.sigmaEst);
          }
        )
        if( Abs(state.secular) <= eps*state.relErrorBound ) 
        {
            // We have converged
            break;
        }
        if( info.numIterations >= ctrl.maxIterations )
        {
            RuntimeError
            ("Secular solver did not converge in ",ctrl.maxIterations,
             " iterations");
        }

        if( state.secular <= zero )
            state.sigmaRelLowerBound =
              Max( state.sigmaRelLowerBound, state.sigmaRelEst );
        else
            state.sigmaRelUpperBound =
              Min( state.sigmaRelUpperBound, state.sigmaRelEst );

        // Decide the next step
        SecularUpdate( k, d, rho, u, state, initialize, ctrl, info );
        ++info.numIterations;
   
        // This strategy was described in Ren-Cang Li's LAWN 89
        if( state.secular*state.secularOld > zero &&
            Abs(state.secular) > Abs(state.secularOld)*ctrl.sufficientDecay )
        {
            state.alternateStrategy = !state.alternateStrategy;
            ++info.numAlternations;
        }
    }
    info.singularValue = state.sigmaEst;
    return info;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularSingularValueInfo<Real>
SecularLast
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        LastState<Real>& state,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichSingularValue;
    const Int n = d.Height();
    DEBUG_ONLY(
      if( k != n-1 )
          LogicError("SecularLast meant for largest singular value");
      if( n <= 2 )
          LogicError("SecularLast meant for n > 2");
    )

    SecularSingularValueInfo<Real> info;
    state.dMinusShift.Resize(n,1);
    state.dPlusShift.Resize(n,1);

    SecularInitialGuessLast( d, rho, u, state, ctrl );
    ++info.numIterations;

    if( Abs(state.secular) <= eps*state.relErrorBound )
    {
        info.singularValue = state.sigmaEst;
        return info;
    }

    // Calculate the first update
    bool initialize = true;
    SecularUpdateLast( rho, u, state, initialize, ctrl );
    ++info.numIterations;

    initialize = false;
    while( true )
    {
        if( Abs(state.secular) <= eps*state.relErrorBound ) 
        {
            // We have converged
            break;
        }
        if( info.numIterations >= ctrl.maxIterations )
        {
            RuntimeError
            ("Secular solver did not converge in ",ctrl.maxIterations,
             " iterations");
        }

        // Decide the next step
        SecularUpdateLast( rho, u, state, initialize, ctrl );
        ++info.numIterations;
    }
    info.singularValue = state.sigmaEst;
    return info;
}

} // namespace secular_svd

// Compute a single singular value corresponding to the square-root of the 
// eigenvalue of the diagonal plus rank one matrix
//
//     diag(d)^2 + rho u u^T,
//
// where || u ||_2 = 1, with
//
//     0 <= d(0) < d(1) < ... < d(n-1)
//
// and rho > 0.
//
// This routine loosely corresponds to LAPACK's {s,d}lasd4 [CITATION].
//

template<typename Real,typename>
SecularSingularValueInfo<Real>
SecularSingularValue
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int k = whichSingularValue;
    const Int n = d.Height();
    DEBUG_ONLY(
      const Real zero(0);
      if( k < 0 || k >= n )
          LogicError("Invalid singular value request");
      if( d(0) < zero )
          LogicError("Assumption that d(0) >= 0 was broken");
      for( Int j=0; j<n-1; ++j )
          if( d(j) >= d(j+1) )
              LogicError("Assumption that d(j) < d(j+1) broken");
      if( rho <= zero )
          LogicError("Assumption that rho > 0 was broken");
      if( u.Height() != n )
          LogicError("u was not the correct height");
      // TODO(poulson): Check the assumption that || u ||_2 = 1
    )

    SecularSingularValueInfo<Real> info;
    if( n == 1 )
    {
        info.singularValue = Sqrt(d(0)*d(0) + rho*u(0)*u(0));
        return info;
    }
    else if( n == 2 )
    {
        info.singularValue =
          secular_svd::TwoByTwo( k, d(0), d(1), rho, u(0), u(1) );
        return info;
    }

    if( k < n-1 )
    {
        secular_svd::State<Real> state;
        info = secular_svd::SecularInner( k, d, rho, u, state, ctrl );
    }
    else
    {
        secular_svd::LastState<Real> state;
        info = secular_svd::SecularLast( k, d, rho, u, state, ctrl );
    }

    return info;
}

template<typename Real,typename>
SecularSingularValueInfo<Real>
SecularSingularValue
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        Matrix<Real>& dMinusShift,
        Matrix<Real>& dPlusShift,
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int k = whichSingularValue;
    const Int n = d.Height();
    DEBUG_ONLY(
      const Real zero(0);
      if( k < 0 || k >= n )
          LogicError("Invalid singular value request");
      if( d(0) < zero )
          LogicError("Assumption that d(0) >= 0 was broken");
      for( Int j=0; j<n-1; ++j )
          if( d(j) >= d(j+1) )
              LogicError("Assumption that d(j) < d(j+1) broken");
      if( rho <= zero )
          LogicError("Assumption that rho > 0 was broken");
      if( u.Height() != n )
          LogicError("u was not the correct height");
      // TODO(poulson): Check the assumption that || u ||_2 = 1
    )

    SecularSingularValueInfo<Real> info;
    dMinusShift.Resize(n,1);
    dPlusShift.Resize(n,1);
    if( n == 1 )
    {
        info.singularValue = Sqrt(d(0)*d(0) + rho*u(0)*u(0));
        // TODO(poulson): Make this computation more accurate for completeness
        dMinusShift(0) = d(0) - info.singularValue;
        dPlusShift(0) = d(0) + info.singularValue;
        return info;
    }
    else if( n == 2 )
    {
        info.singularValue =
          secular_svd::TwoByTwo
          ( k, d(0), d(1), rho, u(0), u(1),
            dMinusShift(0), dMinusShift(1),
            dPlusShift(0), dPlusShift(1) );
        return info;
    }

    if( k < n-1 )
    {
        secular_svd::State<Real> state;
        info = secular_svd::SecularInner( k, d, rho, u, state, ctrl );
        dPlusShift = state.dPlusShift;
        dMinusShift = state.dMinusShift;
    }
    else
    {
        secular_svd::LastState<Real> state;
        info = secular_svd::SecularLast( k, d, rho, u, state, ctrl );
        dPlusShift = state.dPlusShift;
        dMinusShift = state.dMinusShift;
    }

    return info;
}

#define PROTO(Real) \
  template SecularSingularValueInfo<Real> \
  SecularSingularValue \
  ( Int whichSingularValue, \
    const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& u, \
    const SecularSingularValueCtrl<Real>& ctrl ); \
  template SecularSingularValueInfo<Real> \
  SecularSingularValue \
  ( Int whichSingularValue, \
    const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& u, \
          Matrix<Real>& dMinusShift, \
          Matrix<Real>& dPlusShift, \
    const SecularSingularValueCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT

#include <El/macros/Instantiate.h>

} // namespace El
