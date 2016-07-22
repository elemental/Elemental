/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_SVD_DC_SECULAR_HPP
#define EL_BIDIAG_SVD_DC_SECULAR_HPP

#include "./TwoByTwoSecular.hpp"

namespace El {
namespace bidiag_svd {
namespace dc {

template<typename Real,typename=EnableIf<IsReal<Real>>>
struct CubicSecularInfo
{
    Real root;
    Int numIterations = 0;
    bool converged = true;
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
CubicSecularInfo<Real> CubicSecular
( bool rightRoot,
  const Real& rho,
  const Matrix<Real>& z,
  const Matrix<Real>& d,
  const Real& originEval,
  bool initialize )
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
    CubicSecularInfo<Real> info;

    // Cf. LAPACK's {s,d}laed6 for this choice [CITATION]
    const Int maxIterations = 40;

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

        rootEst = SolveQuadratic( a, bNeg, c );

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
        if( info.numIterations >= maxIterations )
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

        Real eta = SolveQuadratic( a, bNeg, c );
        if( secular*eta >= zero )
        {
            // Fall back to a Newton step
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
void SecularInitialGuess
( const Int k,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
  const Real& diagSqDiff,
  Matrix<Real>& dPlusShift,
  Matrix<Real>& dMinusShift,
  Real& sigmaEst,
  Real& sigmaRelEst,
  Real& sigmaRelLowerBound,
  Real& sigmaRelUpperBound,
  Int& origin,
  bool& originOnLeft,
  bool& geometricAverage )
{
    DEBUG_CSE
    const Real zero(0), one(1);
    const Int n = d.Height();
    const Real diagSqDiffHalf = diagSqDiff / 2;
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
        dPlusShift(j) = (d(j) + d(k)) + shiftedCenterRoot;
        dMinusShift(j) = (d(j) - d(k)) - shiftedCenterRoot;
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
    Real psiMinus = zero;
    for( Int j=0; j<k; ++j )
        psiMinus += u(j)*u(j) / (dPlusShift(j)*dMinusShift(j));
    Real phiMinus = zero; 
    for( Int j=k+2; j<n; ++j )
        phiMinus += u(j)*u(j) / (dPlusShift(j)*dMinusShift(j));

    // The secular secular equation can now be expressed as
    //
    //   f(center) = (1/rho) + psiMinus + phiMinus +
    //       u(k)^2 / (d(k)^2 - center) + u(k+1)^2 / (d(k+1)^2 - center),
    //
    // where the top row will turn out to be the 'a' coefficient of a 
    // quadratic equation a x^2 + b x + c = 0 that we will solve to compute
    // the initial guess for the sought-after eigenvalue. But recall that
    // we should carefully compute (d(j)^2 - center).
    const Real a = one/rho + psiMinus + phiMinus;
    const Real secularCenter = a +
      u(k)*u(k) / (dPlusShift(k)*dMinusShift(k)) +
      u(k+1)*u(k+1) / (dPlusShift(k+1)*dMinusShift(k+1));

    geometricAverage = false; // TODO(poulson): Documentation?
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
        originOnLeft = true;
        origin = k;
        sigmaRelLowerBound = zero; 
        sigmaRelUpperBound = shiftedCenterRoot;
        const Real bNeg = a*diagSqDiff + u(k)*u(k) + u(k+1)*u(k+1);
        const Real c = u(k)*u(k)*diagSqDiff;

        // Compute tau = sigmaEst^2 - d(k)^2
        Real tau = SolveQuadratic( a, bNeg, c ); 
 
        sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue( tau, d(origin) );

        // TODO(poulson): Document this LAPACK heuristic
        const Real eps = limits::Epsilon<Real>();
        const Real sqrtEps = Sqrt(eps);
        if( (d(k) <= sqrtEps*d(k+1)) &&
            (Abs(u(k)) <= sqrtEps) &&
            (d(k) > Real(0)) )
        {
            geometricAverage = true;
            sigmaRelEst = Min( 10*d(k), sigmaRelUpperBound );
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
        originOnLeft = false;
        origin = k+1;

        // Again follow LAPACK's {s,d}lasd4's [CITATION] lead in carefully
        // computing centerRoot - d(k+1).
        sigmaRelLowerBound = -diagSqDiffHalf / (d(origin) + centerRoot);
        sigmaRelUpperBound = 0; 
        const Real bNeg = -a*diagSqDiff + u(k)*u(k) + u(k+1)*u(k+1);
        const Real c = -u(k)*u(k)*diagSqDiff;

        const Real tau = SolveQuadratic( a, bNeg, c );

        sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue( tau, d(origin) );
    }

    sigmaEst = sigmaRelEst + d(origin);
    for( Int j=0; j<n; ++j ) 
    {
        dPlusShift(j) = (d(j) + d(origin)) + sigmaRelEst;
        dMinusShift(j) = (d(j) - d(origin)) - sigmaRelEst;
    }
}

// For the case where k = n-1. 
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularInitialGuessLast
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
  const Real& diagSqDiff,
  Matrix<Real>& dPlusShift,
  Matrix<Real>& dMinusShift,
  Real& sigmaEst,
  Real& sigmaRelEst )
{
    DEBUG_CSE
    const Int n = d.Height();
    const Real zero(0), one(1);

    // Since the largest eigenvalue must live within
    //
    //   d(n-1)^2 + (0, rho ||u||_2^2) = d(n-1)^2 + (0, rho),
    //
    // we use the center as our initial guess. Therefore, put
    // 
    //   tau = (rho ||u||_2^2) / 2 = rho / 2 = sigmaEst^2 - d(n-1)^2.
    //
    Real tau = rho / 2;

    sigmaRelEst = RelativeEigenvalueToRelativeSingularValue( tau, d(n-1) );

    sigmaEst = sigmaRelEst + d(n-1);
    for( Int j=0; j<n; ++j ) 
    {
        dPlusShift(j) = (d(j) + d(n-1)) + sigmaRelEst;
        dMinusShift(j) = (d(j) - d(n-1)) - sigmaRelEst;
    }

    const Int origin = n-2;

    // Form psi_{m-1}
    Real psiMinus = zero;
    for( Int j=0; j<n-2; ++j )
        psiMinus += u(j)*u(j) / (dPlusShift(j)*dMinusShift(j));

    const Real phi = u(n-1)*u(n-1) / (dPlusShift(n-1)*dMinusShift(n-1));

    const Real a = one/rho + psiMinus;
    const Real secularMinus = a + phi;
    const Real secular = secularMinus +
      u(origin)*u(origin) / (dPlusShift(origin)*dMinusShift(origin));

    if( secular <= zero )
    {
        const Real singValUpperBound = Sqrt(d(n-1)*d(n-1) + rho);
        // TODO(poulson): Document this branch from LAPACK's {s,d}lasd4
        // [CITATION]
        const Real temp =
          u(n-2)*u(n-2) /
           ((d(n-2) + singValUpperBound)*
            (d(n-1) - d(n-2) + rho/(d(n-1) + singValUpperBound))) +
          u(n-1)*u(n-1) / rho;

        if( a <= temp )
        {
            tau = rho;
        }
        else
        {
            const Real bNeg = -a*diagSqDiff + u(n-2)*u(n-2) + u(n-1)*u(n-1);
            const Real c = -u(n-1)*u(n-1)*diagSqDiff;
            tau = SolveQuadratic( a, bNeg, c );
            sigmaRelEst =
              RelativeEigenvalueToRelativeSingularValue( tau, d(n-1) );
        }
    }
    else
    {
        const Real bNeg = -a*diagSqDiff + u(n-2)*u(n-2) + u(n-1)*u(n-1);
        const Real c = -u(n-1)*u(n-1)*diagSqDiff;
        tau = SolveQuadratic( a, bNeg, c );
        sigmaRelEst = RelativeEigenvalueToRelativeSingularValue( tau, d(n-1) );
    }

    sigmaEst = sigmaRelEst + d(n-1);
    for( Int j=0; j<n; ++j ) 
    {
        dPlusShift(j) = (d(j) + d(n-1)) + sigmaRelEst;
        dMinusShift(j) = (d(j) - d(n-1)) - sigmaRelEst;
    }
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void EvaluateSecular
( const Real& rho,
  const Matrix<Real>& u,
  Int origin,
  const Matrix<Real>& dPlusShift,
  const Matrix<Real>& dMinusShift,
  Real& secular,
  Real& secularDeriv,
  Real& secularMinus,
  Real& secularMinusDeriv,
  Real& psiMinus,
  Real& psiMinusDeriv,
  Real& phi,
  Real& phiDeriv,
  Real& relErrorBound )
{
    DEBUG_CSE
    const Real one(1);
    const Int n = u.Height();
    const Real rhoInv = one / rho;
    Real temp;

    // Compute psi_m and psi'_m, where m is the origin index, and an
    // approximation of the error
    // (see Ren-Cang Li, "Solving Secular Equations Stably and Efficiently",
    // LAPACK Working Note 89, 1993 [CITATION], as well as LAPACK's 
    // {s,d}lasd4 [CITATION]). The loop direction is chosen to heuristically
    // sum from small to large components.
    psiMinus = 0;
    psiMinusDeriv = 0;
    relErrorBound = 0;
    for( Int j=0; j<origin; ++j )
    {
        // Compute the j'th term divided by u(j)
        temp = u(j) / (dPlusShift(j)*dMinusShift(j));
        psiMinus += u(j)*temp;
        psiMinusDeriv += temp*temp;
        relErrorBound += psiMinus;
    }
    relErrorBound = Abs(relErrorBound);

    // Compute phi, its derivative, and accumulate an approximation of the
    // error in both psi_{m-1} and phi_m, where m is the origin index. The loop
    // direction is chosen to heuristically sum from small to large components.
    phi = 0;
    phiDeriv = 0;
    for( Int j=n-1; j>origin; --j )
    {
        // Compute the j'th term divided by u(j)
        temp = u(j) / (dPlusShift(j)*dMinusShift(j));
        phi += u(j)*temp;
        phiDeriv += temp*temp;
        relErrorBound += phi;
    }

    // Compute the secular function with the origin term removed
    // (and its derivative)
    secularMinus = rhoInv + psiMinus + phi;
    secularMinusDeriv = psiMinusDeriv + phiDeriv; 

    // Cf. LAPACK's {s,d}lasd4 [CITATION] for this computational strategy
    temp = u(origin) / (dPlusShift(origin)*dMinusShift(origin));
    secularDeriv = (psiMinusDeriv + temp*temp) + phiDeriv;
    temp *= u(origin);
    secular = secularMinus + temp;
    relErrorBound += 8*(phi-psiMinus) + 2*rhoInv + 3*Abs(temp);
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void EvaluateSecularLast
( const Real& rho,
  const Matrix<Real>& u,
  const Matrix<Real>& dPlusShift,
  const Matrix<Real>& dMinusShift,
  Real& secular,
  Real& psiDeriv,
  Real& phiDeriv,
  Real& relErrorBound )
{
    DEBUG_CSE
    const Real one(1);
    const Int n = u.Height();
    const Real rhoInv = one / rho;
    Real temp;

    Real psi = 0;
    psiDeriv = 0;
    relErrorBound = 0;
    for( Int j=0; j<n-1; ++j )
    {
        // Compute the j'th term divided by u(j)
        temp = u(j) / (dPlusShift(j)*dMinusShift(j));
        psi += u(j)*temp;
        psiDeriv += temp*temp;
        relErrorBound += psi;
    }
    relErrorBound = Abs(relErrorBound);

    Real phi = 0;
    phiDeriv = 0;
    // Compute the (n-1)'th term divided by u(n-1)
    temp = u(n-1) / (dPlusShift(n-1)*dMinusShift(n-1));
    phi += u(n-1)*temp;
    phiDeriv += temp*temp;
    relErrorBound -= phi;

    // Compute the secular function with the origin term removed
    // (and its derivative)
    secular = rhoInv + psi + phi;

    relErrorBound += 8*(-phi-psi) + rhoInv;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
struct SecularInfo
{
    Real singularValue;
    Int numIterations = 0;
    Int numAlternations = 0;
    Int numCubicFailures = 0;
};

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SecularUpdate
( Int origin,
  bool originOnLeft,
  const Matrix<Real>& d,
  const Matrix<Real>& u,
  const Real& diagSqDiff,
  const Matrix<Real>& dPlusShift,
  const Matrix<Real>& dMinusShift,
  const Real& secular,
  const Real& secularDeriv,
  const Real& secularMinus,
  const Real& secularMinusDeriv,
  const Real& psiMinus,
  const Real& psiMinusDeriv,
  const Real& phi,
  const Real& phiDeriv,
  const Real& sigmaEst,
  const Real& sigmaRelEst,
  const Real& sigmaRelLowerBound,
  const Real& sigmaRelUpperBound,
  bool initialize,
  bool geometricAverage,
  bool alternateStrategy,
  bool& useThreePoles,
  SecularInfo<Real>& info )
{
    DEBUG_CSE
    const Real zero(0);
    const Int k = origin;
    Real eta;

    if( useThreePoles )
    {
        // Use the "Hybrid Scheme" described in subsection 3.4 of LAWN 89 and
        // implemented within {s,d}lasd4.
 
        // Carefully compute
        //
        //   leftGap = d(origin-1)^2 - sigmaEst^2, and
        //   rightGap = d(origin+1)^2 - sigmaEst^2.
        //
        const Real leftGap = dPlusShift(origin-1)*dMinusShift(origin-1);
        const Real rightGap = dPlusShift(origin+1)*dMinusShift(origin+1);

        // Carefully compute doubleGap = d(origin+1)^2 - d(origin-1)^2
        const Real doubleGap =
          (d(origin+1)-d(origin-1))*(d(origin+1)+d(origin-1));

        Real a;
        Matrix<Real> zCubic(3,1);
        if( alternateStrategy )
        {
            a = secularMinus - leftGap*psiMinusDeriv - rightGap*phiDeriv;
            zCubic(0) = leftGap*leftGap*psiMinusDeriv;
            zCubic(1) = u(origin)*u(origin);
            zCubic(2) = rightGap*rightGap*phiDeriv;
        }
        else
        {
            if( originOnLeft )
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
                const Real a = (secularMinus-rightGap*secularMinusDeriv) +
                  leftDerivTerm*doubleGap;

                const Real psiMinusTwoDeriv =
                  Max( psiMinusDeriv-leftDerivTerm, zero );
    
                zCubic(0) = u(origin-1)*u(origin-1);
                zCubic(1) = u(origin)*u(origin);
                zCubic(2) = rightGap*rightGap*psiMinusTwoDeriv;
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
                const Real a = (secularMinus-leftGap*secularMinusDeriv) -
                  mPlusDerivTerm*doubleGap;

                const Real phiMinusDeriv = Max( phiDeriv-mPlusDerivTerm, zero );
            
                zCubic(0) = leftGap*leftGap*phiMinusDeriv;
                zCubic(1) = u(origin)*u(origin);
                zCubic(2) = u(origin+1)*u(origin+1);
            }
        }
        Matrix<Real> dCubic(3,1);
        dCubic(0) = leftGap;
        dCubic(1) = dPlusShift(origin)*dMinusShift(origin);
        dCubic(2) = rightGap;

        auto cubicInfo =
          CubicSecular( originOnLeft, a, zCubic, dCubic, secular, initialize );
        if( cubicInfo.converged )
        {
            eta = cubicInfo.root;
        }
        else
        {
            ++info.numCubicFailures;
            useThreePoles = false;
            const Real kGap = dPlusShift(k)*dMinusShift(k);
            const Real kp1Gap = dPlusShift(k+1)*dMinusShift(k+1);
            Real a;
            if( originOnLeft )
            {
                const Real temp = u(k) / kGap;
                a = secular - kp1Gap*secularDeriv + diagSqDiff*(temp*temp);
            }
            else
            {
                const Real temp = u(k+1) / kp1Gap;
                a = secular - kGap*secularDeriv + diagSqDiff*(temp*temp);
            }
            Real bNeg = (kGap+kp1Gap)*secular - kGap*kp1Gap*secularDeriv;
            const Real c = kGap*kp1Gap*secular;
            if( a == zero && bNeg == zero )
            {
                // a x^2 + b x + c = 0 has collapsed to c = 0, which is
                // nonsense. Follow LAPACK's {s,d}lasd4 [CITATION] in handling
                // this breakdown.
                if( originOnLeft )
                {
                    bNeg = u(k)*u(k) + kp1Gap*kp1Gap*secularMinusDeriv;
                }
                else
                {
                    bNeg = u(k+1)*u(k+1) + kGap*kGap*secularMinusDeriv; 
                }
            }
            eta = SolveQuadratic( a, bNeg, c );
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
        const Real kGap = dPlusShift(k)*dMinusShift(k);
        const Real kp1Gap = dPlusShift(k+1)*dMinusShift(k+1);

        Real a;
        if( alternateStrategy )
        {
            // TODO(poulson): Document this strategy
            if( originOnLeft )
            {
                const Real temp = u(k) / kGap;
                a = secular - kGap*(psiMinusDeriv+temp*temp) - kp1Gap*phiDeriv;
            }
            else
            {
                const Real temp = u(k+1) / kp1Gap;
                a = secular - kGap*psiMinusDeriv - kp1Gap*(phiDeriv+temp*temp);
            }
        }
        else
        {
            if( originOnLeft )
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
                a = secular - kp1Gap*secularDeriv + (temp*temp)*diagSqDiff;
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
                a = secular - kGap*secularDeriv - (temp*temp)*diagSqDiff;
            }
        }

        // See Proposition 3 from LAWN 89 [CITATION] for these formulae for the
        // coefficients of the quadratic equation a x^2 - bNeg x + c = 0.
        Real bNeg = (kGap+kp1Gap)*secular - kGap*kp1Gap*secularDeriv;
        Real c = kGap*kp1Gap*secular;

        if( a == zero && bNeg == zero )
        {
            // a eta^2 + b eta + c = 0 has collapsed to c = 0, which is
            // nonsensical to solve for eta. We therefore follow the lead of
            // LAPACK's {s,d}lasd4 [CITATION] and use a mysterious patch-up.
            //
            // TODO(poulson): Provide motivation for these formulae.
            if( alternateStrategy )
            {
                if( originOnLeft )
                {
                    const Real temp = u(k) / kGap;
                    bNeg = kGap*kGap*(psiMinusDeriv+temp*temp) +
                           kp1Gap*kp1Gap*phiDeriv;
                }
                else
                {
                    const Real temp = u(k+1) / kp1Gap;
                    bNeg = kGap*kGap*psiMinusDeriv +
                           kp1Gap*kp1Gap*(phiDeriv+temp*temp);
                }
            }
            else
            {
                if( originOnLeft )
                    bNeg = u(k)*u(k) + kp1Gap*kp1Gap*secularMinusDeriv;
                else
                    bNeg = u(k+1)*u(k+1) + kGap*kGap*secularMinusDeriv;
            }
        }
        eta = SolveQuadratic( a, bNeg, c );
    }

    if( secular*eta >= zero )
    {
        // Use a Newton step instead
        eta = -secular / secularDeriv;
    }

    // Convert eta from an eigenvalue update to singular value update
    eta = RelativeEigenvalueToRelativeSingularValue( eta, sigmaEst );

    // Since sigmaRelEst = sigmaEst - d(origin), the following sets
    // sigmaRelEstNew = sigmaEstNew - d(origin).
    const Real sigmaRelEstNew = sigmaRelEst + eta;
    if( sigmaRelEstNew > sigmaRelUpperBound ||
        sigmaRelEstNew < sigmaRelLowerBound )
    {
        // We again follow LAPACK's {s,d}lasd4 for how to handle this breakdown
        if( secular < zero )
            eta = (sigmaRelUpperBound - sigmaRelEst) / 2;
        else
            eta = (sigmaRelLowerBound - sigmaRelEst) / 2;
        if( geometricAverage )
        {
            if( secular < zero )
            {
                if( sigmaRelEst > zero )
                    eta = Sqrt(sigmaRelUpperBound*sigmaRelEst) - sigmaRelEst;
            }
            else
            {
                if( sigmaRelLowerBound > zero )
                    eta = Sqrt(sigmaRelLowerBound*sigmaRelEst) - sigmaRelEst;
            }
        }
    }

    return eta;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real SecularUpdateLast
( const Real& rho,
  const Matrix<Real>& dPlusShift,
  const Matrix<Real>& dMinusShift,
  const Real& secular,
  const Real& psiDeriv,
  const Real& phiDeriv,
  const Real& sigmaEst,
  bool initialize )
{
    DEBUG_CSE
    const Real zero(0);
    const Int n = dPlusShift.Height();

    const Real kGap = dPlusShift(n-1)*dMinusShift(n-1);
    const Real km1Gap = dPlusShift(n-2)*dMinusShift(n-2);

    Real a = secular - km1Gap*psiDeriv - kGap*phiDeriv;

    Real eta;
    if( initialize )
    {
        // TODO(poulson): Explain this special case

        if( a < zero )
            a = -a;

        if( a == zero )
        {
            eta = rho - sigmaEst*sigmaEst;
        }
        else
        {
            const Real bNeg =
              (kGap+km1Gap)*secular - kGap*km1Gap*(psiDeriv+phiDeriv);
            const Real c = kGap*km1Gap*secular;
            eta = SolveQuadratic( a, bNeg, c );
        }
    }
    else
    {
        const Real bNeg =
          (kGap+km1Gap)*secular - kGap*km1Gap*(psiDeriv+phiDeriv);
        const Real c = kGap*km1Gap*secular;
        eta = SolveQuadratic( a, bNeg, c );
    }

    if( secular*eta >= zero )
    {
        // Use a Newton step instead
        eta = -secular / (psiDeriv+phiDeriv);
    }

    if( initialize )
    {
        // eta = Min( eta, rho+kGap )
        if( eta-kGap > rho )
            eta = rho + kGap;
    }
    else
    {
        if( eta-kGap < zero )
            eta /= 2;
    }

    // Convert eta from an eigenvalue update to singular value update
    eta = RelativeEigenvalueToRelativeSingularValue( eta, sigmaEst );

    return eta;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularInfo<Real>
SecularInner
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        Matrix<Real>& dMinusShift,
        Matrix<Real>& dPlusShift )
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

    const Int maxIterations = 400; // This is LAPACK's limit

    SecularInfo<Real> info;
    dMinusShift.Resize(n,1);
    dPlusShift.Resize(n,1);

    const Real diagSum = d(k+1) + d(k);
    const Real diagDiff = d(k+1) - d(k);
    const Real diagSqDiff = diagSum*diagDiff; // d(k+1)^2 - d(k)^2

    Real sigmaEst, sigmaRelEst, sigmaRelLowerBound, sigmaRelUpperBound;
    Int origin;
    bool originOnLeft;
    bool geometricAverage;
    SecularInitialGuess
    ( k, d, rho, u, diagSqDiff,
      dPlusShift, dMinusShift,
      sigmaEst, sigmaRelEst, sigmaRelLowerBound, sigmaRelUpperBound,
      origin, originOnLeft, geometricAverage );

    Real secular, secularDeriv,
         secularMinus, secularMinusDeriv,
         psiMinus, psiMinusDeriv, phi, phiDeriv, relErrorBound;
    EvaluateSecular
    ( rho, u, origin, dPlusShift, dMinusShift,
      secular, secularDeriv,
      secularMinus, secularMinusDeriv,
      psiMinus, psiMinusDeriv,
      phi, phiDeriv,
      relErrorBound );
    ++info.numIterations;

    // See the discussion in LAWN 89 [CITATION] on the usage of three poles
    bool useThreePoles =
      (originOnLeft ? secularMinus<zero : secularMinus>zero);
    if( origin == 0 || origin == n-1 )
        useThreePoles = false;

    // Check if we have already converged
    if( Abs(secular) <= eps*relErrorBound )
    {
        info.singularValue = sigmaEst;
        return info;
    }

    if( secular <= zero )
    {
        // Mark the sought eigenvalue as right of our current estimate
        sigmaRelLowerBound = Max( sigmaRelLowerBound, sigmaRelEst );
    }
    else
    {
        // Mark the sought eigenvalue as left of our current estimate
        sigmaRelUpperBound = Min( sigmaRelUpperBound, sigmaRelEst );
    }

    // Compute the first update to our estimate
    bool initialize = true;
    bool alternateStrategy = false;
    Real eta = SecularUpdate
      ( origin, originOnLeft,
        d, u, diagSqDiff, dPlusShift, dMinusShift,
        secular, secularDeriv,
        secularMinus, secularMinusDeriv,
        psiMinus, psiMinusDeriv, phi, phiDeriv,
        sigmaEst, sigmaRelEst, sigmaRelLowerBound, sigmaRelUpperBound,
        initialize, geometricAverage, alternateStrategy, useThreePoles, info );
    ++info.numIterations;

    sigmaEst += eta;
    sigmaRelEst += eta;
    Real secularOld = secular;
    for( Int j=0; j<n; ++j )
    {
        dPlusShift(j) += eta;
        dMinusShift(j) -= eta;
    }

    EvaluateSecular
    ( rho, u, origin, dPlusShift, dMinusShift,
      secular, secularDeriv,
      secularMinus, secularMinusDeriv,
      psiMinus, psiMinusDeriv,
      phi, phiDeriv,
      relErrorBound );
    
    // This strategy was described in Ren-Cang Li's LAWN 89
    if( originOnLeft )
    {
        if( -secular > Abs(secularOld) / 10 )
        {
            alternateStrategy = true;        
            ++info.numAlternations;
        }
    }
    else
    {
        if( secular > Abs(secularOld) / 10 )
        {
            alternateStrategy = true;
            ++info.numAlternations;
        }
    }

    initialize = false;
    while( true )
    {
        if( Abs(secular) <= eps*relErrorBound ) 
        {
            break;
        }
        if( info.numIterations >= maxIterations )
        {
            RuntimeError
            ("Secular solver did not converge in ",maxIterations," iterations");
        }

        if( secular <= zero )
            sigmaRelLowerBound = Max( sigmaRelLowerBound, sigmaRelEst );
        else
            sigmaRelUpperBound = Min( sigmaRelUpperBound, sigmaRelEst );

        // Decide the next step
        eta = SecularUpdate
          ( origin, originOnLeft,
            d, u, diagSqDiff, dPlusShift, dMinusShift,
            secular, secularDeriv,
            secularMinus, secularMinusDeriv,
            psiMinus, psiMinusDeriv, phi, phiDeriv,
            sigmaEst, sigmaRelEst, sigmaRelLowerBound, sigmaRelUpperBound,
            initialize,
            geometricAverage, alternateStrategy, useThreePoles, info );
        ++info.numIterations;

        sigmaEst += eta;
        sigmaRelEst += eta;
        secularOld = secular;
        for( Int j=0; j<n; ++j )
        {
            dPlusShift(j) += eta;
            dMinusShift(j) -= eta;
        }

        EvaluateSecular
        ( rho, u, origin, dPlusShift, dMinusShift,
          secular, secularDeriv,
          secularMinus, secularMinusDeriv,
          psiMinus, psiMinusDeriv,
          phi, phiDeriv,
          relErrorBound );
    
        // This strategy was described in Ren-Cang Li's LAWN 89
        if( secular*secularOld > zero && Abs(secular) > Abs(secularOld)/10 )
        {
            alternateStrategy = !alternateStrategy;
            ++info.numAlternations;
        }
    }
    info.singularValue = sigmaEst;
    return info;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularInfo<Real>
SecularLast
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        Matrix<Real>& dMinusShift,
        Matrix<Real>& dPlusShift )
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

    // Cf. LAPACK's {s,d}lasd4 for this choice
    const Int maxIterations = 400;

    SecularInfo<Real> info;
    dMinusShift.Resize(n,1);
    dPlusShift.Resize(n,1);

    const Real diagSum = d(k+1) + d(k);
    const Real diagDiff = d(k+1) - d(k);
    const Real diagSqDiff = diagSum*diagDiff; // d(k+1)^2 - d(k)^2

    Real sigmaEst, sigmaRelEst;
    SecularInitialGuessLast
    ( d, rho, u, diagSqDiff, dPlusShift, dMinusShift, sigmaEst, sigmaRelEst );

    Real secular, psiDeriv, phiDeriv, relErrorBound;
    EvaluateSecularLast
    ( rho, u, dPlusShift, dMinusShift,
      secular, psiDeriv, phiDeriv, relErrorBound );
    ++info.numIterations;

    if( Abs(secular) <= eps*relErrorBound )
    {
        info.singularValue = sigmaEst;
        return info;
    }

    // Calculate the first update
    bool initialize = true;
    Real eta =
      SecularUpdateLast
      ( rho, dPlusShift, dMinusShift,
        secular, psiDeriv, phiDeriv, sigmaEst, initialize );
    ++info.numIterations;

    sigmaEst += eta;
    sigmaRelEst += eta;
    Real secularOld = secular;
    for( Int j=0; j<n; ++j )
    {
        dPlusShift(j) += eta;
        dMinusShift(j) -= eta;
    }

    EvaluateSecularLast
    ( rho, u, dPlusShift, dMinusShift,
      secular, psiDeriv, phiDeriv, relErrorBound );

    initialize = false;
    while( true )
    {
        if( Abs(secular) <= eps*relErrorBound ) 
        {
            break;
        }
        if( info.numIterations >= maxIterations )
        {
            RuntimeError
            ("Secular solver did not converge in ",maxIterations," iterations");
        }

        // Decide the next step
        eta = SecularUpdateLast
          ( rho, dPlusShift, dMinusShift,
            secular, psiDeriv, phiDeriv, sigmaEst, initialize );
        ++info.numIterations;

        sigmaEst += eta;
        sigmaRelEst += eta;
        secularOld = secular;
        for( Int j=0; j<n; ++j )
        {
            dPlusShift(j) += eta;
            dMinusShift(j) -= eta;
        }

        EvaluateSecularLast
        ( rho, u, dPlusShift, dMinusShift,
          secular, psiDeriv, phiDeriv, relErrorBound );
    }
    info.singularValue = sigmaEst;
    return info;
}

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
template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularInfo<Real>
Secular
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        Matrix<Real>& dMinusShift,
        Matrix<Real>& dPlusShift )
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

    SecularInfo<Real> info;
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
          TwoByTwoSecular
          ( k, d(0), d(1), rho, u(0), u(1),
            dMinusShift(0), dMinusShift(1),
            dPlusShift(0), dPlusShift(1) );
        return info;
    }

    if( k < n-1 )
    {
        return SecularInner( k, d, rho, u, dMinusShift, dPlusShift );
    }
    else
    {
        return SecularLast( k, d, rho, u, dMinusShift, dPlusShift );
    }
}

} // namespace dc
} // namespace bidiag_svd
} // namespace El

#endif // ifndef EL_BIDIAG_SVD_DC_SECULAR_HPP
