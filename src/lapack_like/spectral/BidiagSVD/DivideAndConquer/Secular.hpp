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
struct SecularInfo
{
    Real singularValue;
    Int numIterations = 0;
};

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
    const Real zero(0), one(1), two(2), three(3), four(4), eight(8), ten(10);
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichSingularValue;
    const Int n = d.Height();
    DEBUG_ONLY(
      if( k == n-1 )
          LogicError("SecularInner meant for inner singular values");
      if( n <= 2 )
          LogicError("SecularInner meant for n > 2");
    )

    SecularInfo<Real> info;
    dMinusShift.Resize(n,1);
    dPlusShift.Resize(n,1);

    const Real rhoInv = one / rho;
    const Real diagSum = d(k+1) + d(k);
    const Real diagDiff = d(k+1) - d(k);
    const Real diagSqDiff = diagSum*diagDiff; // d(k+1)^2 - d(k)^2
    const Real diagSqDiffHalf = diagSqDiff / two;

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
    const Real center = (d(k)*d(k) + d(k+1)*d(k+1)) / two; 
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
    const Real a = rhoInv + psiMinus + phiMinus;
    const Real secularCenter = a +
      u(k)*u(k) / (dPlusShift(k)*dMinusShift(k)) +
      u(k+1)*u(k+1) / (dPlusShift(k+1)*dMinusShift(k+1));

    bool originOnLeft;
    Int origin; 
    Real sigmaRelLowerBound, sigmaRelUpperBound;
    Real sigmaRelEst;
    bool geometricAverage = false; // TODO(poulson): Documentation?
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
        //   b = -a Delta - (u(k)^2 + u(k+1)^2),
        //   c = u(k)^2 Delta,
        //
        // with Delta = d(k+1)^2 - d(k)^2.
        originOnLeft = true;
        origin = k;
        sigmaRelLowerBound = zero; 
        sigmaRelUpperBound = shiftedCenterRoot;
        const Real bNeg = a*diagSqDiff + u(k)*u(k) + u(k+1)*u(k+1);
        const Real c = u(k)*u(k)*diagSqDiff;

        // Just as in the case of TwoByTwoSecular, we differ from LAPACK
        // (in this case, {s,d}lasd4 [CITATION]) by clipping the
        // discriminant at zero rather than taking its absolute value.
        const Real discrim = Max( bNeg*bNeg-four*a*c, zero );

        // Compute tau = sigmaEst^2 - d(k)^2
        Real tau;
        if( bNeg <= zero )
        {
            // Use the standard quadratic equation formula
            tau = (bNeg - Sqrt(discrim)) / (two*a); 
        } 
        else
        {
            // Use the inverted quadratic equation formula
            tau = 2*c / (bNeg + Sqrt(discrim));
        }

        // See the notes in TwoByTwoSecular for a step-by-step explanation
        // of how this converts tau = sigmaEst^2 - d(k)^2 into
        // sigmaEst - d(k) (though it is fairly trivial).
        sigmaRelEst = tau / (d(origin) + Sqrt(d(origin)*d(origin) + tau));

        // TODO(poulson): Document this LAPACK heuristic
        const Real sqrtEps = Sqrt(eps);
        if( (d(k) <= sqrtEps*d(k+1)) &&
            (Abs(u(k)) <= sqrtEps) &&
            (d(k) > zero) )
        {
            geometricAverage = true;
            sigmaRelEst = Min( ten*d(k), sigmaRelUpperBound );
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
        //   b = a Delta - (u(k)^2 + u(k+1)^2),
        //   c = -u(k)^2 Delta,
        //
        // with Delta = d(k+1)^2 - d(k)^2.
        originOnLeft = false;
        origin = k+1;
        // Again follow LAPACK's {s,d}lasd4's [CITATION] lead in carefully
        // computing centerRoot - d(k+1).
        sigmaRelLowerBound = -diagSqDiffHalf / (d(origin) + centerRoot);
        sigmaRelUpperBound = zero; 
        const Real b = a*diagSqDiff - u(k)*u(k) - u(k+1)*u(k+1);
        const Real cNeg = u(k)*u(k)*diagSqDiff;
        // Just as in the case of TwoByTwoSecular, we differ from LAPACK
        // (in this case, {s,d}lasd4 [CITATION]) by clipping the
        // discriminant at zero rather than taking its absolute value.
        const Real discrim = Max( b*b+four*a*cNeg, zero );

        Real tau;
        if( b >= zero )
        {
            // Use the standard quadratic formula
            tau = (-b - Sqrt(discrim)) / (two*a);
        }
        else
        {
            // Use the inverted quadratic formula
            tau = two*cNeg / (b - Sqrt(discrim));
        }

        // See the notes in TwoByTwoSecular for a step-by-step explanation
        // of how this converts tau = sigmaEst^2 - d(k+1)^2 into
        // sigmaEst - d(k+1) (though it is fairly trivial).
        sigmaRelEst = tau / (d(origin) + Sqrt(d(origin)*d(origin) + tau));
    }

    Real sigmaEst = sigmaRelEst + d(origin);
    for( Int j=0; j<n; ++j ) 
    {
        dPlusShift(j) = (d(j) + d(origin)) + sigmaRelEst;
        dMinusShift(j) = (d(j) - d(origin)) - sigmaRelEst;
    }

    // Compute psiMinus, its derivative, and an approximation of the error
    // (see Ren-Cang Li, "Solving Secular Equations Stably and Efficiently",
    // LAPACK Working Note 89, 1993 [CITATION], as well as LAPACK's 
    // {s,d}lasd4 [CITATION]). The loop direction is chosen to heuristically
    // sum from small to large components.
    psiMinus = zero;
    Real psiMinusDeriv = zero;
    Real relErrorBound = zero;
    for( Int j=0; j<origin; ++j )
    {
        // Compute the j'th term divided by u(j)
        const Real temp = u(j) / (dPlusShift(j)*dMinusShift(j));
        psiMinus += u(j)*temp;
        psiMinusDeriv += temp*temp;
        relErrorBound += psiMinus;
    }
    relErrorBound = Abs(relErrorBound);

    // Compute phi, its derivative, and accumulate an approximation of the
    // error in both psiMinus and phi. The loop direction is chosen to 
    // heuristically sum from small to large components.
    Real phi = zero;
    Real phiDeriv = zero;
    for( Int j=n-1; j>origin; --j )
    {
        // Compute the j'th term divided by u(j)
        const Real temp = u(j) / (dPlusShift(j)*dMinusShift(j));
        phi += u(j)*temp;
        phiDeriv += temp*temp;
        relErrorBound += phi;
    }

    // Compute the secular function with the origin term removed
    const Real secularMinus = rhoInv + psiMinus + phi;

    // See the discussion in LAWN 89 [CITATION] on the usage of three poles
    bool useThreePoles =
      (originOnLeft ? secularMinus<zero : secularMinus>zero);
    if( origin == 0 || origin == n-1 )
        useThreePoles = false;

    // Cf. LAPACK's {s,d}lasd4 [CITATION] for this computational strategy
    Real temp = u(origin) / (dPlusShift(origin)*dMinusShift(origin));
    const Real secularDeriv = (psiMinusDeriv + temp*temp) + phiDeriv;
    temp *= u(origin);
    const Real secular = secularMinus + temp;
    relErrorBound += eight*(phi-psiMinus) + two*rhoInv + three*Abs(temp);

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
    ++info.numIterations;
    Real eta;
    if( useThreePoles )
    {
        // Use the "Hybrid Scheme" described in subsection 3.4 of LAWN 89 and
        // implemented within {s,d}lasd4.
 
        // Carefully compute
        //
        //   originMinusGap = d(origin-1)^2 - sigmaEst^2, and
        //   originPlusGap = d(origin+1)^2 - sigmaEst^2.
        //
        const Real originMinusGap = dPlusShift(origin-1)*dMinusShift(origin-1);
        const Real originPlusGap = dPlusShift(origin+1)*dMinusShift(origin+1);

        // Carefully compute doubleGap = d(origin+1)^2 - d(origin-1)^2
        const Real doubleGap =
          (d(origin+1)-d(origin-1))*(d(origin+1)+d(origin-1));

        // Quickly compute f_m'(y) = f(y) - u(m)^2 / (d(m)^2 - y), where m is
        // k if 'originOnLeft' and k+1 otherwise.
        const Real secularMinusDeriv = psiMinusDeriv + phiDeriv;

        Real c, numerator0, numerator1, numerator2;
        if( originOnLeft )
        {
            // Since the shift origin, m, is k, We will interpolate the secular
            // equation as
            //
            //  Q(x; a, s, S) = a + u(m-1)^2 / (d(m-1)^2 - x) + 
            //                      u(m  )^2 / (d(m  )^2 - x) + 
            //                      S        / (d(m+1)^2 - x),
            //
            // which can be effected by applying the Middle Way, the Fixed
            // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme to
            //
            //   f_m(y) = f(y) - u(m)^2 / (d(m)^2 - x).
            //
            // We balance between the Fixed Weight Method and the Middle Way;
            // their formulae for computing a are identical. In both cases,
            // 
            //  a = f_m-rightGap*f_m'+(u(m-1)/leftGap)^2*(d(m+1)^2-d(m-1)^2).
            //
            // The computation of S primarily follows the Fixed Weight Method
            // but enforces the non-negativity of psi_{k-2}' when computing
            //
            //   f'_{m-1,m}(y) = psi'_{m-2}(y) + phi'_m(y).
            //
            // The two subscripts on f denote the m-1 and m terms of the secular
            // equation being left out.
            //
            const Real leftRatio = u(m-1) / originMinusGap;
            const Real mMinusDerivTerm = leftRatio*leftRatio;
            const Real c = (secularMinus-rightGap*secularMinusDeriv) +
              mMinusDerivTerm*doubleGap;

            const Real psiMinusTwoDeriv =
              Max( psiMinusDeriv-mMinusDerivTerm, zero );

            numerator0 = u(origin-1)*u(origin-1);
            numerator1 = u(origin)*u(origin);
            numerator2 = rightGap*rightGap*psiMinusTwoDeriv;
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
            // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme to
            //
            //   f_m(y) = f(y) - u(m)^2 / (d(m)^2 - x).
            //
            // We balance between the Fixed Weight Method and the Middle Way;
            // their formulae for computing a are identical. In both cases,
            // 
            //  a = f_m-leftGap*f_m'-(u(m+1)/rightGap)^2*(d(m+1)^2-d(m-1)^2),
            //
            // The computation of s primarily follows the Fixed Weight Method
            // but enforces the non-negativity of phi_{m+1}' when computing
            //
            //   f'_{m,m+1}(y) = psi'_{m-1}(y) + phi'_{m+1}(y).
            //
            // The two subscripts on f denote the m and m+1 terms of the secular
            // equation being left out.
            //
            const Real rightRatio = u(origin+1) / originPlusGap;
            const Real mPlusDerivTerm = rightRatio*rightRatio;
            const Real c = (secularMinus-leftGap*secularMinusDeriv) -
              mPlusDerivTerm*doubleGap;

            const Real phiMinusDeriv = Max( phiDeriv-mPlusDerivTerm, zero );
            
            numerator0 = leftGap*leftGap*phiMinusDeriv;
            numerator1 = u(origin)*u(origin);
            numerator2 = u(origin+1)*u(origin+1);
        }
        const Real pole0 = leftGap;
        const Real pole1 = dPlusShift(origin)*dMinusShift(origin);
        const Real pole2 = rightGap;

        // TODO: Call Borges/Gragg/Thornton/Warner scheme here
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
        const Real leftGap = dPlusShift(k)*dMinusShift(k);
        const Real rightGap = dPlusShift(k+1)*dMinusShift(k+1);

        Real a;
        if( originOnLeft )
        {
            // Apply Eq'ns (32) through (34) from LAWN 89 to set
            //
            //  a = f - rightGap f' + (u(k)^2/leftGap^2) (d(k+1)^2 - d(k)^2),
            //
            // and implicitly set
            //
            //  s = u(k)^2,
            //  S = rightGap^2 (f' - u(k)^2/leftGap^2),
            //
            // where y is the current estimate of the k'th eigenvalue and 
            // f and f' are the secular equation and its derivative evaluated
            // at y.
            const Real leftRatio = u(k) / leftGap;
            a = secular - rightGap*secularDeriv +
              (leftRatio*leftRatio)*diagSqDiff;
        }
        else
        {
            // Apply Eq'ns (35) through (37) from LAWN 89 to set
            //
            //  a = f - leftGap f' - (u(k+1)^2/rightGap^2) (d(k+1)^2 - d(k)^2),
            //
            // and implicitly set
            //
            //  s = leftGap^2 (f' - u(k+1)^2/rightGap^2),
            //  S = u(k+1)^2.
            //
            const Real rightRatio = u(k+1) / rightGap;
            a = secular - leftGap*secularDeriv -
              (rightRatio*rightRatio)*diagSqDiff;
        }

        // See Proposition 3 from LAWN 89 [CITATION] for these formulae for the
        // coefficients of the quadratic equation a x^2 - bNeg x + c = 0.
        Real bNeg = (leftGap+rightGap)*secular - leftGap*rightGap*secularDeriv;
        Real c = leftGap*rightGap*secular;

        if( a == zero )
        {
            // a eta^2 + b eta + c = 0 has collapsed to b eta + c = 0, which is 
            // obviously uniquely solved by eta = -c/b.
            if( bNeg == zero )
            {
                // b eta + c = 0 has collapsed to c = 0, which is nonsensical to
                // solve for eta. We therefore follow the lead of LAPACK's
                // {s,d}lasd4 [CITATION] and use a mysterious patch-up.
                //
                // TODO(poulson): Provide motivation for these formulae.
                if( originOnLeft )
                    bNeg = u(k)*u(k) +
                      rightGap*rightGap*(psiMinusDeriv+phiDeriv);
                else
                    bNeg = u(k+1)*u(k+1) +
                      leftGap*leftGap*(psiMinusDeriv+phiDeriv);
            }
            // Set eta = -c/b
            eta = c / bNeg; 
        }
        else
        {
            // As noted in many places, LAPACK's {s,d}lasd4 [CITATION] prefers
            // to compute the absolute value of b^2 - 4 a c, but we prefer to
            // clip up to zero.
            const Real discrim = Max( bNeg*bNeg-four*a*c, zero );
            if( bNeg <= zero )
            {
                // Apply the standard quadratic formula
                eta = (bNeg - Sqrt(discrim)) / (two*a);
            }
            else
            {
                // Apply the inverse quadratic formula to avoid cancellation 
                eta = two*c / (bNeg + Sqrt(discrim));
            }
        }
    }

    // TODO(poulson)

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
    const Real zero(0), one(1), two(2), three(3), four(4), eight(8), ten(10);
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichSingularValue;
    const Int n = d.Height();
    DEBUG_ONLY(
      if( k != n-1 )
          LogicError("SecularLast meant for largest singular value");
      if( n <= 2 )
          LogicError("SecularLast meant for n > 2");
    )

    SecularInfo<Real> info;
    dMinusShift.Resize(n,1);
    dPlusShift.Resize(n,1);

    const Real rhoInv = one / rho;
    const Real diagSum = d(k+1) + d(k);
    const Real diagDiff = d(k+1) - d(k);
    const Real diagSqDiff = diagSum*diagDiff; // d(k+1)^2 - d(k)^2
    const Real diagSqDiffHalf = diagSqDiff / two;

    // TODO(poulson)

    return info;
}


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
