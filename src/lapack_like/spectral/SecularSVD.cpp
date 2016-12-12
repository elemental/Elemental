/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// As discussed in FP_Consistency_070816.pdf, which is linked from
//
// <https://software.intel.com/en-us/articles/consistency-of-floating-point-results-using-the-intel-compiler>,
//
// the Intel compilers do not preserve evaluation order by default and may issue
// an FMA that ignores parentheses. The following pragmas avoid this behavior.
//
#ifdef __INTEL_COMPILER
#pragma float_control (precise, on)
#pragma float_control (source, on)
#endif

#include "./SecularSVD/TwoByTwo.hpp"

namespace El {

namespace secular_svd {

template<typename Real>
struct State
{
    Int origin;
    bool originOnLeft;

    Matrix<Real> dPlusShift;
    Matrix<Real> dMinusShift;

    Real diagSqDiff;

    Real rootRelEst;

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

    Real rootRelEst;

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

    // The origin term of psi and its derivative
    Real psiOrigin;
    Real psiOriginDeriv;

    Real relErrorBound;
};

template<typename Real,typename=EnableIf<IsReal<Real>>>
void EvaluateSecular
( const Real& rho,
  const Matrix<Real>& z,
        State<Real>& state,
  bool penalizeDerivative=false )
{
    EL_DEBUG_CSE
    const Real one(1);
    const Int n = z.Height();
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
        // Compute the j'th term divided by z(j)
        temp = z(j) / (state.dPlusShift(j)*state.dMinusShift(j));
        state.psiMinus += z(j)*temp; // This should be a negative contribution
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
        // Compute the j'th term divided by z(j)
        temp = z(j) / (state.dPlusShift(j)*state.dMinusShift(j));
        state.phi += z(j)*temp;
        state.phiDeriv += temp*temp;
        state.relErrorBound += state.phi;
    }

    // Compute the secular function with the origin term removed
    // (and its derivative)
    state.secularMinus = rhoInv + state.psiMinus + state.phi;
    state.secularMinusDeriv = state.psiMinusDeriv + state.phiDeriv;

    // Cf. LAPACK's {s,d}lasd4 [CITATION] for this computational strategy
    temp = z(origin) / (state.dPlusShift(origin)*state.dMinusShift(origin));
    state.secularDeriv = (state.psiMinusDeriv + temp*temp) + state.phiDeriv;
    temp *= z(origin);
    state.secular = state.secularMinus + temp;
    state.relErrorBound +=
      8*(state.phi-state.psiMinus) + 2*rhoInv + 3*Abs(temp);
    if( penalizeDerivative )
        state.relErrorBound += Abs(state.rootRelEst)*state.secularDeriv;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void EvaluateSecularLast
( const Real& rho,
  const Matrix<Real>& z,
        LastState<Real>& state,
  bool penalizeDerivative=false )
{
    EL_DEBUG_CSE
    const Real one(1);
    const Int n = z.Height();
    const Int origin = n-1;
    const Real rhoInv = one / rho;
    Real temp;

    state.psiMinus = 0;
    state.psiMinusDeriv = 0;
    state.relErrorBound = 0;
    for( Int j=0; j<origin; ++j )
    {
        // Compute the j'th term divided by z(j)
        temp = z(j) / (state.dPlusShift(j)*state.dMinusShift(j));
        state.psiMinus += z(j)*temp;
        state.psiMinusDeriv += temp*temp;
        state.relErrorBound += state.psiMinus; // This should be negative
    }
    state.relErrorBound = Abs(state.relErrorBound); // This should be negation

    // Compute the origin term divided by z(origin)
    temp = z(origin) / (state.dPlusShift(origin)*state.dMinusShift(origin));
    state.psiOrigin = z(origin)*temp;
    state.psiOriginDeriv = temp*temp;
    state.secularDeriv = state.psiMinusDeriv + state.psiOriginDeriv;
    state.relErrorBound -= state.psiOrigin;

    // Compute the secular function with the origin term removed
    // (and its derivative)
    const Real psi = state.psiMinus + state.psiOrigin;
    state.secular = rhoInv + psi;

    state.relErrorBound += rhoInv - 8*psi;
    if( penalizeDerivative )
        state.relErrorBound += Abs(state.rootRelEst)*state.secularDeriv;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularInitialGuess
( const Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        State<Real>& state,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0), one(1);
    const Int k = whichValue;
    const Int n = d.Height();
    const Real diagSqDiffHalf = state.diagSqDiff / 2;
    if( k >= n-1 )
        LogicError("Assumption of inner singular value broken");

    // We can test whether the k'th eigenvalue is in the first or second
    // half of (d(k)^2, d(k+1)^2) by testing the sign of the secular
    // equation evaluated at the square-root of the center point,
    //
    //   center = (d(k)^2 + d(k+1)^2)/2.
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
    //   f(x) = (1/rho) + sum_{j=0}^{n-1} z(j)^2 / (d(j)^2 - x)
    //        = (1/rho) + psi_k(x) + phi_k(x),
    //
    // where
    //
    //   psi_k(x) = sum_{j=0}^k       z(j)^2 / (d(j)^2 - x),
    //   phi_k(x) = sum_{j=k+1}^{n-1} z(j)^2 / (d(j)^2 - x),
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
          z(j)*z(j) / (state.dPlusShift(j)*state.dMinusShift(j));
    Real phiMinus = zero; 
    for( Int j=k+2; j<n; ++j )
        phiMinus +=
          z(j)*z(j) / (state.dPlusShift(j)*state.dMinusShift(j));

    // The secular equation can now be expressed as
    //
    //   f(center) = (1/rho) + psiMinus + phiMinus +
    //       z(k)^2 / (d(k)^2 - center) + z(k+1)^2 / (d(k+1)^2 - center),
    //
    // where the top row will turn out to be the 'a' coefficient of a 
    // quadratic equation a x^2 + b x + c = 0 that we will solve to compute
    // the initial guess for the sought-after eigenvalue. But recall that
    // we should carefully compute (d(j)^2 - center).
    const Real a = one/rho + state.psiMinus + phiMinus;
    const Real secularCenter = a +
      z(k)*z(k) / (state.dPlusShift(k)*state.dMinusShift(k)) +
      z(k+1)*z(k+1) / (state.dPlusShift(k+1)*state.dMinusShift(k+1));

    state.geometricAverage = false; // TODO(poulson): Documentation?
    if( secularCenter >= zero )
    {
        // The eigenvalue lives in (d(k)^2, center), so we solve the 
        // quadratic equation
        //
        //   a + z(k)^2 / (d(k)^2 - x) + z(k+1)^2 / (d(k+1)^2 - x) = 0,
        //
        // directly for rootRelEst = x - d(k)^2, which yields the quadratic
        // equation a rootRelEst^2 + b rootRelEst + c = 0, with
        //
        //   b = -a gap - (z(k)^2 + z(k+1)^2),
        //   c = z(k)^2 gap,
        //
        // with gap = d(k+1)^2 - d(k)^2.
        state.originOnLeft = true;
        state.origin = k;
        state.sigmaRelLowerBound = zero; 
        state.sigmaRelUpperBound = shiftedCenterRoot;
        const Real bNeg = a*state.diagSqDiff + z(k)*z(k) + z(k+1)*z(k+1);
        const Real c = z(k)*z(k)*state.diagSqDiff;

        // Compute rootRelEst = sigmaEst^2 - d(k)^2
        state.rootRelEst = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix ); 
        state.sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue
          ( state.rootRelEst, d(state.origin) );

        // TODO(poulson): Document this LAPACK heuristic
        const Real eps = limits::Epsilon<Real>();
        const Real sqrtEps = Sqrt(eps);
        if( (d(k) <= sqrtEps*d(k+1)) &&
            (Abs(z(k)) <= sqrtEps) &&
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
        //   a + z(k)^2 / (d(k)^2 - x) + z(k+1)^2 / (d(k+1)^2 - x) = 0,
        //
        // directly for rootRelEst = x - d(k+1)^2, which yields the quadratic
        // equation a rootRelEst^2 + b rootRelEst + c = 0, with
        //
        //   b = a gap - (z(k)^2 + z(k+1)^2),
        //   c = -z(k)^2 gap,
        //
        // with gap = d(k+1)^2 - d(k)^2.
        state.originOnLeft = false;
        state.origin = k+1;

        // Again follow LAPACK's {s,d}lasd4's [CITATION] lead in carefully
        // computing centerRoot - d(k+1).
        state.sigmaRelLowerBound =
          -diagSqDiffHalf / (d(state.origin) + centerRoot);
        state.sigmaRelUpperBound = zero; 
        const Real bNeg = -a*state.diagSqDiff + z(k)*z(k) + z(k+1)*z(k+1);
        const Real c = -z(k+1)*z(k+1)*state.diagSqDiff;

        state.rootRelEst = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix );
        state.sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue
          ( state.rootRelEst, d(state.origin) );
    }
    if( ctrl.progress )
        Output
        ("Initial relative interval is [",state.sigmaRelLowerBound,",",
         state.sigmaRelUpperBound,"]");

    state.sigmaEst = state.sigmaRelEst + d(state.origin);
    for( Int j=0; j<n; ++j ) 
    {
        state.dPlusShift(j) = (d(j) + d(state.origin)) + state.sigmaRelEst;
        state.dMinusShift(j) = (d(j) - d(state.origin)) - state.sigmaRelEst;
    }

    EvaluateSecular( rho, z, state, ctrl.penalizeDerivative );
}

// For seeking the last root of the secular equation in
//
//   (d(n-1)^2, d(n-1)^2+rho).
//
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularInitialGuessLast
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        LastState<Real>& state,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = d.Height();
    const Real zero(0), one(1);
    const Real rhoInv = one / rho;

    const Int origin = n-1;

    // Since the largest eigenvalue must live within
    //
    //   d(n-1)^2 + (0, rho ||z||_2^2) = d(n-1)^2 + (0, rho),
    //
    // we use the center as our initial guess. Therefore, put
    // 
    //   rootRelEst = (rho ||z||_2^2) / 2 = rho / 2 = sigmaEst^2 - d(n-1)^2.
    //
    const Real shiftedCenter = rho / 2;
    const Real shiftedCenterRoot =
      RelativeEigenvalueToRelativeSingularValue( shiftedCenter, d(origin) );

    state.rootRelEst = shiftedCenter;
    state.sigmaRelEst = shiftedCenterRoot;
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
          z(j)*z(j) / (state.dPlusShift(j)*state.dMinusShift(j));

    state.psiMinus = psiDoubleMinus +
      z(origin-1)*z(origin-1) /
      (state.dPlusShift(origin-1)*state.dMinusShift(origin-1));

    state.psiOrigin = z(origin)*z(origin) /
      (state.dPlusShift(origin)*state.dMinusShift(origin));

    const Real a = rhoInv + psiDoubleMinus;
    state.secularMinus = rhoInv + state.psiMinus;
    state.secular = state.secularMinus + state.psiOrigin;

    if( state.secular <= zero )
    {
        // Evaluate the negative of the two terms of the secular equation
        // near the origin at d(origin)^2 + rho to see if they would overpower
        // the remaining terms. If so, we accept d(origin)^2 + rho as the
        // initial guess, otherwise we solve the implied quadratic equation.
        //
        // This strategy originated in LAPACK's {s,d}lasd4 [CITATION].
        //
        const Real singValUpperBound = Sqrt(d(origin)*d(origin) + rho);
        const Real temp =
          z(origin-1)*z(origin-1) /
           ((d(origin-1) + singValUpperBound)*
            (d(origin) - d(origin-1) + rho/(d(origin) + singValUpperBound))) +
          z(origin)*z(origin) / rho;

        if( a <= temp )
        {
            state.rootRelEst = rho;
        }
        else
        {
            const Real bNeg = -a*state.diagSqDiff + z(origin-1)*z(origin-1) +
              z(origin)*z(origin);
            const Real c = -z(origin)*z(origin)*state.diagSqDiff;
            state.rootRelEst =
              SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
            state.sigmaRelEst =
              RelativeEigenvalueToRelativeSingularValue
              ( state.rootRelEst, d(origin) );
        }
        state.sigmaRelLowerBound = shiftedCenterRoot;
        state.sigmaRelUpperBound =
          RelativeEigenvalueToRelativeSingularValue( rho, d(origin) );
    }
    else
    {
        // We will approximate the secular equation with
        //
        //  f_{n-1,n-2}(y) + z(n-2)^2 / (d(n-2)^2 - x) +
        //                   z(n-1)^2 / (d(n-1)^2 - x),
        //
        // where f_{n-1,n-2}(y) is the secular equation, with the (n-1)'th and 
        // (n-2)'th terms removed, evaluated at the current estimate. We solve
        // for a root of this equation in terms of rotRelEst = x - d(n-1)^2.
        // In particular, we pick the '-' branch of the quadratic equation.

        const Real bNeg = -a*state.diagSqDiff + z(origin-1)*z(origin-1) +
          z(origin)*z(origin);
        const Real c = -z(origin)*z(origin)*state.diagSqDiff;
        state.rootRelEst = SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
        state.sigmaRelEst =
          RelativeEigenvalueToRelativeSingularValue
          ( state.rootRelEst, d(origin) );
        state.sigmaRelLowerBound = 0;
        state.sigmaRelUpperBound = shiftedCenterRoot;
    }
    if( ctrl.progress )
        Output
        ("Initial relative interval is [",state.sigmaRelLowerBound,",",
         state.sigmaRelUpperBound,"]");

    state.sigmaEst = state.sigmaRelEst + d(origin);
    for( Int j=0; j<n; ++j ) 
    {
        state.dPlusShift(j) = (d(j) + d(origin)) + state.sigmaRelEst;
        state.dMinusShift(j) = (d(j) - d(origin)) - state.sigmaRelEst;
    }
    EvaluateSecularLast( rho, z, state, ctrl.penalizeDerivative );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularUpdate
( bool  initialize,
  Int   whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        State<Real>& state,
        SecularSVDInfo& info,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Int n = d.Height();
    const Int k = whichValue;
    const Int origin = state.origin;
    Real eta;

    if( state.secular <= zero )
        state.sigmaRelLowerBound =
          Max( state.sigmaRelLowerBound, state.sigmaRelEst );
    else
        state.sigmaRelUpperBound =
          Min( state.sigmaRelUpperBound, state.sigmaRelEst );
    if( ctrl.progress )
        Output
        ("Relative interval is [",state.sigmaRelLowerBound,",",
         state.sigmaRelUpperBound,"]");

    if( state.useThreePoles )
    {
        // Use the "Hybrid Scheme" described in subsection 3.4 of LAWN 89
        // [CITATION] and implemented within {s,d}lasd4 [CITATION].
 
        // Carefully compute
        //
        //   leftGap = d(origin-1)^2 - sigmaEst^2, and
        //   rightGap = d(origin+1)^2 - sigmaEst^2.
        //
        const Real leftGap =
          state.dPlusShift(origin-1)*state.dMinusShift(origin-1);
        const Real rightGap =
          state.dPlusShift(origin+1)*state.dMinusShift(origin+1);

        Real a;
        Matrix<Real> zCubic(3,1);
        if( state.alternateStrategy )
        {
            a = state.secularMinus - leftGap*state.psiMinusDeriv -
              rightGap*state.phiDeriv;
            zCubic(0) = leftGap*leftGap*state.psiMinusDeriv;
            zCubic(1) = z(origin)*z(origin);
            zCubic(2) = rightGap*rightGap*state.phiDeriv;
        }
        else
        {
            // Carefully compute doubleGap = d(origin+1)^2 - d(origin-1)^2
            const Real doubleGap =
              (d(origin+1)-d(origin-1))*(d(origin+1)+d(origin-1));

            if( state.originOnLeft )
            {
                // Since the shift origin, m, is k, We will interpolate the
                // secular equation as
                //
                //  Q(x; a, s, S) = a + z(m-1)^2 / (d(m-1)^2 - x) + 
                //                      z(m  )^2 / (d(m  )^2 - x) + 
                //                      S        / (d(m+1)^2 - x),
                //
                // which can be effected by applying the Middle Way, the Fixed
                // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme
                //  to
                //
                //   f_m(y) = f(y) - z(m)^2 / (d(m)^2 - x).
                //
                // We balance between the Fixed Weight Method and the Middle
                // Way; their formulae for computing a are identical. In both
                // cases,
                // 
                //  a = f_m - rightGap*f'_m +
                //      (z(m-1)/leftGap)^2*(d(m+1)^2-d(m-1)^2).
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
                const Real leftRatio = z(origin-1) / leftGap;
                const Real leftDerivTerm = leftRatio*leftRatio;
                a =
                  (state.secularMinus-rightGap*state.secularMinusDeriv) +
                  leftDerivTerm*doubleGap;

                // We could either flip or clip here, but LAPACK tends to flip
                // and clips here, so we will follow suit so that requesting
                // flips leads to mirroring LAPACK.
                const Real psiDoubleMinusDeriv =
                  Max( state.psiMinusDeriv-leftDerivTerm, zero );
    
                zCubic(0) = z(origin-1)*z(origin-1);
                zCubic(1) = z(origin)*z(origin);
                zCubic(2) =
                  rightGap*rightGap*(psiDoubleMinusDeriv+state.phiDeriv);
            }
            else
            {
                // Since the shift origin, m, is k+1, we will interpolate the
                // secular equation as
                //
                //  Q(x; a, s, S) = a + s        / (d(m-1)^2 - x) +
                //                      z(m)^2   / (d(m  )^2 - x) +
                //                      z(m+1)^2 / (d(m+1)^2 - x),
                //
                // which can be effected by applying the Middle Way, the Fixed
                // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme
                //  to
                //
                //   f_m(y) = f(y) - z(m)^2 / (d(m)^2 - x).
                //
                // We balance between the Fixed Weight Method and the Middle
                // Way; their formulae for computing a are identical. In both
                // cases,
                // 
                //  a = f_m - leftGap*f'_m -
                //      (z(m+1)/rightGap)^2*(d(m+1)^2-d(m-1)^2),
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
                const Real rightRatio = z(origin+1) / rightGap;
                const Real mPlusDerivTerm = rightRatio*rightRatio;
                a =
                  (state.secularMinus-leftGap*state.secularMinusDeriv) -
                  mPlusDerivTerm*doubleGap;

                // We could either flip or clip here, but LAPACK tends to flip
                // and clips here, so we will follow suit so that requesting
                // flips leads to mirroring LAPACK.
                const Real phiMinusDeriv =
                  Max( state.phiDeriv-mPlusDerivTerm, zero );
            
                zCubic(0) = leftGap*leftGap*(state.psiMinusDeriv+phiMinusDeriv);
                zCubic(1) = z(origin)*z(origin);
                zCubic(2) = z(origin+1)*z(origin+1);
            }
        }
        Matrix<Real> dCubic(3,1);
        dCubic(0) = leftGap;
        dCubic(1) = state.dPlusShift(origin)*state.dMinusShift(origin);
        dCubic(2) = rightGap;

        auto cubicInfo =
          CubicSecular
          ( initialize, state.originOnLeft, a, zCubic, dCubic, state.secular,
            eta, ctrl.cubicCtrl );
        info.numCubicIterations += cubicInfo.numIterations;
        if( cubicInfo.converged )
        {
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
                const Real temp = z(k) / kGap;
                a = state.secular - kp1Gap*state.secularDeriv +
                  state.diagSqDiff*(temp*temp);
            }
            else
            {
                const Real temp = z(k+1) / kp1Gap;
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
                    bNeg = z(k)*z(k) + kp1Gap*kp1Gap*state.secularMinusDeriv;
                }
                else
                {
                    bNeg = z(k+1)*z(k+1) + kGap*kGap*state.secularMinusDeriv; 
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
                const Real temp = z(k) / kGap;
                a = state.secular - kGap*(state.psiMinusDeriv+temp*temp) -
                  kp1Gap*state.phiDeriv;
            }
            else
            {
                const Real temp = z(k+1) / kp1Gap;
                a = state.secular - kGap*state.psiMinusDeriv -
                  kp1Gap*(state.phiDeriv+temp*temp);
            }
        }
        else
        {
            if( state.originOnLeft )
            {
                // Apply Eq'ns (32) through (34) from LAWN 89 [CITATION] to set
                //
                //  a = f - kp1Gap f' + (z(k)^2/kGap^2) (d(k+1)^2 - d(k)^2),
                //
                // and implicitly set
                //
                //  s = z(k)^2,
                //  S = kp1Gap^2 (f' - z(k)^2/kGap^2),
                //
                // where y is the current estimate of the k'th eigenvalue and 
                // f and f' are the secular equation and its derivative
                // evaluated at y.
                const Real temp = z(k) / kGap;
                a = state.secular - kp1Gap*state.secularDeriv +
                  (temp*temp)*state.diagSqDiff;
            }
            else
            {
                // Apply Eq'ns (35) through (37) from LAWN 89 [CITATION] to set
                //
                //  a = f - kGap f' - (z(k+1)^2/kp1Gap^2) (d(k+1)^2 - d(k)^2),
                //
                // and implicitly set
                //
                //  s = kGap^2 (f' - z(k+1)^2/kp1Gap^2),
                //  S = z(k+1)^2.
                //
                const Real temp = z(k+1) / kp1Gap;
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
                    bNeg = z(k)*z(k) + kGap*kGap*state.psiMinusDeriv +
                           kp1Gap*kp1Gap*state.phiDeriv;
                }
                else
                {
                    bNeg = z(k+1)*z(k+1) + kGap*kGap*state.psiMinusDeriv +
                           kp1Gap*kp1Gap*state.phiDeriv;
                }
            }
            else
            {
                if( state.originOnLeft )
                    bNeg = z(k)*z(k) + kp1Gap*kp1Gap*state.secularMinusDeriv;
                else
                    bNeg = z(k+1)*z(k+1) + kGap*kGap*state.secularMinusDeriv;
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
    state.rootRelEst = -state.dPlusShift(origin)*state.dMinusShift(origin);
    EvaluateSecular( rho, z, state, ctrl.penalizeDerivative );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularUpdateLast
( bool  initialize,
  const Real& rho,
  const Matrix<Real>& z,
        LastState<Real>& state,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Int n = state.dPlusShift.Height();
    const Int origin = n-1;

    if( state.secular <= zero )
        state.sigmaRelLowerBound =
          Max( state.sigmaRelLowerBound, state.sigmaRelEst );
    else
        state.sigmaRelUpperBound =
          Min( state.sigmaRelUpperBound, state.sigmaRelEst );
    if( ctrl.progress )
        Output
        ("Relative interval is [",state.sigmaRelLowerBound,",",
         state.sigmaRelUpperBound,"]");

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
            //
            // Note that LAPACK's {s,d}lasd4 [CITATION] uses the equivalent of
            // the following, but LAPACK's {s,d}laed4 [CITATION] uses the 
            // best known upper bound so far rather than rho.
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
        // While LAPACK's {s,d}lasd4 uses this approach for the inner singular
        // values, it does not for outer singular values and instead uses
        // trivial lower and upper bounds. We part ways with LAPACK here.

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
    }

    state.sigmaEst += eta;
    state.sigmaRelEst += eta;
    state.secularOld = state.secular;
    for( Int j=0; j<n; ++j )
    {
        state.dPlusShift(j) += eta;
        state.dMinusShift(j) -= eta;
    }
    state.rootRelEst = -state.dPlusShift(origin)*state.dMinusShift(origin);
    EvaluateSecularLast( rho, z, state, ctrl.penalizeDerivative );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularSVDInfo
SecularInner
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        State<Real>& state,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
      if( k == n-1 )
          LogicError("SecularInner meant for inner singular values");
      if( n <= 2 )
          LogicError("SecularInner meant for n > 2");
    )

    SecularSVDInfo info;
    state.dMinusShift.Resize(n,1);
    state.dPlusShift.Resize(n,1);

    const Real diagSum = d(k+1) + d(k);
    const Real diagDiff = d(k+1) - d(k);
    state.diagSqDiff = diagSum*diagDiff; // d(k+1)^2 - d(k)^2

    SecularInitialGuess( k, d, rho, z, state, ctrl );
    ++info.numIterations;

    // See the discussion in LAWN 89 [CITATION] on the usage of three poles
    state.useThreePoles =
      (state.originOnLeft ? state.secularMinus<zero : state.secularMinus>zero);
    if( state.origin == 0 || state.origin == n-1 )
        state.useThreePoles = false;

    // Check if we have already converged
    if( Abs(state.secular) <= eps*state.relErrorBound )
        return info;

    // Compute the first update to our estimate
    bool initialize = true;
    state.alternateStrategy = false;
    SecularUpdate( initialize, k, d, rho, z, state, info, ctrl );
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
        EL_DEBUG_ONLY(
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

        // Decide the next step
        SecularUpdate( initialize, k, d, rho, z, state, info, ctrl );
        ++info.numIterations;
   
        // This strategy was described in Ren-Cang Li's LAWN 89
        if( state.secular*state.secularOld > zero &&
            Abs(state.secular) > Abs(state.secularOld)*ctrl.sufficientDecay )
        {
            state.alternateStrategy = !state.alternateStrategy;
            ++info.numAlternations;
        }
    }
    return info;
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularSVDInfo
SecularLast
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        LastState<Real>& state,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
      if( k != n-1 )
          LogicError("SecularLast meant for largest singular value");
      if( n <= 2 )
          LogicError("SecularLast meant for n > 2");
    )

    SecularSVDInfo info;
    state.dMinusShift.Resize(n,1);
    state.dPlusShift.Resize(n,1);
    const Real diagSum = d(n-1) + d(n-2);
    const Real diagDiff = d(n-1) - d(n-2);
    state.diagSqDiff = diagSum*diagDiff; // d(n-1)^2 - d(n-2)^2

    SecularInitialGuessLast( d, rho, z, state, ctrl );
    ++info.numIterations;

    if( Abs(state.secular) <= eps*state.relErrorBound )
        return info;

    // Calculate the first update
    bool initialize = true;
    SecularUpdateLast( initialize, rho, z, state, ctrl );
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
        SecularUpdateLast( initialize, rho, z, state, ctrl );
        ++info.numIterations;
    }
    return info;
}

} // namespace secular_svd

// Compute a single singular value corresponding to the square-root of the 
// eigenvalue of the diagonal plus rank one matrix
//
//     diag(d)^2 + rho z z^T,
//
// where || z ||_2 = 1, with
//
//     0 <= d(0) < d(1) < ... < d(n-1)
//
// and rho > 0.
//
// This routine loosely corresponds to LAPACK's {s,d}lasd4 [CITATION].
//

template<typename Real,typename>
SecularSVDInfo
SecularSingularValue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& singularValue,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
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
      if( z.Height() != n )
          LogicError("z was not the correct height");
      // TODO(poulson): Check the assumption that || z ||_2 = 1
    )

    SecularSVDInfo info;
    if( n == 1 )
    {
        singularValue = Sqrt(d(0)*d(0) + rho);
        return info;
    }
    else if( n == 2 )
    {
        singularValue = secular_svd::TwoByTwo( k, d(0), d(1), rho, z(0), z(1) );
        return info;
    }

    if( k < n-1 )
    {
        secular_svd::State<Real> state;
        info = secular_svd::SecularInner( k, d, rho, z, state, ctrl );
        singularValue = state.sigmaEst;
    }
    else
    {
        secular_svd::LastState<Real> state;
        info = secular_svd::SecularLast( k, d, rho, z, state, ctrl );
        singularValue = state.sigmaEst;
    }

    return info;
}

template<typename Real,typename>
SecularSVDInfo
SecularSingularValue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& singularValue,
        Matrix<Real>& dMinusShift,
        Matrix<Real>& dPlusShift,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
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
      if( z.Height() != n )
          LogicError("z was not the correct height");
      // TODO(poulson): Check the assumption that || z ||_2 = 1
    )

    SecularSVDInfo info;
    dMinusShift.Resize(n,1);
    dPlusShift.Resize(n,1);
    if( n == 1 )
    {
        singularValue = Sqrt(d(0)*d(0) + rho);
        // TODO(poulson): Make this computation more accurate for completeness?
        dMinusShift(0) = d(0) - singularValue;
        dPlusShift(0) = d(0) + singularValue;
        return info;
    }
    else if( n == 2 )
    {
        singularValue =
          secular_svd::TwoByTwo
          ( k, d(0), d(1), rho, z(0), z(1),
            dMinusShift(0), dMinusShift(1),
            dPlusShift(0), dPlusShift(1) );
        return info;
    }

    if( k < n-1 )
    {
        secular_svd::State<Real> state;
        info = secular_svd::SecularInner( k, d, rho, z, state, ctrl );
        singularValue = state.sigmaEst;
        dPlusShift = state.dPlusShift;
        dMinusShift = state.dMinusShift;
    }
    else
    {
        secular_svd::LastState<Real> state;
        info = secular_svd::SecularLast( k, d, rho, z, state, ctrl );
        singularValue = state.sigmaEst;
        dPlusShift = state.dPlusShift;
        dMinusShift = state.dMinusShift;
    }

    return info;
}

template<typename Real,typename>
SecularSVDInfo
SecularSVD
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Matrix<Real>& U,
        Matrix<Real>& s,
        Matrix<Real>& V,
  const SecularSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = d.Height();
    SecularSVDInfo info;

    s.Resize( n, 1 );
    if( n == 0 )
    {
        U.Resize( n, n );
        V.Resize( n, n );
        return info;
    }

    // TODO(poulson): Batch secular equation solvers?
    //
    // Compute all of the singular values and the vector r ~= sqrt(rho) z which
    // would produce the given singular values to high relative accuracy.
    //
    // We accumulate the products involved in computing r in an online manner
    // and take the square-root of the absolute value (and the correct sign)
    // at the end. The key is to recognize that the only term left out of entry
    // i of the corrected vector in Eq. (3.6) of Gu/Eisenstat in the product
    //
    //    prod_{k=0}^{n-1} (sigma_k^2 - d(i)^2) / (d(k)^2 - d(i)^2)
    //
    // is 1 / (d(i)^2 - d(i)^2) (and we emphasize that the numerator is kept).
    // It is thus easy to see how the first term and second product in
    //
    //    r(i)^2 = (sigma_{n-1}^2 - d(i)^2) *
    //      prod_{k=0}^{i-1} (sigma_k^2 - d(i)^2) / (d(k  )^2 - d(i)^2) *
    //      prod_{k=i}^{n-2} (sigma_k^2 - d(i)^2) / (d(k+1)^2 - d(i)^2)
    //
    // can be reorganized to yield
    //
    //    r(i)^2 = (sigma_{i}^2 - d(i)^2) *
    //      prod_{k=0  }^{i-1} (sigma_k^2 - d(i)^2) / (d(k)^2 - d(i)^2) *
    //      prod_{k=i+1}^{n-1} (sigma_k^2 - d(i)^2) / (d(k)^2 - d(i)^2).
    //
    // Given that we have available only one sigma_j at a time, we greedily
    // update all r(i) products given each sigma_j
    // (Cf. LAPACK's {s,d}lasd8 [CITATION] for this approach).
    //
    U.Resize( n, n );
    V.Resize( n, n );
    Matrix<Real> r;
    Ones( r, n, 1 );
    auto vScratch = V(ALL,IR(0));
    for( Int j=0; j<n; ++j )
    {
        // While we will temporarily store dMinusShift and dPlusShift in the
        // j'th columns of U and V, respectively, it is worth noting that we 
        // only require access to their Hadamard product after this loop ends.
        auto u = U(ALL,IR(j));
        auto valueInfo =
          SecularSingularValue( j, d, rho, z, s(j), u, vScratch, ctrl );
        
        info.numIterations += valueInfo.numIterations;
        info.numAlternations += valueInfo.numAlternations;
        info.numCubicIterations += valueInfo.numCubicIterations;
        info.numCubicFailures += valueInfo.numCubicFailures;

        // u currently hold d-s(j) and vScratch currently holds d+s(j).
        // Overwrite u with their element-wise product since that is all we 
        // require from here on out.
        for( Int k=0; k<n; ++k )
            u(k) *= vScratch(k);
      
        r(j) *= u(j);
        for( Int k=0; k<n; ++k )
        {
            if( k == j )
                continue;
            r(k) *= u(k) / ((d(j)+d(k))*(d(j)-d(k)));
        }
    }
    for( Int j=0; j<n; ++j )
        r(j) = Sgn(z(j),false) * Sqrt(Abs(r(j)));

    for( Int j=0; j<n; ++j )
    {
        // Compute the j'th left and right singular vectors via
        // Eqs. (3.4) and (3.3), respectively.
        auto u = U(ALL,IR(j));
        auto v = V(ALL,IR(j));
        {
            const Real deltaSqMinusShiftSq = u(0);
            u(0) = -1;
            v(0) = r(0) / deltaSqMinusShiftSq;
        }
        for( Int i=1; i<n; ++i )
        {
            const Real deltaSqMinusShiftSq = u(i);
            v(i) = r(i) / deltaSqMinusShiftSq;
            u(i) = d(i) * v(i);
        }
        u *= Real(1) / FrobeniusNorm( u );
        v *= Real(1) / FrobeniusNorm( v );
    }

    return info;
}

#define PROTO(Real) \
  template SecularSVDInfo \
  SecularSingularValue \
  ( Int whichValue, \
    const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& z, \
          Real& singularValue, \
    const SecularSVDCtrl<Real>& ctrl ); \
  template SecularSVDInfo \
  SecularSingularValue \
  ( Int whichValue, \
    const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& z, \
          Real& singularValue, \
          Matrix<Real>& dMinusShift, \
          Matrix<Real>& dPlusShift, \
    const SecularSVDCtrl<Real>& ctrl ); \
  template SecularSVDInfo \
  SecularSVD \
  ( const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& z, \
          Matrix<Real>& U, \
          Matrix<Real>& s, \
          Matrix<Real>& V, \
    const SecularSVDCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT

#include <El/macros/Instantiate.h>

} // namespace El
