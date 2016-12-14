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

#include "./SecularEVD/TwoByTwo.hpp"

namespace El {

namespace secular_evd {

template<typename Real>
struct State
{
    Int origin;
    bool originOnLeft;

    Matrix<Real> dMinusShift;

    Real diagDiff;

    Real rootEst; 
    Real rootRelEst;
    Real rootRelLowerBound;
    Real rootRelUpperBound;

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
    bool alternateStrategy;
};

template<typename Real>
struct LastState
{
    Matrix<Real> dMinusShift;

    Real diagDiff;

    Real rootEst; 
    Real rootRelEst;
    Real rootRelLowerBound;
    Real rootRelUpperBound;

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
  bool penalizeDerivative=true )
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
    // {s,d}laed4 [CITATION]). The loop direction is chosen to heuristically
    // sum from small to large components.
    state.psiMinus = 0;
    state.psiMinusDeriv = 0;
    state.relErrorBound = 0;
    for( Int j=0; j<origin; ++j )
    {
        // Compute the j'th term divided by z(j)
        temp = z(j) / state.dMinusShift(j);
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
        temp = z(j) / state.dMinusShift(j);
        state.phi += z(j)*temp;
        state.phiDeriv += temp*temp;
        state.relErrorBound += state.phi;
    }

    // Compute the secular function with the origin term removed
    // (and its derivative)
    state.secularMinus = rhoInv + state.psiMinus + state.phi;
    state.secularMinusDeriv = state.psiMinusDeriv + state.phiDeriv;

    // Cf. LAPACK's {s,d}laed4 [CITATION] for this computational strategy
    temp = z(origin) / state.dMinusShift(origin);
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
  bool penalizeDerivative=true )
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
        temp = z(j) / state.dMinusShift(j);
        state.psiMinus += z(j)*temp;
        state.psiMinusDeriv += temp*temp;
        state.relErrorBound += state.psiMinus; // This should be negative
    }
    state.relErrorBound = Abs(state.relErrorBound); // This should be negation

    // Compute the origin term divided by z(origin)
    temp = z(origin) / state.dMinusShift(origin);
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
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0), one(1);
    const Int k = whichValue;
    const Int n = d.Height();
    if( k >= n-1 )
        LogicError("Assumption of inner eigenvalue broken");

    // We can test whether the k'th eigenvalue is in the first or second
    // half of (d(k) d(k+1)) by testing the sign of the secular
    // equation evaluated at the center point,
    //
    //   center = (d(k) + d(k+1))/2.
    //
    // Follow LAPACK's {s,d}laed4's [CITATION] lead in carefully computing
    //
    //   shiftedCenter = (d(k) + d(k+1))/2 - d(k) = (d(k+1) - d(k))/2.
    //
    const Real shiftedCenter = state.diagDiff / 2;

    // Again use LAPACK's {s,d}laed4's [CITATION] suggestion that the
    // diagonal entries minus the center roots be computed in a safe
    // (but somewhat obscured) manner.
    for( Int j=0; j<n; ++j ) 
        state.dMinusShift(j) = (d(j) - d(k)) - shiftedCenter;

    // Given the partition of the secular equation as
    //
    //   f(x) = (1/rho) + sum_{j=0}^{n-1} z(j)^2 / (d(j) - x)
    //        = (1/rho) + psi_k(x) + phi_k(x),
    //
    // where
    //
    //   psi_k(x) = sum_{j=0}^k       z(j)^2 / (d(j) - x),
    //   phi_k(x) = sum_{j=k+1}^{n-1} z(j)^2 / (d(j) - x),
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
        state.psiMinus += z(j)*z(j) / state.dMinusShift(j);
    Real phiMinus = zero; 
    for( Int j=k+2; j<n; ++j )
        phiMinus += z(j)*z(j) / state.dMinusShift(j);

    // The secular equation can now be expressed as
    //
    //   f(center) = (1/rho) + psiMinus + phiMinus +
    //       z(k)^2 / (d(k) - center) + z(k+1)^2 / (d(k+1) - center),
    //
    // where the top row will turn out to be the 'a' coefficient of a 
    // quadratic equation a x^2 + b x + c = 0 that we will solve to compute
    // the initial guess for the sought-after eigenvalue. But recall that
    // we should carefully compute (d(j) - center).
    const Real a = one/rho + state.psiMinus + phiMinus;
    const Real secularCenter = a +
      z(k)*z(k) / state.dMinusShift(k) + z(k+1)*z(k+1) / state.dMinusShift(k+1);

    if( secularCenter >= zero )
    {
        // The eigenvalue lives in (d(k), center), so we solve the 
        // quadratic equation
        //
        //   a + z(k)^2 / (d(k) - x) + z(k+1)^2 / (d(k+1) - x) = 0,
        //
        // directly for rootRelEst = x - d(k), which yields the quadratic
        // equation a rootRelEst^2 + b rootRelEst + c = 0, with
        //
        //   b = -a gap - (z(k)^2 + z(k+1)^2),
        //   c = z(k)^2 gap,
        //
        // with gap = d(k+1) - d(k).
        state.originOnLeft = true;
        state.origin = k;
        state.rootRelLowerBound = zero; 
        state.rootRelUpperBound = shiftedCenter;
        const Real bNeg = a*state.diagDiff + z(k)*z(k) + z(k+1)*z(k+1);
        const Real c = z(k)*z(k)*state.diagDiff;
        state.rootRelEst = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix ); 
    }
    else
    {
        // The eigenvalue lives in [center, d(k+1)), so we solve the
        // quadratic equation
        //
        //   a + z(k)^2 / (d(k) - x) + z(k+1)^2 / (d(k+1) - x) = 0,
        //
        // directly for rootRelEst = x - d(k+1), which yields the quadratic
        // equation a rootRelEst^2 + b rootRelEst + c = 0, with
        //
        //   b = a gap - (z(k)^2 + z(k+1)^2),
        //   c = -z(k)^2 gap,
        //
        // with gap = d(k+1) - d(k).
        state.originOnLeft = false;
        state.origin = k+1;
        state.rootRelLowerBound = -shiftedCenter;
        state.rootRelUpperBound = 0; 
        const Real bNeg = -a*state.diagDiff + z(k)*z(k) + z(k+1)*z(k+1);
        const Real c = -z(k+1)*z(k+1)*state.diagDiff;
        state.rootRelEst = SolveQuadraticMinus( a, bNeg, c, ctrl.negativeFix );
    }
    if( ctrl.progress )
        Output
        ("Initial relative interval is [",state.rootRelLowerBound,",",
         state.rootRelUpperBound,"]");

    state.rootEst = state.rootRelEst + d(state.origin);
    for( Int j=0; j<n; ++j ) 
        state.dMinusShift(j) = (d(j) - d(state.origin)) - state.rootRelEst;

    EvaluateSecular( rho, z, state, ctrl.penalizeDerivative );
}

// For seeking the last root of the secular equation in
//
//   (d(n-1), d(n-1)+rho).
//
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularInitialGuessLast
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        LastState<Real>& state,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = d.Height();
    const Real zero(0), one(1);
    const Real rhoInv = one / rho;

    const Int origin = n-1;

    // Since the largest eigenvalue must live within
    //
    //   d(n-1) + (0, rho ||z||_2^2) = d(n-1) + (0, rho),
    //
    // we use the center as our initial guess. Therefore, put
    // 
    //   rootRelEst = (rho ||z||_2^2) / 2 = rho / 2 = rootEst - d(n-1).
    //
    const Real shiftedCenter = rho / 2;
    state.rootRelEst = shiftedCenter;
    state.rootEst = state.rootRelEst + d(origin);
    for( Int j=0; j<n; ++j ) 
    {
        state.dMinusShift(j) = (d(j) - d(origin)) - state.rootRelEst;
    }

    // Form psi_{origin-2}
    Real psiDoubleMinus = zero;
    psiDoubleMinus = zero;
    for( Int j=0; j<origin-1; ++j )
        psiDoubleMinus += z(j)*z(j) / state.dMinusShift(j);

    state.psiMinus = psiDoubleMinus +
      z(origin-1)*z(origin-1) / state.dMinusShift(origin-1);

    state.psiOrigin = z(origin)*z(origin) / state.dMinusShift(origin);

    const Real a = rhoInv + psiDoubleMinus;
    state.secularMinus = rhoInv + state.psiMinus;
    state.secular = state.secularMinus + state.psiOrigin;

    if( state.secular <= zero )
    {
        // Evaluate the negative of the two terms of the secular equation
        // near the origin at d(origin) + rho to see if they would overpower
        // the remaining terms. If so, we accept d(origin) + rho as the initial
        // guess, otherwise we solve the implied quadratic equation.
        //
        // This strategy originated in LAPACK's {s,d}laed4 [CITATION].
        //
        const Real temp =
          z(origin-1)*z(origin-1) / (d(origin) - d(origin-1) + rho) +
          z(origin)*z(origin) / rho;

        if( a <= temp )
        {
            state.rootRelEst = rho;
        }
        else
        {
            const Real bNeg = -a*state.diagDiff + z(origin-1)*z(origin-1) +
              z(origin)*z(origin);
            const Real c = -z(origin)*z(origin)*state.diagDiff;
            state.rootRelEst =
              SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
        }
        state.rootRelLowerBound = shiftedCenter;
        state.rootRelUpperBound = rho;
    }
    else
    {
        // We will approximate the secular equation with
        //
        //  f_{n-1,n-2}(y) + z(n-2)^2 / (d(n-2) - x) +
        //                   z(n-1)^2 / (d(n-1) - x),
        //
        // where f_{n-1,n-2}(y) is the secular equation, with the (n-1)'th and 
        // (n-2)'th terms removed, evaluated at the current estimate. We solve
        // for a root of this equation in terms of tau = x - d(n-1).
        // In particular, we pick the '-' branch of the quadratic equation.

        const Real bNeg = -a*state.diagDiff + z(origin-1)*z(origin-1) +
          z(origin)*z(origin);
        const Real c = -z(origin)*z(origin)*state.diagDiff;
        state.rootRelEst = SolveQuadraticPlus( a, bNeg, c, ctrl.negativeFix );
        state.rootRelLowerBound = zero;
        state.rootRelUpperBound = shiftedCenter;
    }
    if( ctrl.progress )
        Output
        ("Initial relative interval is [",state.rootRelLowerBound,",",
         state.rootRelUpperBound,"]");

    state.rootEst = state.rootRelEst + d(origin);
    for( Int j=0; j<n; ++j ) 
        state.dMinusShift(j) = (d(j) - d(origin)) - state.rootRelEst;
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
        SecularEVDInfo& info,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Int n = d.Height();
    const Int k = whichValue;
    const Int origin = state.origin;
    Real eta;

    if( state.secular <= zero )
        state.rootRelLowerBound =
          Max( state.rootRelLowerBound, state.rootRelEst );
    else
        state.rootRelUpperBound =
          Min( state.rootRelUpperBound, state.rootRelEst );
    if( ctrl.progress )
        Output
        ("Relative interval is [",state.rootRelLowerBound,",",
         state.rootRelUpperBound,"]");

    if( state.useThreePoles )
    {
        // Use the "Hybrid Scheme" described in subsection 3.4 of LAWN 89
        // [CITATION] and implemented within {s,d}laed4 [CITATION].
 
        // Carefully compute
        //
        //   leftGap = d(origin-1) - rootEst, and
        //   rightGap = d(origin+1) - rootEst.
        //
        const Real leftGap = state.dMinusShift(origin-1);
        const Real rightGap = state.dMinusShift(origin+1);

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
            const Real doubleGap = d(origin+1) - d(origin-1);

            if( state.originOnLeft )
            {
                // Since the shift origin, m, is k, We will interpolate the
                // secular equation as
                //
                //  Q(x; a, s, S) = a + z(m-1)^2 / (d(m-1) - x) + 
                //                      z(m  )^2 / (d(m  ) - x) + 
                //                      S        / (d(m+1) - x),
                //
                // which can be effected by applying the Middle Way, the Fixed
                // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme
                //  to
                //
                //   f_m(y) = f(y) - z(m)^2 / (d(m) - x).
                //
                // We balance between the Fixed Weight Method and the Middle
                // Way; their formulae for computing a are identical. In both
                // cases,
                // 
                //  a = f_m - rightGap*f'_m +
                //      (z(m-1)/leftGap)^2*(d(m+1)-d(m-1)).
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
                //  Q(x; a, s, S) = a + s        / (d(m-1) - x) +
                //                      z(m)^2   / (d(m  ) - x) +
                //                      z(m+1)^2 / (d(m+1) - x),
                //
                // which can be effected by applying the Middle Way, the Fixed
                // Weight Method, or Borges/Gragg/Thornton/Warner's cubic scheme
                //  to
                //
                //   f_m(y) = f(y) - z(m)^2 / (d(m) - x).
                //
                // We balance between the Fixed Weight Method and the Middle
                // Way; their formulae for computing a are identical. In both
                // cases,
                // 
                //  a = f_m - leftGap*f'_m -
                //      (z(m+1)/rightGap)^2*(d(m+1)-d(m-1)),
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
        dCubic(1) = state.dMinusShift(origin);
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
            // For some reason, LAPACK's {s,d}laed4 [CITATION] choose to fail
            // if the cubic iteration did not succeed, but {s,d}lasd4 [CITATION]
            // uses a recovery process analogous to what follows. We will part 
            // ways with LAPACK here and run the recovery process rather than
            // simply failing.
            if( ctrl.progress )
                Output("Cubic did *not* converge");
            ++info.numCubicFailures;
            state.useThreePoles = false;
            const Real kGap = state.dMinusShift(k);
            const Real kp1Gap = state.dMinusShift(k+1);
            Real a;
            if( state.originOnLeft )
            {
                const Real temp = z(k) / kGap;
                a = state.secular - kp1Gap*state.secularDeriv +
                  state.diagDiff*(temp*temp);
            }
            else
            {
                const Real temp = z(k+1) / kp1Gap;
                a = state.secular - kGap*state.secularDeriv +
                  state.diagDiff*(temp*temp);
            }
            Real bNeg = (kGap+kp1Gap)*state.secular -
              kGap*kp1Gap*state.secularDeriv;
            const Real c = kGap*kp1Gap*state.secular;
            if( a == zero && bNeg == zero )
            {
                // a x^2 + b x + c = 0 has collapsed to c = 0, which is
                // nonsense. Follow LAPACK's {s,d}laed4 [CITATION] in handling
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
        //   Q(x; a, s, S) = a + s / (d(k) - x) + S / (d(k+1) - x).
        //
        // As noted after Proposition 3 of LAWN 89, the iteration formula does
        // not depend upon either s or S!

        // Safely compute the distance from our current eigenvalue estimate to
        // the surrounding grid points, d(k)^2 and d(k+1)^2.
        const Real kGap = state.dMinusShift(k);
        const Real kp1Gap = state.dMinusShift(k+1);

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
                //  a = f - kp1Gap f' + (z(k)^2/kGap^2) (d(k+1) - d(k)),
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
                  (temp*temp)*state.diagDiff;
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
                  (temp*temp)*state.diagDiff;
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
            // LAPACK's {s,d}laed4 [CITATION] and use a mysterious patch-up.
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

    // Since rootRelEst = rootEst - d(origin), the following sets
    // rootRelEstProp = rootEstProp - d(origin).
    const Real rootRelEstProp = state.rootRelEst + eta;
    if( rootRelEstProp > state.rootRelUpperBound ||
        rootRelEstProp < state.rootRelLowerBound )
    {
        if( ctrl.progress )
            Output("Stepped out of bounds");
        // We again follow LAPACK's {s,d}laed4 for how to handle this breakdown

        if( state.secular < zero )
        {
            // Since the secular equation implies that the root is to the right
            // of the un-updated estimate, but our update pushed us out of
            // bounds, we default to the update pushing us to the center of the
            // feasible domain, [rootRelEst, rootRelUpperBound).
            eta = (state.rootRelUpperBound - state.rootRelEst) / 2;
        }
        else
        {
            // Since the secular equation implies that the root is to the left
            // of the un-updated estimate, but our update pushed us out of
            // bounds, we default to the update being in the center of the
            // feasible domain, (rootRelLowerBound,rootRelEst].
            eta = (state.rootRelLowerBound - state.rootRelEst) / 2;
        }
    }

    state.rootEst += eta;
    state.rootRelEst += eta;
    state.secularOld = state.secular;
    for( Int j=0; j<n; ++j )
        state.dMinusShift(j) -= eta;
    EvaluateSecular( rho, z, state, ctrl.penalizeDerivative );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void SecularUpdateLast
( bool  initialize,
  const Real& rho,
  const Matrix<Real>& z,
        LastState<Real>& state,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Int n = state.dMinusShift.Height();
    const Int origin = n-1;

    if( state.secular <= zero )
        state.rootRelLowerBound =
          Max( state.rootRelLowerBound, state.rootRelEst );
    else
        state.rootRelUpperBound =
          Min( state.rootRelUpperBound, state.rootRelEst );
    if( ctrl.progress )
        Output
        ("Relative interval is [",state.rootRelLowerBound,",",
         state.rootRelUpperBound,"]");

    const Real kGap = state.dMinusShift(origin);
    const Real km1Gap = state.dMinusShift(origin-1);

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
            // TODO(poulson): Explain this special case.
            //
            // Note that LAPACK's {s,d}lasd4 [CITATION] uses the equivalent of 
            //
            //     eta = rho - rootEst,
            //
            // while LAPACK's {s,d}laed4 [CITATION] has such an updated
            // commented out and
            // instead performs the following.
            eta = state.rootRelUpperBound - state.rootRelEst;
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

    const Real rootRelEstProp = state.rootRelEst + eta;
    if( rootRelEstProp > state.rootRelUpperBound ||
        rootRelEstProp < state.rootRelLowerBound )
    {
        if( ctrl.progress )
            Output("Stepped out of bounds");
        // Move halfway towards the bound instead of exceeding it. Note that 
        // LAPACK's {s,d}laed4 [CITATION] follows this strategy, but LAPACK's
        // {s,d}lasd4 uses trivial upper and lower bounds rather than
        // continually updating them.
        if( state.secular < zero )
            eta = (state.rootRelUpperBound - state.rootRelEst) / 2;
        else 
            eta = (state.rootRelLowerBound - state.rootRelEst) / 2;
    }

    state.rootEst += eta;
    state.rootRelEst += eta;
    state.secularOld = state.secular;
    for( Int j=0; j<n; ++j )
        state.dMinusShift(j) -= eta;
    EvaluateSecularLast( rho, z, state, ctrl.penalizeDerivative );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularEVDInfo
SecularInner
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        State<Real>& state,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
      if( k == n-1 )
          LogicError("SecularInner meant for inner eigenvalues");
      if( n <= 2 )
          LogicError("SecularInner meant for n > 2");
    )

    SecularEVDInfo info;
    state.dMinusShift.Resize(n,1);
    state.diagDiff = d(k+1) - d(k);

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
          if( !limits::IsFinite(state.rootEst) )
          {
              RuntimeError("Produced non-finite rootEst=",state.rootEst);
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
SecularEVDInfo
SecularLast
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        LastState<Real>& state,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real eps = limits::Epsilon<Real>();
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
      if( k != n-1 )
          LogicError("SecularLast meant for largest eigenvalue");
      if( n <= 2 )
          LogicError("SecularLast meant for n > 2");
    )

    SecularEVDInfo info;
    state.dMinusShift.Resize(n,1);
    state.diagDiff = d(n-1) - d(n-2);

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

} // namespace secular_evd

// Compute an eigenvalue corresponding to the diagonal plus rank-one matrix
//
//     diag(d) + rho z z^T,
//
// where || z ||_2 = 1, with
//
//     d(0) < d(1) < ... < d(n-1)
//
// and rho > 0.
//
// This routine loosely corresponds to LAPACK's {s,d}laed4 [CITATION].
//

template<typename Real,typename>
SecularEVDInfo
SecularEigenvalue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& eigenvalue,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
      const Real zero(0);
      if( k < 0 || k >= n )
          LogicError("Invalid eigenvalue request");
      for( Int j=0; j<n-1; ++j )
          if( d(j) >= d(j+1) )
              LogicError("Assumption that d(j) < d(j+1) broken");
      if( rho <= zero )
          LogicError("Assumption that rho > 0 was broken");
      if( z.Height() != n )
          LogicError("z was not the correct height");
      // TODO(poulson): Check the assumption that || z ||_2 = 1
    )

    SecularEVDInfo info;
    if( n == 1 )
    {
        eigenvalue = d(0) + rho;
        return info;
    }
    else if( n == 2 )
    {
        eigenvalue = secular_evd::TwoByTwo( k, d(0), d(1), rho, z(0), z(1) );
        return info;
    }

    if( k < n-1 )
    {
        secular_evd::State<Real> state;
        info = secular_evd::SecularInner( k, d, rho, z, state, ctrl );
        eigenvalue = state.rootEst;
    }
    else
    {
        secular_evd::LastState<Real> state;
        info = secular_evd::SecularLast( k, d, rho, z, state, ctrl );
        eigenvalue = state.rootEst;
    }

    return info;
}

template<typename Real,typename>
SecularEVDInfo
SecularEigenvalue
( Int whichValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Real& eigenvalue,
        Matrix<Real>& dMinusShift,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int k = whichValue;
    const Int n = d.Height();
    EL_DEBUG_ONLY(
      const Real zero(0);
      if( k < 0 || k >= n )
          LogicError("Invalid eigenvalue request");
      for( Int j=0; j<n-1; ++j )
          if( d(j) >= d(j+1) )
              LogicError("Assumption that d(j) < d(j+1) broken");
      if( rho <= zero )
          LogicError("Assumption that rho > 0 was broken");
      if( z.Height() != n )
          LogicError("z was not the correct height");
      // TODO(poulson): Check the assumption that || z ||_2 = 1
    )

    SecularEVDInfo info;
    dMinusShift.Resize(n,1);
    if( n == 1 )
    {
        eigenvalue = d(0) + rho;
        // TODO(poulson): Make this computation more accurate for completeness?
        dMinusShift(0) = d(0) - eigenvalue;
        return info;
    }
    else if( n == 2 )
    {
        eigenvalue =
          secular_evd::TwoByTwo
          ( k, d(0), d(1), rho, z(0), z(1), dMinusShift(0), dMinusShift(1) );
        return info;
    }

    if( k < n-1 )
    {
        secular_evd::State<Real> state;
        info = secular_evd::SecularInner( k, d, rho, z, state, ctrl );
        eigenvalue = state.rootEst;
        dMinusShift = state.dMinusShift;
    }
    else
    {
        secular_evd::LastState<Real> state;
        info = secular_evd::SecularLast( k, d, rho, z, state, ctrl );
        eigenvalue = state.rootEst;
        dMinusShift = state.dMinusShift;
    }

    return info;
}

template<typename Real,typename>
SecularEVDInfo
SecularEVD
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const SecularEVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = d.Height();
    SecularEVDInfo info;

    w.Resize( n, 1 );
    if( n == 0 )
    {
        Q.Resize( n, n );
        return info;
    }

    // TODO(poulson): Batch secular equation solvers?
    //
    // Compute all of the eigenvalues and the vector r ~= sqrt(rho) z which
    // would produce the given eigenvalues to high relative accuracy.
    //
    // We accumulate the products involved in computing r in an online manner
    // and take the square-root of the absolute value (and the correct sign)
    // at the end. The key is to recognize that the only term left out of entry
    // i of the corrected vector in Eq. (3.6) of Gu/Eisenstat in the product
    //
    //    prod_{k=0}^{n-1} (lambda_k - d(i)) / (d(k) - d(i))
    //
    // is 1 / (d(i) - d(i)) (and we emphasize that the numerator is kept).
    // It is thus easy to see how the first term and second product in
    //
    //    r(i)^2 = (lambda_{n-1} - d(i)) *
    //      prod_{k=0}^{i-1} (lambda_k - d(i)) / (d(k  ) - d(i)) *
    //      prod_{k=i}^{n-2} (lambda_k - d(i)) / (d(k+1) - d(i))
    //
    // can be reorganized to yield
    //
    //    r(i)^2 = (lambda_{i} - d(i)) *
    //      prod_{k=0  }^{i-1} (lambda_k - d(i)) / (d(k) - d(i)) *
    //      prod_{k=i+1}^{n-1} (lambda_k - d(i)) / (d(k) - d(i)).
    //
    // Given that we have available only one lambda_j at a time, we greedily
    // update all r(i) products given each lambda_j
    // (Cf. LAPACK's {s,d}lasd8 [CITATION] for this approach).
    //
    Q.Resize( n, n );
    Matrix<Real> r;
    Ones( r, n, 1 );
    for( Int j=0; j<n; ++j )
    {
        // While we will temporarily store dMinusShift in the j'th column of Q.
        auto q = Q(ALL,IR(j));
        auto valueInfo = SecularEigenvalue( j, d, rho, z, w(j), q, ctrl );
        
        info.numIterations += valueInfo.numIterations;
        info.numAlternations += valueInfo.numAlternations;
        info.numCubicIterations += valueInfo.numCubicIterations;
        info.numCubicFailures += valueInfo.numCubicFailures;

        r(j) *= q(j);
        for( Int k=0; k<n; ++k )
        {
            if( k == j )
                continue;
            r(k) *= q(k) / (d(j)-d(k));
        }
    }
    for( Int j=0; j<n; ++j )
        r(j) = Sgn(z(j),false) * Sqrt(Abs(r(j)));

    for( Int j=0; j<n; ++j )
    {
        // Compute the j'th eigenvectors via Eqs. (3.4) and (3.3), respectively.
        auto q = Q(ALL,IR(j));
        for( Int i=0; i<n; ++i )
        {
            q(i) = r(i) / q(i);
        }
        q *= Real(1) / FrobeniusNorm( q );
    }

    return info;
}

#define PROTO(Real) \
  template SecularEVDInfo \
  SecularEigenvalue \
  ( Int whichValue, \
    const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& z, \
          Real& eigenvalue, \
    const SecularEVDCtrl<Real>& ctrl ); \
  template SecularEVDInfo \
  SecularEigenvalue \
  ( Int whichValue, \
    const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& z, \
          Real& eigenvalue, \
          Matrix<Real>& dMinusShift, \
    const SecularEVDCtrl<Real>& ctrl ); \
  template SecularEVDInfo \
  SecularEVD \
  ( const Matrix<Real>& d, \
    const Real& rho, \
    const Matrix<Real>& z, \
          Matrix<Real>& w, \
          Matrix<Real>& Q, \
    const SecularEVDCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT

#include <El/macros/Instantiate.h>

} // namespace El
