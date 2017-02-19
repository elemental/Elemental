/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_SVD_QR_HPP
#define EL_BIDIAG_SVD_QR_HPP

namespace El {
namespace bidiag_svd {

// The estimate of sigma_min(B) follows
//
//   N. J. Higham,
//   "Efficient algorithms for computing the condition number of a
//    tridiagonal matrix",
//   SIAM J. Sci. Stat. Comput. 7 (1986), pp. 150-65 [CITATION],
//
// but with the modification of working with the inverses of the entries of
// the solution to the triangular system of equations.

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real InverseInfinityNormOfBidiagInverse
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag )
{
    EL_DEBUG_CSE
    const Int n = mainDiag.Height();
    const Real zero(0);
    EL_DEBUG_ONLY(
      if( n == 0 )
          LogicError("Requested inverse of norm of empty matrix");
    )
    if( n == 0 )
        return zero;

    Real lambda = Abs(mainDiag(n-1));
    if( lambda == zero )
        return zero;
    Real inverseInfNorm = lambda;
    for( Int j=n-2; j>=0; --j )
    {
        lambda = Abs(mainDiag(j))*(lambda/(lambda+Abs(superDiag(j))));
        if( lambda == zero )
            return zero;
        inverseInfNorm = Min( inverseInfNorm, lambda );
    }
    return inverseInfNorm;
}

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real InverseOneNormOfBidiagInverse
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag )
{
    EL_DEBUG_CSE
    const Int n = mainDiag.Height();
    const Real zero(0);
    EL_DEBUG_ONLY(
      if( n == 0 )
          LogicError("Requested inverse of norm of empty matrix");
    )
    if( n == 0 )
        return zero;

    Real mu = Abs(mainDiag(0));
    if( mu == zero )
        return zero;
    Real inverseOneNorm = mu;
    for( Int j=1; j<n; ++j )
    {
        mu = Abs(mainDiag(j))*(mu/(mu+Abs(superDiag(j-1))));
        if( mu == zero )
            return zero;
        inverseOneNorm = Min( inverseOneNorm, mu );
    }
    return inverseOneNorm;
}

// LAPACK's {s,d}bdsqr only make use of the || inv(B) ||_1 and
// produce a lower bound that is tight by a factor of n rather than
// sqrt(n). Given that Demmel and Kahan recommend also using || inv(B) ||_oo,
// it would be worth investigating why the factor of sqrt(n) was forfeited.
// In the mean time, we will enable the looser LAPACK approach via the
// boolean 'looseBound'.
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real MinSingularValueEstimateOfBidiag
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag,
  bool looseBound )
{
    EL_DEBUG_CSE
    const Real invOneNormOfInv =
      InverseOneNormOfBidiagInverse(mainDiag,superDiag);
    if( looseBound )
    {
        const Int n = mainDiag.Height();
        // This lower bound is tight by a factor of n
        return invOneNormOfInv / Sqrt(Real(n));
    }
    else
    {
        const Real invInfNormOfInv =
          InverseInfinityNormOfBidiagInverse(mainDiag,superDiag);
        // This lower bound is tight by a factor of sqrt(n)
        return Min( invOneNormOfInv, invInfNormOfInv );
    }
}

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real MaxSingularValueEstimateOfBidiag
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag )
{
    EL_DEBUG_CSE
    return Max( MaxNorm(mainDiag), MaxNorm(superDiag) );
}

namespace qr {

// Cf. LAPACK's {s,d}bdsqr for these sweep strategies.
template<typename Field>
void Sweep
(       Matrix<Base<Field>>& mainDiag,
        Matrix<Base<Field>>& superDiag,
        Matrix<Field>& U,
        Matrix<Field>& V,
  const Base<Field>& shift,
        ForwardOrBackward direction,
        Matrix<Base<Field>>& cUList,
        Matrix<Base<Field>>& sUList,
        Matrix<Base<Field>>& cVList,
        Matrix<Base<Field>>& sVList,
  const BidiagSVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = mainDiag.Height();
    const Real zero(0), one(1);

    // NOTE: These should never actually map down to malloc since the calling
    // routine should have initially set their size to an upper bound of the
    // window size
    if( ctrl.wantU )
    {
        cUList.Resize( n-1, 1 );
        sUList.Resize( n-1, 1 );
    }
    if( ctrl.wantV )
    {
        cVList.Resize( n-1, 1 );
        sVList.Resize( n-1, 1 );
    }

    // TODO(poulson): Optimized versions of the following (especially to avoid
    // temporaries to improve the performance for heap scalars)
    Real cU, sU, cV, sV, rho, eta;
    if( direction == FORWARD )
    {
        if( shift == zero )
        {
            // Run the zero-shift sweep
            cU = cV = one;
            for( Int i=0; i<n-1; ++i )
            {
                rho = Givens( mainDiag(i)*cV, superDiag(i), cV, sV );
                if( i > 0 )
                    superDiag(i-1) = sU*rho;

                mainDiag(i) = Givens( cU*rho, mainDiag(i+1)*sV, cU, sU );

                if( ctrl.wantU )
                {
                    cUList(i) = cU;
                    sUList(i) = sU;
                }
                if( ctrl.wantV )
                {
                    cVList(i) = cV;
                    sVList(i) = sV;
                }
            }
            eta = mainDiag(n-1)*cV;
            mainDiag(n-1) = eta*cU;
            superDiag(n-2) = eta*sU;
        }
        else
        {
            // Run the classical sweep

            // B' B has a diagonal of
            //
            //    alpha_0^2, beta_0^2 + alpha_1^2, beta_1^2 + alpha_2^2, ...
            //
            // and an off-diagonal of
            //
            //    alpha_0 beta_0, alpha_1 beta_1, ...
            //
            // so the Givens rotation to rotate the (0,1) entry of
            // B' B - shift^2 into the (0,0) entry is of the form
            //
            //    |  cV,  sV | | alpha_0^2 - shift^2 | = | phi |.
            //    | -sV,  cV | |    alpha_0 beta_0   |   |  0  |
            //
            // But the same rotations would be generated by
            //
            //    | alpha_0 - shift^2/alpha_0 |,
            //    |           beta_0          |
            //
            // and the top term can be safely computed as
            //
            //    (|alpha_0|-shift) (sgn(alpha_0)+shift/alpha_0).
            //
            // If alpha_0 was zero, we should have forced a zero-shift iteration
            // to push it down the diagonal.
            //
            Real f = (Abs(mainDiag(0))-shift)*
                     (Sgn(mainDiag(0),false)+shift/mainDiag(0));
            Real g = superDiag(0);

            for( Int i=0; i<n-1; ++i )
            {
                rho = Givens( f, g, cV, sV );
                if( i > 0 )
                    superDiag(i-1) = rho;

                f = cV*mainDiag(i) + sV*superDiag(i);
                superDiag(i) = cV*superDiag(i) - sV*mainDiag(i);
                g = sV*mainDiag(i+1);
                mainDiag(i+1) *= cV;
                mainDiag(i) = Givens( f, g, cU, sU );

                f = cU*superDiag(i) + sU*mainDiag(i+1);
                mainDiag(i+1) = cU*mainDiag(i+1) - sU*superDiag(i);
                if( i < n-2 )
                {
                    g = sU*superDiag(i+1);
                    superDiag(i+1) *= cU;
                }

                if( ctrl.wantU )
                {
                    cUList(i) = cU;
                    sUList(i) = sU;
                }
                if( ctrl.wantV )
                {
                    cVList(i) = cV;
                    sVList(i) = sV;
                }
            }
            superDiag(n-2) = f;
        }
        if( ctrl.wantU )
        {
            ApplyGivensSequence
            ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD, cUList, sUList, U );
        }
        if( ctrl.wantV )
        {
            ApplyGivensSequence
            ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD, cVList, sVList, V );
        }
    }
    else
    {
        if( shift == zero )
        {
            // Run the zero-shift sweep
            cU = cV = one;
            for( Int i=n-1; i>0; --i )
            {
                rho = Givens( mainDiag(i)*cU, superDiag(i-1), cU, sU );
                if( i < n-1 )
                    superDiag(i) = sV*rho;

                mainDiag(i) = Givens( cV*rho, mainDiag(i-1)*sU, cV, sV );

                if( ctrl.wantU )
                {
                    cUList(i-1) = cU;
                    sUList(i-1) = -sU;
                }
                if( ctrl.wantV )
                {
                    cVList(i-1) = cV;
                    sVList(i-1) = -sV;
                }
            }
            eta = mainDiag(0)*cU;
            mainDiag(0) = eta*cV;
            superDiag(0) = eta*sV;
        }
        else
        {
            // Run the classical sweep

            // B B' has a diagonal of
            //
            //   beta_0^2 + alpha_0^2, beta_1^2 + alpha_1^2, ..., alpha_{n-1}^2,
            //
            // and an off-diagonal of
            //
            //    alpha_1 beta_0, alpha_2 beta_1, ...
            //
            // so the Givens rotation to rotate the (end,end-1) entry of
            // B B' - shift^2 into the (end,end) entry is of the form
            //
            //    |  cU,  sU | | alpha_{n-1}^2 - shift^2 | = | phi |.
            //    | -sU,  cU | |  alpha_{n-1} beta_{n-2} |   |  0  |
            //
            // But the same rotations would be generated by
            //
            //    | alpha_{n-1} - shift^2/alpha_{n-1} |,
            //    |            beta_{n-2}             |
            //
            // and the top term can be safely computed as
            //
            //    (|alpha_{n-1}|-shift) (sgn(alpha_{n-1})+shift/alpha_{n-1}).
            //
            // If alpha_{n-1} was zero, we should have forced a zero-shift
            // iteration to push it down the diagonal.
            Real f = (Abs(mainDiag(n-1))-shift)*
                     (Sgn(mainDiag(n-1),false)+shift/mainDiag(n-1));
            Real g = superDiag(n-2);

            for( Int i=n-1; i>0; --i )
            {
                rho = Givens( f, g, cU, sU );
                if( i < n-1 )
                    superDiag(i) = rho;

                f = cU*mainDiag(i) + sU*superDiag(i-1);
                superDiag(i-1) = cU*superDiag(i-1) - sU*mainDiag(i);
                g = sU*mainDiag(i-1);
                mainDiag(i-1) *= cU;
                mainDiag(i) = Givens( f, g, cV, sV );

                f = cV*superDiag(i-1) + sV*mainDiag(i-1);
                mainDiag(i-1) = cV*mainDiag(i-1) - sV*superDiag(i-1);
                if( i > 1 )
                {
                    g = sV*superDiag(i-2);
                    superDiag(i-2) *= cV;
                }

                if( ctrl.wantU )
                {
                    cUList(i-1) = cU;
                    sUList(i-1) = -sU;
                }
                if( ctrl.wantV )
                {
                    cVList(i-1) = cV;
                    sVList(i-1) = -sV;
                }
            }
            superDiag(0) = f;
        }
        if( ctrl.wantU )
        {
            ApplyGivensSequence
            ( RIGHT, VARIABLE_GIVENS_SEQUENCE, BACKWARD, cUList, sUList, U );
        }
        if( ctrl.wantV )
        {
            ApplyGivensSequence
            ( RIGHT, VARIABLE_GIVENS_SEQUENCE, BACKWARD, cVList, sVList, V );
        }
    }
}

// This implementation follows the sketch provided in
//
//   James Demmel and W. Kahan,
//   "Accurate Singular Values of Bidiagonal Matrices",
//   LAPACK Working Note 03, 1990 [CITATION]
//
// and implemented within LAPACK's {s,d}bdsqr [CITATION]. But note that DQDS
// is purposefully not called, as Elemental will make this decision at a higher
// level once a DQDS implementation is available.
//

template<typename Field>
bidiag_svd::QRInfo
Helper
( Matrix<Base<Field>>& mainDiag,
  Matrix<Base<Field>>& superDiag,
  Matrix<Field>& U,
  Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = mainDiag.Height();
    const Int mU = U.Height();
    const Int mV = V.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real zeroShiftFudge = n;
    const Real zero(0), negOne(-1), two(2);
    const bool relativeToSelfTol =
      (ctrl.tolType == RELATIVE_TO_SELF_SING_VAL_TOL);
    bidiag_svd::QRInfo info;

    const Real maxSingValEst =
      MaxSingularValueEstimateOfBidiag(mainDiag,superDiag);
    if( ctrl.progress )
        Output("Estimated || B ||_2 ~= ",maxSingValEst);

    Real tol = ctrl.tol;
    if( tol == zero )
    {
        // Demmel and Kahan recommend 100*eps for double-precision,
        // but the following choice, which is used within LAPACK's
        // {s,d}bdsqr [CITATION], interpolates between 10*eps and 100*eps,
        // depending upon whether eps is larger than 10^-16.
        tol = eps*Max( Real(10), Min(Real(100),Pow(eps,Real(-1)/Real(8))) );

        if( ctrl.tolType == ABSOLUTE_SING_VAL_TOL )
        {
            // No value was set, so choose an absolute tolerance that would be
            // equal to the default relative-to-max tolerance
            tol *= maxSingValEst;
        }
    }

    const Int maxInnerLoops = ctrl.qrCtrl.maxIterPerVal*n*n;
    if( ctrl.progress )
        Output("Set maxInnerLoops to ",maxInnerLoops);
    Real threshold = safeMin*maxInnerLoops;
    if( relativeToSelfTol )
    {
        const Real minSingValEst =
          MinSingularValueEstimateOfBidiag
          ( mainDiag, superDiag, ctrl.qrCtrl.looseMinSingValEst );
        if( ctrl.progress )
            Output("Estimated sigmaMin(B) ~= ",minSingValEst);
        // Cf. Demmel and Kahan for this choice
        threshold = Max( threshold, tol*minSingValEst );
    }
    else if( ctrl.tolType == RELATIVE_TO_MAX_SING_VAL_TOL )
    {
        threshold = Max( threshold, tol*maxSingValEst );
    }
    else
    {
        threshold = Max( threshold, tol );
    }
    if( ctrl.progress )
        Output("Set threshold to ",threshold);

    // 'winEnd' will always point just past the last unconverged singular value
    // and each iteration will select the largest possible unreduced bidiagonal
    //  window of the form [winBeg,winEnd).
    Int winEnd = n;
    Int oldWinBeg=-1, oldWinEnd=-1;
    ForwardOrBackward direction = FORWARD;
    Matrix<Real> cUList(n,1), sUList(n,1), cVList(n,1), sVList(n,1);
    Matrix<Real> mainDiagSub, superDiagSub;
    Matrix<Field> USub, VSub;
    while( winEnd > 0 )
    {
        if( info.numInnerLoops > maxInnerLoops )
        {
            if( ctrl.qrCtrl.demandConverged )
                LogicError("Did not converge all singular values");
            info.numUnconverged = winEnd;
            break;
        }

        if( !relativeToSelfTol && Abs(mainDiag(winEnd-1)) <= threshold )
        {
            if( ctrl.progress )
                Output("  Zeroing mainDiag(",winEnd-1,")=",mainDiag(winEnd-1));
            mainDiag(winEnd-1) = zero;
        }

        Real winMaxSingValEst = Abs(mainDiag(winEnd-1));
        Int winBeg = 0;
        for( Int j=winEnd-2; j>=0; --j )
        {
            const Real alphaAbs = Abs(mainDiag(j));
            if( !relativeToSelfTol && alphaAbs <= threshold )
            {
                if( ctrl.progress )
                    Output("  Zeroing mainDiag(",j,")=",mainDiag(j));
                mainDiag(j) = zero;
            }
            const Real betaAbs = Abs(superDiag(j));
            if( betaAbs <= threshold )
            {
                if( ctrl.progress )
                    Output("  Zeroing superDiag(",j,")=",superDiag(j));
                superDiag(j) = zero;
                winBeg = j+1;
                break;
            }
            winMaxSingValEst = Max( winMaxSingValEst, alphaAbs );
            winMaxSingValEst = Max( winMaxSingValEst, betaAbs );
        }
        if( ctrl.progress )
        {
            Output
            ("  After ",info.numInnerLoops," inner loops: window is [",
             winBeg,",",winEnd,")");
        }

        if( winBeg+1 == winEnd )
        {
            // The bottom singular value converged
            if( ctrl.progress )
                Output("  mainDiag(",winBeg,")=",mainDiag(winBeg)," converged");
            winEnd -= 1;
            continue;
        }
        else if( winBeg+2 == winEnd )
        {
            // Deflate the disconnected bottom-right 2x2 block
            // NOTE: We will worry about enforcing the positivity of the
            //       singular values later
            Real sigmaMin, sigmaMax;
            if( ctrl.wantU || ctrl.wantV )
            {
                Real sgnMax, sgnMin, cU, sU, cV, sV;
                svd::TwoByTwoUpper
                ( mainDiag(winBeg), superDiag(winBeg), mainDiag(winBeg+1),
                  sigmaMax, sgnMax, sigmaMin, sgnMin, cU, sU, cV, sV );
                sigmaMax *= sgnMax; // The signs will be fixed at the end
                sigmaMin *= sgnMin; // The signs will be fixed at the end
                if( ctrl.wantU )
                {
                    blas::Rot
                    ( mU, U.Buffer(0,winBeg  ), 1,
                          U.Buffer(0,winBeg+1), 1, cU, sU );
                }
                if( ctrl.wantV )
                {
                    blas::Rot
                    ( mV, V.Buffer(0,winBeg  ), 1,
                          V.Buffer(0,winBeg+1), 1, cV, sV );
                }
            }
            else
            {
                svd::TwoByTwoUpper
                ( mainDiag(winBeg), superDiag(winBeg), mainDiag(winBeg+1),
                  sigmaMax, sigmaMin );
            }

            mainDiag(winBeg) = sigmaMax;
            superDiag(winBeg) = zero;
            mainDiag(winBeg+1) = sigmaMin;
            if( ctrl.progress )
            {
                Output
                ("  mainDiag(",winBeg,")=",mainDiag(winBeg)," and ",
                 "mainDiag(",winBeg+1,")=",mainDiag(winBeg+1)," converged");
            }
            winEnd -= 2;
            continue;
        }

        if( winBeg >= oldWinEnd || winEnd <= oldWinBeg )
        {
            // We have a disjoint window, so choose a direction based upon the
            // simplest possible grading test
            if( Abs(mainDiag(winBeg)) >= Abs(mainDiag(winEnd-1)) )
            {
                // We are roughly graded downward, so chase a bulge downward
                direction = FORWARD;
                if( ctrl.progress )
                    Output("  Will chase FORWARD");
            }
            else
            {
                // We are roughly graded upward, so chase a bulge upward
                direction = BACKWARD;
                if( ctrl.progress )
                    Output("  Will chase BACKWARD");
            }
        }

        // Check for convergence
        // TODO(poulson): Document the Demmel/Kahan tests employed here
        Real winMinSingValEst = zero;
        if( direction == FORWARD )
        {
            if( Abs(superDiag(winEnd-2)) <= tol*Abs(mainDiag(winEnd-1)) ||
                (!relativeToSelfTol && Abs(superDiag(winEnd-2)) <= threshold) )
            {
                superDiag(winEnd-2) = zero;
                continue;
            }
            // Run an ad-hoc variant of Higham's || inv(B) ||_1 estimator
            Real mu = Abs(mainDiag(winBeg));
            winMinSingValEst = mu;
            bool deflated = false;
            for( Int j=winBeg; j<winEnd-1; ++j )
            {
                if( Abs(superDiag(j)) <= tol*mu )
                {
                    superDiag(j) = zero;
                    deflated = true;
                    break;
                }
                mu = Abs(mainDiag(j+1))*(mu/(mu+Abs(superDiag(j))));
                winMinSingValEst = Min( winMinSingValEst, mu );
            }
            if( deflated )
                continue;
        }
        else
        {
            if( Abs(superDiag(winBeg)) <= tol*Abs(mainDiag(winBeg)) ||
                (!relativeToSelfTol && Abs(superDiag(winBeg)) <= threshold) )
            {
                superDiag(winBeg) = zero;
                continue;
            }
            // Run an ad-hoc variant of Higham's || inv(B) ||_oo estimator
            Real lambda = Abs(mainDiag(winEnd-1));
            winMinSingValEst = lambda;
            bool deflated = false;
            for( Int j=winEnd-2; j>=winBeg; --j )
            {
                if( Abs(superDiag(j)) <= tol*lambda )
                {
                    superDiag(j) = zero;
                    deflated = true;
                    break;
                }
                lambda =
                  Abs(mainDiag(j))*(lambda/(lambda+Abs(superDiag(j))));
                winMinSingValEst = Min( winMinSingValEst, lambda );
            }
            if( deflated )
                continue;
        }
        const Real winCondEst = winMaxSingValEst / winMinSingValEst;

        // No deflation checks succeeded, so save the window before shifting
        oldWinBeg = winBeg;
        oldWinEnd = winEnd;
        info.numInnerLoops += winEnd-winBeg;
        ++info.numIterations;

        Real shift = zero;
        // We cannot simply use winMinSingValEst to additionally guard the
        // non high-relative-accuracy case since it was not actually computed
        // in said instance. Instead, we simply check if the relevant starting
        // diagonal entry is zero to avoid a divide-by-zero when applying the
        // implicit Q theorem within the upcoming sweep.
        const bool zeroStartingValue =
          direction == FORWARD ?
          mainDiag(winBeg) == zero :
          mainDiag(winEnd-1) == zero;
        const bool poorlyConditioned =
          zeroShiftFudge*tol/winCondEst <= Max(eps,tol/100);
        const bool extremelyPoorlyConditioned =
          zeroShiftFudge*tol/winCondEst <= eps;
        if( zeroStartingValue ||
            (relativeToSelfTol && poorlyConditioned) ||
            (!relativeToSelfTol && extremelyPoorlyConditioned) )
        {
            shift = zero;
        }
        else
        {
            Real alphaAbs, sigmaMax, sigmaMin;
            if( direction == FORWARD )
            {
                alphaAbs = Abs(mainDiag(winBeg));
                svd::TwoByTwoUpper
                ( mainDiag(winEnd-2), superDiag(winEnd-2), mainDiag(winEnd-1),
                  sigmaMax, sigmaMin );
            }
            else
            {
                alphaAbs = Abs(mainDiag(winEnd-1));
                svd::TwoByTwoUpper
                ( mainDiag(winBeg), superDiag(winBeg), mainDiag(winBeg+1),
                  sigmaMax, sigmaMin );
            }
            shift = sigmaMin;

            // Test if the shift would produce the same results as a zero shift
            // (if so, set to zero to use the faster zero-shift algorithm)
            if( alphaAbs != zero && Pow(shift/alphaAbs,two) < eps )
            {
                shift = zero;
            }
        }
        if( shift == zero )
        {
            if( direction == FORWARD )
            {
                info.numZeroShiftForwardInnerLoops += winEnd-winBeg;
                ++info.numZeroShiftForwardIterations;
            }
            else
            {
                info.numZeroShiftBackwardInnerLoops += winEnd-winBeg;
                ++info.numZeroShiftBackwardIterations;
            }
        }
        else
        {
            if( direction == FORWARD )
            {
                info.numNonzeroShiftForwardInnerLoops += winEnd-winBeg;
                ++info.numNonzeroShiftForwardIterations;
            }
            else
            {
                info.numNonzeroShiftBackwardInnerLoops += winEnd-winBeg;
                ++info.numNonzeroShiftBackwardIterations;
            }
        }
        if( ctrl.progress )
            Output("  Set shift to ",shift);

        // TODO(poulson): Decide if it is worthwhile to avoid the cost of these
        // views
        View( mainDiagSub, mainDiag, IR(winBeg,winEnd), ALL );
        View( superDiagSub, superDiag, IR(winBeg,winEnd-1), ALL );
        if( ctrl.wantU )
        {
            View( USub, U, ALL, IR(winBeg,winEnd) );
        }
        if( ctrl.wantV )
        {
            View( VSub, V, ALL, IR(winBeg,winEnd) );
        }
        Sweep
        ( mainDiagSub, superDiagSub, USub, VSub, shift, direction,
          cUList, sUList, cVList, sVList, ctrl );

        // Test for convergence of the last off-diagonal of the sweep
        if( direction == FORWARD )
        {
            if( Abs(superDiag(winEnd-2)) <= threshold )
                superDiag(winEnd-2) = zero;
        }
        else
        {
            if( Abs(superDiag(winBeg)) <= threshold )
                superDiag(winBeg) = zero;
        }
    }

    // Force the singular values to be positive (absorbing signs into V)
    for( Int j=0; j<info.numUnconverged; ++j )
        mainDiag(j) = Real(-1);
    for( Int j=info.numUnconverged; j<n; ++j )
    {
        if( mainDiag(j) < zero )
        {
            mainDiag(j) = -mainDiag(j);
            if( ctrl.wantV )
            {
                blas::Scal( mV, negOne, V.Buffer(0,j), 1 );
            }
        }
    }

    return info;
}

template<typename Real,
         typename=EnableIf<IsBlasScalar<Real>>>
bidiag_svd::QRInfo
LAPACKHelper
( Matrix<Real>& mainDiag,
  Matrix<Real>& superDiag,
  const BidiagSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = mainDiag.Height();
    bidiag_svd::QRInfo info;

    // While directly calling DQDS rather than BidiagQR appears to break the
    // conventions of Elemental's algorithm heirarchy, LAPACK's bidiagonal QR
    // algorithm immediately calls DQDS if singular vectors are not desired

    // lapack::BidiagDQDS expects superDiag to be of length n
    if( superDiag.LDim() < n )
    {
        Matrix<Real> superDiagPlus(n,1);
        superDiagPlus = superDiag;
        lapack::BidiagDQDS( n, mainDiag.Buffer(), superDiagPlus.Buffer() );
    }
    else
    {
        lapack::BidiagDQDS( n, mainDiag.Buffer(), superDiag.Buffer() );
    }

    return info;
}

template<typename Real,
         typename=DisableIf<IsBlasScalar<Real>>,
         typename=void>
bidiag_svd::QRInfo
LAPACKHelper
( Matrix<Real>& mainDiag,
  Matrix<Real>& superDiag,
  const BidiagSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    bidiag_svd::QRInfo info;
    LogicError("LAPACK does not support this datatype");
    return info;
}

} // namespace qr

template<typename Real>
bidiag_svd::QRInfo
QRAlg
( Matrix<Real>& mainDiag,
  Matrix<Real>& superDiag,
  const BidiagSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( superDiag.Height() != mainDiag.Height()-1 )
        LogicError("Invalid superDiag length");
    if( IsBlasScalar<Real>::value && ctrl.qrCtrl.useLAPACK )
    {
        return qr::LAPACKHelper( mainDiag, superDiag, ctrl );
    }

    auto ctrlMod( ctrl );
    ctrlMod.wantU = false;
    ctrlMod.wantV = false;
    Matrix<Real> U, V;
    return qr::Helper( mainDiag, superDiag, U, V, ctrlMod );
}

namespace qr {

#ifdef EL_HAVE_FLA_BSVD
template<typename Field,
         typename=EnableIf<IsBlasScalar<Field>>>
bidiag_svd::QRInfo
FLAMEHelper
( Matrix<Base<Field>>& mainDiag,
  Matrix<Base<Field>>& superDiag,
  Matrix<Field>& U,
  Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = mainDiag.Height();
    bidiag_svd::QRInfo info;

    flame::BidiagSVD
    ( n, U.Height(), V.Height(),
      mainDiag.Buffer(), superDiag.Buffer(),
      U.Buffer(), U.LDim(),
      V.Buffer(), V.LDim() );

    return info;
}

template<typename Field,
         typename=DisableIf<IsBlasScalar<Field>>,
         typename=void>
bidiag_svd::QRInfo
FLAMEHelper
( Matrix<Base<Field>>& mainDiag,
  Matrix<Base<Field>>& superDiag,
  Matrix<Field>& U,
  Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    bidiag_svd::QRInfo info;
    LogicError("libFLAME does not support this datatype");
    return info;
}
#endif // ifdef EL_HAVE_FLA_BSVD

template<typename Field,
         typename=EnableIf<IsBlasScalar<Field>>>
bidiag_svd::QRInfo
LAPACKHelper
( Matrix<Base<Field>>& mainDiag,
  Matrix<Base<Field>>& superDiag,
  Matrix<Field>& U,
  Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = mainDiag.Height();
    bidiag_svd::QRInfo info;

    Matrix<Field> VAdj;
    Adjoint( V, VAdj );

    // lapack::BidiagSVDQRAlg expects superDiag to be of length n
    if( superDiag.LDim() < n )
    {
        Matrix<Real> superDiagPlus(n,1);
        superDiagPlus = superDiag;
        lapack::BidiagSVDQRAlg
        ( UPPER, n, VAdj.Width(), U.Height(),
          mainDiag.Buffer(), superDiagPlus.Buffer(),
          VAdj.Buffer(), VAdj.LDim(),
          U.Buffer(), U.LDim() );
    }
    else
    {
        lapack::BidiagSVDQRAlg
        ( UPPER, n, VAdj.Width(), U.Height(),
          mainDiag.Buffer(), superDiag.Buffer(),
          VAdj.Buffer(), VAdj.LDim(),
          U.Buffer(), U.LDim() );
    }

    Adjoint( VAdj, V );

    return info;
}

template<typename Field,
         typename=DisableIf<IsBlasScalar<Field>>,
         typename=void>
bidiag_svd::QRInfo
LAPACKHelper
( Matrix<Base<Field>>& mainDiag,
  Matrix<Base<Field>>& superDiag,
  Matrix<Field>& U,
  Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    bidiag_svd::QRInfo info;
    LogicError("LAPACK does not support this datatype");
    return info;
}

} // namespace qr

template<typename Field>
bidiag_svd::QRInfo
QRAlg
( Matrix<Base<Field>>& mainDiag,
  Matrix<Base<Field>>& superDiag,
  Matrix<Field>& U,
  Matrix<Field>& V,
  const BidiagSVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = mainDiag.Height();
    if( superDiag.Height() != n-1 )
        LogicError("Invalid superDiag length");

    if( !ctrl.wantU && !ctrl.wantV )
    {
        return QRAlg( mainDiag, superDiag, ctrl );
    }

    if( ctrl.wantU )
    {
        if( ctrl.accumulateU )
        {
            if( U.Width() != n )
                LogicError("U was an invalid width");
        }
        else
        {
            Identity( U, n, n );
        }
    }

    if( ctrl.wantV )
    {
        if( ctrl.accumulateV )
        {
            if( V.Width() != n )
                LogicError("V was an invalid width");
        }
        else
        {
            Identity( V, n, n );
        }
    }

#ifdef EL_HAVE_FLA_BSVD
    if( IsBlasScalar<Field>::value && ctrl.qrCtrl.useFLAME )
    {
        return qr::FLAMEHelper( mainDiag, superDiag, U, V, ctrl );
    }
#endif
    if( IsBlasScalar<Field>::value && ctrl.qrCtrl.useLAPACK )
    {
        return qr::LAPACKHelper( mainDiag, superDiag, U, V, ctrl );
    }

    return qr::Helper( mainDiag, superDiag, U, V, ctrl );
}

} // namespace bidiag_svd
} // namespace El

#endif // ifndef EL_BIDIAG_SVD_QR_HPP
