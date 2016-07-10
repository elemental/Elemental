/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVD_BIDIAG_QR_HPP
#define EL_SVD_BIDIAG_QR_HPP

namespace El {
namespace svd {

// The estimate of sigma_min(B) follows 
//
//   N. J. Higham,
//   "Efficient algorithms for computing the condition number of a
//    tridiagonal matrix",
//   SIAM J. Sci. Stat. Comput. 7 (1986), pp. 150-65 [CITATION],
//
// but with the modification of working with the inverses of the entries of
// the solution to the triangular system of equations.

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real InverseInfinityNormOfBidiagInverse
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag )
{
    DEBUG_CSE
    const Int n = mainDiag.Height();
    const Real zero(0);
    DEBUG_ONLY(
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

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real InverseOneNormOfBidiagInverse
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag )
{
    DEBUG_CSE
    const Int n = mainDiag.Height();
    const Real zero(0);
    DEBUG_ONLY(
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
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MinSingularValueEstimateOfBidiag
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag,
  bool looseBound )
{
    DEBUG_CSE
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

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MaxSingularValueEstimateOfBidiag
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag )
{
    DEBUG_CSE
    return Max( MaxNorm(mainDiag), MaxNorm(superDiag) );
}

namespace bidiag_qr {

template<typename F>
void QRSweep
(       Matrix<Base<F>>& mainDiag,
        Matrix<Base<F>>& superDiag,
        Matrix<F>& U,
        Matrix<F>& V,
  const Base<F>& shift,
        ForwardOrBackward direction,  
        Matrix<Base<F>>& cUList,
        Matrix<Base<F>>& sUList,
        Matrix<Base<F>>& cVList,
        Matrix<Base<F>>& sVList,
  const BidiagQRCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
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
                mainDiag(i+1) = cU*mainDiag(i+1) - sU*superDiag(i+1);
                if( i < n-2 )
                {
                    g = sV*superDiag(i+1);
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
                    cUList(i) = cU;
                    sUList(i) = -sU;
                }
                if( ctrl.wantV )
                {
                    cVList(i) = cV;
                    sVList(i) = -sV;
                }
            }
            eta = mainDiag(0)*cU;
            mainDiag(0) = eta*cV;
            superDiag(0) = eta*sV;
        }
        else
        {
            // Run the classical sweep 
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
                g = sU*superDiag(i-1);
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
                    cUList(i) = cU;
                    sUList(i) = -sU;
                }
                if( ctrl.wantV )
                {
                    cVList(i) = cV;
                    sVList(i) = -sV;
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

template<typename F>
BidiagQRInfo
Helper
( Matrix<Base<F>>& mainDiag,
  Matrix<Base<F>>& superDiag,
  Matrix<Base<F>>& s,
  Matrix<F>& U,
  Matrix<F>& V,
  const BidiagQRCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = mainDiag.Height();
    const Int mU = U.Height();
    const Int mV = V.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real zero(0), negOne(-1), two(2);
    BidiagQRInfo info;

    Real tol = ctrl.tol;
    if( tol == zero )
    {
        // Cf. LAPACK's {s,d}bdsqr for this default choice
        tol = eps*Max( Real(10), Min(Real(100),Pow(eps,Real(-1)/Real(8))) );
    }

    const Real maxSingValEst =
      MaxSingularValueEstimateOfBidiag(mainDiag,superDiag);

    const Int maxInnerLoops = ctrl.maxIterPerVal*n*n;
    Real threshold = safeMin*maxInnerLoops;
    if( ctrl.relativeTol )
    {
        const Real minSingValEst =
          MinSingularValueEstimateOfBidiag
          ( mainDiag, superDiag, ctrl.looseMinSingValEst );
        // Cf. Demmel and Kahan for this choice
        threshold = Max( threshold, tol*minSingValEst );
    }
    else
    {
        threshold = Max( threshold, tol*maxSingValEst );
    }

    // 'winEnd' will always point just past the last unconverged singular value
    // and each iteration will select the largest possible unreduced bidiagonal
    //  window of the form [winBeg,winEnd).
    Int innerLoops=0, numZeroForwardInnerLoops=0,
                      numZeroBackwardInnerLoops=0,
                      numNonzeroForwardInnerLoops=0,
                      numNonzeroBackwardInnerLoops=0;
    Int winEnd = n;
    Int oldWinBeg=-1, oldWinEnd=-1;
    ForwardOrBackward direction = FORWARD;
    Matrix<Real> cUList(n,1), sUList(n,1), cVList(n,1), sVList(n,1);
    Matrix<Real> mainDiagSub, superDiagSub;
    Matrix<F> USub, VSub;
    while( winEnd > 0 )
    {
        if( innerLoops > maxInnerLoops )
        {
            if( ctrl.demandConverged )
                LogicError("Did not converge all singular values");
            info.numUnconverged = winEnd;
        }

        if( !ctrl.relativeTol && Abs(mainDiag(winEnd-1)) <= threshold )
        {
            mainDiag(winEnd-1) = zero;
        }

        Real winMaxSingValEst = Abs(mainDiag(winEnd-1));
        Int winBeg = 0;
        for( Int j=winEnd-2; j>=0; --j )
        {
            const Real alphaAbs = Abs(mainDiag(j));
            if( !ctrl.relativeTol && alphaAbs <= threshold )
            {
                mainDiag(j) = zero;
            }
            const Real betaAbs = Abs(superDiag(j));
            if( betaAbs <= threshold )
            {
                superDiag(j) = zero;
                winBeg = j+1;
                break;
            } 
            winMaxSingValEst = Max( winMaxSingValEst, alphaAbs );
            winMaxSingValEst = Max( winMaxSingValEst, betaAbs );
        }

        if( winBeg+1 == winEnd )
        {
            // The bottom singular value converged
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
                TwoByTwoUpper
                ( mainDiag(winBeg), superDiag(winBeg), mainDiag(winBeg+1), 
                  sigmaMax, sgnMax, sigmaMin, sgnMin, cU, sU, cV, sV ); 
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
                TwoByTwoUpper
                ( mainDiag(winBeg), superDiag(winBeg), mainDiag(winBeg+1), 
                  sigmaMax, sigmaMin ); 
            }

            mainDiag(winBeg) = sigmaMax;
            superDiag(winBeg) = zero;
            mainDiag(winBeg+1) = sigmaMin;
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
            }
            else
            {
                // We are roughly graded upward, so chase a bulge upward
                direction = BACKWARD;
            }
        }

        // Check for convergence
        // TODO(poulson): Document the Demmel/Kahan tests employed here
        Real winMinSingValEst = zero;
        if( direction == FORWARD )
        {
            if( Abs(superDiag(winEnd-2)) <= tol*Abs(mainDiag(winEnd-1)) ||
                (!ctrl.relativeTol && Abs(superDiag(winEnd-2)) <= threshold) )
            {
                superDiag(winEnd-2) = zero;
                continue;
            }
            if( ctrl.relativeTol )
            {
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
        }
        else
        {
            if( Abs(superDiag(winBeg)) <= tol*Abs(mainDiag(winBeg)) ||
                (!ctrl.relativeTol && Abs(superDiag(winBeg)) <= threshold) )
            {
                superDiag(winBeg) = zero;
                continue;
            }
            if( ctrl.relativeTol )
            {
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
            }
        }

        // No deflation checks succeeded, so save the window before shifting
        oldWinBeg = winBeg;
        oldWinEnd = winEnd;
        innerLoops += winEnd-winBeg;

        Real shift = zero;
        if( ctrl.relativeTol &&
            n*tol*(winMinSingValEst/winMaxSingValEst) <= Max(eps,tol/100) )
        {
            shift = zero; 
        }
        else
        {
            Real alphaAbs, sigmaMax, sigmaMin;
            if( direction == FORWARD )
            {
                alphaAbs = Abs(mainDiag(winBeg));
                TwoByTwoUpper
                ( mainDiag(winEnd-2), superDiag(winEnd-2), mainDiag(winEnd-1),
                  sigmaMax, sigmaMin );
            }
            else
            {
                alphaAbs = Abs(mainDiag(winEnd-1));
                TwoByTwoUpper
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
                numZeroForwardInnerLoops += winEnd-winBeg;
            else
                numZeroBackwardInnerLoops += winEnd-winBeg;
        }
        else
        {
            if( direction == FORWARD )
                numNonzeroForwardInnerLoops += winEnd-winBeg;
            else
                numNonzeroBackwardInnerLoops += winEnd-winBeg;
        }

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
        QRSweep
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
    for( Int j=0; j<n; ++j )
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

    // TODO(poulson): Move this up a level?
    if( ctrl.wantU || ctrl.wantV )
    {
        auto sortPairs = TaggedSort( s, ASCENDING );
        for( Int j=0; j<n; ++j )
            s(j) = sortPairs[j].value;
        if( ctrl.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrl.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V );
    }
    else
    {
        Sort( s, ASCENDING );
    }

    // TODO: Extend info to include innerLoops
    return info;
}

} // namespace bidiag_qr

template<typename Real,typename>
BidiagQRInfo
BidiagQR
( Matrix<Real>& mainDiag,
  Matrix<Real>& superDiag,
  Matrix<Real>& s,
  const BidiagQRCtrl<Real>& ctrl )
{
    DEBUG_CSE
    auto ctrlMod( ctrl );
    ctrlMod.wantU = false;
    ctrlMod.wantV = false;
    Matrix<Real> U, V;
    return bidiag_qr::Helper( mainDiag, superDiag, s, U, V, ctrlMod );
}

template<typename F>
BidiagQRInfo
BidiagQR
( Matrix<Base<F>>& mainDiag,
  Matrix<Base<F>>& superDiag,
  Matrix<Base<F>>& s,
  Matrix<F>& U,
  Matrix<F>& V,
  const BidiagQRCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    const Int n = mainDiag.Height();
    auto ctrlMod( ctrl );

    ctrlMod.wantU = true;
    if( ctrl.accumulateU )
    {
        if( U.Width() != n )
            LogicError("U was an invalid width");
    }
    else
    {
        Identity( U, n, n );
    }

    ctrlMod.wantV = true;
    if( ctrl.accumulateV )
    {
        if( V.Width() != n )
            LogicError("V was an invalid width");
    }
    else
    {
        Identity( V, n, n );
    }

    return bidiag_qr::Helper( mainDiag, superDiag, s, U, V, ctrlMod );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_BIDIAG_QR_HPP
