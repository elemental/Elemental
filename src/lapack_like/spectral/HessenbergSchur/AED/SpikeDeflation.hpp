/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_AED_SPIKE_DEFLATION_HPP
#define EL_SCHUR_HESS_AED_SPIKE_DEFLATION_HPP

namespace El {
namespace hess_schur {

struct AEDInfo
{
    Int numUnconverged=0;
    Int numShiftCandidates=0;
    Int numDeflated=0;
};

namespace aed {

template<typename Real,typename=EnableIf<IsReal<Real>>>
AEDInfo SpikeDeflation
(       Matrix<Real>& T,
        Matrix<Real>& V,
  const Real& eta,
        Int numUnconverged,
        vector<Real>& work )
{
    DEBUG_CSE

    const Int n = T.Height();
    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);
     
    work.resize( Max(n,work.size()) );

    // Clear the two diagonals below the upper-Hessenberg portion for
    // SchurExchange
    for( Int i=0; i<n-3; ++i )
    {
        T(i+2,i) = zero;
        T(i+3,i) = zero;
    }
    if( n >= 3 )
        T(n-1,n-3) = zero;

    Int winBeg = numUnconverged;
    Int winEnd = n;
    while( winBeg < winEnd )
    {
        const bool twoByTwo =
          ( winEnd==1 ? false : T(winEnd-1,winEnd-2) != zero );

        if( twoByTwo )
        {
            // Follow LAPACK's suit (rather than Braman et al.) and use the 
            // eigenvalue of the 2x2 with largest magnitude in order to 
            // determine if this entry of the spike qualifies for 
            // "nearby-diagonal" deflation. Recall that the 2x2 block is assumed
            // to be in standard form,
            //
            //    | alpha, gamma |, where beta*gamma < 0,
            //    | beta,  alpha |
            //
            // so that the eigenvalues are alpha +- sqrt(beta*gamma) and the
            // spectral radius is |alpha| + sqrt(beta*gamma).
            //
            const Real& alpha = T(winEnd-2,winEnd-2);
            const Real& beta  = T(winEnd-1,winEnd-2); 
            const Real& gamma = T(winEnd-2,winEnd-1);
            const Real spectralRadius =
                Abs(alpha) + Sqrt(Abs(beta))*Sqrt(Abs(gamma));
            const Real scale =
              ( spectralRadius > 0 ? spectralRadius : Abs(eta) );

            // The relevant two entries of spike V^T [eta; zeros(n-1,1)]
            const Real sigma0 = eta*V(0,winEnd-2);
            const Real sigma1 = eta*V(0,winEnd-1);

            if( Max( Abs(sigma0), Abs(sigma1) ) <= Max( smallNum, ulp*scale ) )
            {
                // The two-by-two block satisfies the "nearby-diagonal" test
                winEnd -= 2;
            }
            else
            {
                // Move this undeflatable 2x2 block to the top of the window 
                lapack::SchurExchange
                ( n, &T(0,0), T.LDim(), &V(0,0), V.LDim(),
                  winEnd-2, winBeg, work.data() );
                winBeg += 2;
            }
        }
        else
        {
            const Real spectralRadius = Abs(T(winEnd-1,winEnd-1));
            const Real scale =
              ( spectralRadius > 0 ? spectralRadius : Abs(eta) );

            // The relevant entry of the spike V^T [eta; zeros(n-1,1)]
            if( Abs(eta)*Abs(V(0,winEnd-1)) <= Max( smallNum, ulp*scale ) )
            {
                // The one-by-one block satisfies the "nearby-diagonal" test
                winEnd -= 1;
            }
            else
            {
                // Move the undeflatable 1x1 block to the top of the window
                lapack::SchurExchange
                ( n, &T(0,0), T.LDim(), &V(0,0), V.LDim(),
                  winEnd-1, winBeg, work.data() );
                winBeg += 1;
            }
        }
    }

    AEDInfo info;
    info.numUnconverged = numUnconverged;
    info.numShiftCandidates = winBeg-numUnconverged;
    info.numDeflated = n-winBeg;
    return info;
}

template<typename Real>
AEDInfo SpikeDeflation
(       Matrix<Complex<Real>>& T,
        Matrix<Complex<Real>>& V,
  const Complex<Real>& eta,
        Int numUnconverged,
        vector<Complex<Real>>& work )
{
    DEBUG_CSE
    const Int n = T.Height();
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);
     
    work.resize( Max(n,work.size()) );

    Int winBeg = numUnconverged;
    Int winEnd = n;
    while( winBeg < winEnd )
    {
        const Real spectralRadius = OneAbs(T(winEnd-1,winEnd-1));
        const Real scale =
          ( spectralRadius > 0 ? spectralRadius : OneAbs(eta) );

        // The relevant entry of the spike V' [eta; zeros(n-1,1)]
        if( OneAbs(eta)*OneAbs(V(0,winEnd-1)) <= Max( smallNum, ulp*scale ) )
        {
            // The one-by-one block satisfies the "nearby-diagonal" test
            winEnd -= 1;
        }
        else
        {
            // Move the undeflatable 1x1 block to the top of the window
            lapack::SchurExchange
            ( n, &T(0,0), T.LDim(), &V(0,0), V.LDim(), winEnd-1, winBeg );
            winBeg += 1;
        }
    }

    AEDInfo info;
    info.numUnconverged = numUnconverged;
    info.numShiftCandidates = winBeg-numUnconverged;
    info.numDeflated = n-winBeg;
    return info;
}

} // namespace aed

} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_AED_SPIKE_DEFLATION_HPP
