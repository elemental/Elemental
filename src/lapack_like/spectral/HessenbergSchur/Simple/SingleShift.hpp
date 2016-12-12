/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_SINGLE_SHIFT_HPP
#define EL_HESS_SCHUR_SINGLE_SHIFT_HPP

#include "../Util/MakeSubdiagonalReal.hpp"
#include "./SingleShift/Sweep.hpp"

namespace El {
namespace hess_schur {

template<typename Real>
HessenbergSchurInfo
SingleShift
(       Matrix<Complex<Real>>& H,
        Matrix<Complex<Real>>& w,
        Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> F;
    const Real zero(0), threeFourths=Real(3)/Real(4);
    const Int maxIter=30;    
    const Int n = H.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int windowSize = winEnd - winBeg;
    HessenbergSchurInfo info;

    w.Resize( n, 1 );

    if( windowSize == 0 )
    {
        return info;
    }
    if( windowSize == 1 )
    {
        w(winBeg) = H(winBeg,winBeg);
        return info;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=winBeg; j<winEnd-3; ++j ) 
    {
        H(j+2,j) = zero;
        H(j+3,j) = zero;
    }
    if( winBeg <= winEnd-3 )
        H(winEnd-1,winEnd-3) = zero;
    
    // The optimized sweeping procedure assumes/maintains a real subdiagonal
    util::MakeSubdiagonalReal( H, Z, ctrl );

    // Attempt to converge the eigenvalues one at a time
    auto ctrlSweep( ctrl );
    while( winBeg < winEnd )
    {
        Int iterBeg = winBeg;
        Int iter;
        for( iter=0; iter<maxIter; ++iter )
        {
            auto winInd = IR(iterBeg,winEnd);
            iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
            if( iterBeg > winBeg )
            {
                H(iterBeg,iterBeg-1) = zero;
            }
            if( iterBeg == winEnd-1 )
            {
                w(iterBeg) = H(iterBeg,iterBeg);
                --winEnd;
                break;
            }

            // Pick either the Wilkinson shift or an exceptional shift
            F shift;
            if( iter == maxIter/3 )
            {
                const F diagVal = H(iterBeg,iterBeg);
                const Real subdiagVal = RealPart(H(iterBeg+1,iterBeg));
                shift = diagVal + threeFourths*Abs(subdiagVal);
            } 
            else if( iter == 2*maxIter/3 )
            {
                const F diagVal = H(winEnd-1,winEnd-1);
                const Real subdiagVal = RealPart(H(winEnd-1,winEnd-2));
                shift = diagVal + threeFourths*Abs(subdiagVal);
            }
            else
            {
                auto subInd = IR(iterBeg,winEnd);
                shift = WilkinsonShift( H(subInd,subInd) );
            }

            ctrlSweep.winBeg = iterBeg;
            ctrlSweep.winEnd = winEnd;
            single_shift::SweepOpt( H, shift, Z, ctrlSweep );
            ++info.numIterations;
        }
        if( iter == maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("QR iteration did not converge");
            else
                break;
        }
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_SINGLE_SHIFT_HPP
