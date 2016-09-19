/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_HPP
#define EL_HESS_SCHUR_MULTIBULGE_HPP

#include "./Simple.hpp"
#include "./MultiBulge/HandleTwoByTwo.hpp"
#include "./MultiBulge/ComputeShifts.hpp"
#include "./MultiBulge/Sweep.hpp"

// This is not yet functional but is included for testing reasons
#include "./MultiBulge/IntraBlockChase.hpp"

namespace El {
namespace hess_schur {

template<typename F>
HessenbergSchurInfo
MultiBulge
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE 
    typedef Base<F> Real;
    const Real zero(0);

    const Int n = H.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const Int minMultiBulgeSize = Max( ctrl.minMultiBulgeSize, 4 );
    HessenbergSchurInfo info;

    if( winSize < minMultiBulgeSize )
    {
        return Simple( H, w, Z, ctrl );
    }

    w.Resize( n, 1 );
    Matrix<F> U, W, WAccum;

    auto ctrlShifts( ctrl );
    ctrlShifts.winBeg = 0;
    ctrlShifts.winEnd = END;
    ctrlShifts.fullTriangle = false;

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);

    Int iterBegLast=-1, winEndLast=-1;
    while( winBeg < winEnd )
    {
        if( info.numIterations >= maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("MultiBulge QR iteration did not converge");
            else
                break;
        }

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        Int iterBeg = winBeg;
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
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterBeg == winEnd-2 )
        {
            multibulge::HandleTwoByTwo( H, w, Z, iterBeg, winEnd, ctrl );
            winEnd -= 2;
            numIterSinceDeflation = 0;
            continue;
        }

        const Int iterWinSize = winEnd-iterBeg;
        if( iterWinSize < minMultiBulgeSize )
        {
            // The window is small enough to switch to the simple scheme
            auto ctrlSub( ctrl );
            ctrlSub.winBeg = iterBeg;
            ctrlSub.winEnd = winEnd;
            Simple( H, w, Z, ctrlSub );
            winEnd = iterBeg;
            continue;
        }

        const Int numShiftsRec = ctrl.numShifts( n, iterWinSize );
        if( ctrl.progress )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
            Output("  recommending ",numShiftsRec," shifts");
        }

        const Int shiftBeg = multibulge::ComputeShifts
        ( H, w, iterBeg, winBeg, winEnd, numShiftsRec, numIterSinceDeflation,
          numStaleIterBeforeExceptional, ctrlShifts );
        auto shiftInd = IR(shiftBeg,winEnd);
        auto wShifts = w(shiftInd,ALL);

        // Perform a small-bulge sweep
        auto ctrlSweep( ctrl );
        ctrlSweep.winBeg = iterBeg;
        ctrlSweep.winEnd = winEnd;
        multibulge::Sweep( H, wShifts, Z, U, W, WAccum, ctrlSweep );

        ++info.numIterations;
        if( iterBeg == iterBegLast && winEnd == winEndLast )
            ++numIterSinceDeflation;
        iterBegLast = iterBeg;
        winEndLast = winEnd;
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_HPP
