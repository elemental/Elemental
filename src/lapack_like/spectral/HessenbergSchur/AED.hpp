/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_AED_HPP
#define EL_HESS_SCHUR_AED_HPP

#include "./Simple.hpp"
#include "./MultiBulge/Sweep.hpp"
#include "./AED/UpdateDeflationSize.hpp"
#include "./AED/ModifyShifts.hpp"
#include "./AED/SpikeDeflation.hpp"
#include "./AED/Nibble.hpp"

namespace El {
namespace hess_schur {

// The primary references for the following Aggressive Early Deflation
// implementation are
//
//   Karen Braman, Ralph Byers, and Roy Mathias,
//   "The multishift QR algorithm. Part II: Aggressive Early Deflation",
//   SIAM J. Matrix Anal. Appl., Vol. 23, No. 4, pp. 948--973, 2002
//
// and the LAPACK implementation DLAQR2, which has several distinct differences
// from the suggestions of Braman et al., such as:
//
//   1) Solely using "nearby-diagonal deflation" instead of Braman et al.'s 
//      suggestion of also allowing for "window-Schur deflation".
//
//   2) Using the largest (in magnitude) eigenvalue of a 2x2 Schur block to 
//      determine whether it qualifies for "nearby-diagonal deflation" rather
//      that using the square-root of the absolute value of the determinant
//      (which would correspond to the geometric mean of the eigenvalue
//       magnitudes). 
//
// In both respects, the LAPACK implementation is significantly more
// conservative than the original suggestions of Braman et al.
//

template<typename F>
HessenbergSchurInfo
AED
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE 
    typedef Base<F> Real;
    const Real zero(0);

    const Int n = H.Height();
    const Int minMultiBulgeSize = Max( ctrl.minMultiBulgeSize, 4 );
    HessenbergSchurInfo info;

    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    if( winSize < minMultiBulgeSize )
    {
        return Simple( H, w, Z, ctrl );
    }

    w.Resize( n, 1 );

    const Int numShiftsRec = ctrl.numShifts( n, winSize );
    const Int deflationSizeRec = ctrl.deflationSize( n, winSize, numShiftsRec );
    if( ctrl.progress )
    {
        Output
        ("Recommending ",numShiftsRec," shifts and a deflation window of size ",
         deflationSizeRec);
    }
    Int deflationSize = deflationSizeRec;

    // For multibulge::Sweep
    Matrix<F> U, W, WAccum;
    auto ctrlSub( ctrl );

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);

    Int decreaseLevel = -1;
    while( winBeg < winEnd )
    {
        if( info.numIterations >= maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("AED QR iteration did not converge");
            else
                break;
        }

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        Int iterBeg=winEnd-1;
        for( ; iterBeg>winBeg; --iterBeg )
            if( H(iterBeg,iterBeg-1) == zero ) 
                break;
        if( ctrl.progress )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
        }
        const Int iterWinSize = winEnd-iterBeg;
        aed::UpdateDeflationSize
        ( deflationSize, decreaseLevel, deflationSizeRec, numIterSinceDeflation,
          numStaleIterBeforeExceptional, iterWinSize, winEnd, H );

        // Run AED on the bottom-right window of size deflationSize
        ctrlSub.winBeg = iterBeg;
        ctrlSub.winEnd = winEnd;
        auto deflateInfo = aed::Nibble( H, deflationSize, w, Z, ctrlSub );
        const Int numDeflated = deflateInfo.numDeflated;
        winEnd -= numDeflated;
        Int shiftBeg = winEnd - deflateInfo.numShiftCandidates;

        const Int newIterWinSize = winEnd-iterBeg;
        const Int sufficientDeflation = ctrl.sufficientDeflation(deflationSize);
        // TODO(poulson): Provide an explanation for this strategy for when to
        // avoid sweeps
        if( numDeflated == 0 ||
          (numDeflated <= sufficientDeflation && 
           newIterWinSize >= minMultiBulgeSize) )
        {
            shiftBeg =
              aed::ModifyShifts
              ( numShiftsRec, newIterWinSize, numIterSinceDeflation, 
                numStaleIterBeforeExceptional, winBeg, winEnd, shiftBeg,
                H, w, ctrl );

            // Perform a small-bulge sweep
            auto wSub = w(IR(shiftBeg,winEnd),ALL); 
            ctrlSub.winBeg = iterBeg;
            ctrlSub.winEnd = winEnd;
            multibulge::Sweep( H, wSub, Z, U, W, WAccum, ctrlSub );
        }
        else if( ctrl.progress )
            Output("  Skipping QR sweep");

        ++info.numIterations;
        if( numDeflated > 0 )
            numIterSinceDeflation = 0;
        else
            ++numIterSinceDeflation;
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_AED_HPP
