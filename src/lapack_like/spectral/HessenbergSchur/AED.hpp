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

// The primary references for the sequential implementation of Aggressive Early
// Deflation are
//
//   Karen Braman, Ralph Byers, and Roy Mathias,
//   "The multishift QR algorithm. Part II: Aggressive Early Deflation",
//   SIAM J. Matrix Anal. Appl., Vol. 23, No. 4, pp. 948--973, 2002.
//
// and the LAPACK implementations {S,D}LAQR2; the latter have several distinct
// differences from the suggestions of Braman et al., such as:
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
// The primary references for the parallel implementations are
//
//   Robert Granat, Bo Kagstrom, and Daniel Kressner,
//   "A novel parallel QR algorithm for hybrid distributed memory HPC systems",
//   SIAM J. Sci. Comput., Vol. 32, No. 4, pp. 2345--2378, 2010.
//
// and
//
//   Robert Granat, Bo Kagstrom, Daniel Kressner, and Meiyue Shao,
//   "Parallel library software for the multishift QR algorithm with Aggressive
//   Early Deflation", ACM Transactions on Mathematical Software, Vol. 41,
//   No. 4, 29 pages, 2015.
//

template<typename Field>
HessenbergSchurInfo
AED
( Matrix<Field>& H,
  Matrix<Complex<Base<Field>>>& w,
  Matrix<Field>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
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
    Matrix<Field> U, W, WAccum;
    auto ctrlSub( ctrl );

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);

    Int decreaseLevel = -1;
    Matrix<Field> hSubIter;
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
        // TODO(poulson): Describe why we don't use DetectSmallSubdiagonal
        Int iterBeg=winEnd-1;
        for( ; iterBeg>winBeg; --iterBeg )
            if( H(iterBeg,iterBeg-1) == Field(0) )
                break;
        if( ctrl.progress )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
        }
        // TODO(poulson): Describe why we don't switch to a simpler algorithm
        // here if the iteration window size is sufficiently small

        auto HIter = H( IR(iterBeg,winEnd), IR(iterBeg,winEnd) );
        GetDiagonal( HIter, hSubIter, -1 );
        aed::UpdateDeflationSize
        ( deflationSize, decreaseLevel, deflationSizeRec, numIterSinceDeflation,
          numStaleIterBeforeExceptional, hSubIter );

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

template<typename Field>
HessenbergSchurInfo
AED
( DistMatrix<Field,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Base<Field>>,STAR,STAR>& w,
  DistMatrix<Field,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();

    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;

    Int minMultiBulgeSize = Max( ctrl.minMultiBulgeSize, 4 );
    // TODO(poulson): Implement a more reasonable/configurable means of deciding
    // when to call the sequential implementation
    minMultiBulgeSize = Max( minMultiBulgeSize, 2*blockSize );
    // This maximum is meant to account for parallel overheads and needs to be
    // more principled (and perhaps based upon the number of workers and the
    // cluster characteristics)
    minMultiBulgeSize = Max( minMultiBulgeSize, ctrl.minDistMultiBulgeSize );

    HessenbergSchurInfo info;

    w.Resize( n, 1 );
    if( winSize < minMultiBulgeSize )
    {
        return multibulge::RedundantlyHandleWindow( H, w, Z, ctrl );
    }

    const Int numShiftsRec = ctrl.numShifts( n, winSize );
    const Int deflationSizeRec = ctrl.deflationSize( n, winSize, numShiftsRec );
    if( ctrl.progress && grid.Rank() == 0 )
    {
        Output
        ("Recommending ",numShiftsRec," shifts and a deflation window of size ",
         deflationSizeRec);
    }
    Int deflationSize = deflationSizeRec;

    // For multibulge::Sweep
    auto ctrlSub( ctrl );

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);

    Int decreaseLevel = -1;
    DistMatrix<Field,STAR,STAR> hMainWin(grid), hSubWin(grid);
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
        // TODO(poulson): Describe why we don't use DetectSmallSubdiagonal
        // TODO(poulson): Avoid the gather of the main diagonal in cases where
        // the QR sweep will be skipped
        util::GatherBidiagonal( H, IR(winBeg,winEnd), hMainWin, hSubWin );
        Int iterBeg=winEnd-1;
        for( ; iterBeg>winBeg; --iterBeg )
            if( hSubWin.GetLocal(iterBeg-1,0) == Field(0) )
                break;
        if( ctrl.progress && grid.Rank() == 0 )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
        }
        if( winEnd-iterBeg < minMultiBulgeSize )
        {
            // The window is small enough to switch to the simple scheme
            if( ctrl.progress && grid.Rank() == 0 )
                Output("Redundantly handling window [",iterBeg,",",winEnd,"]");
            auto ctrlIter( ctrl );
            ctrlIter.winBeg = iterBeg;
            ctrlIter.winEnd = winEnd;
            multibulge::RedundantlyHandleWindow( H, w, Z, ctrlIter );

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
            continue;
        }

        auto hSubIter = hSubWin( IR(iterBeg-winBeg,winEnd-1), ALL );
        aed::UpdateDeflationSize
        ( deflationSize, decreaseLevel, deflationSizeRec, numIterSinceDeflation,
          numStaleIterBeforeExceptional, hSubIter.Matrix() );

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
                H, hMainWin.Matrix(), hSubWin.Matrix(), w, ctrl );

            // Perform a small-bulge sweep
            auto wSub = w(IR(shiftBeg,winEnd),ALL);
            ctrlSub.winBeg = iterBeg;
            ctrlSub.winEnd = winEnd;
            multibulge::Sweep( H, wSub, Z, ctrlSub );
        }
        else if( ctrl.progress && grid.Rank() == 0 )
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
