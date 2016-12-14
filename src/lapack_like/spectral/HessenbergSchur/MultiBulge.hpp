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
#include "./MultiBulge/TwoByTwo.hpp"
#include "./MultiBulge/ComputeShifts.hpp"
#include "./MultiBulge/RedundantlyHandleWindow.hpp"
#include "./MultiBulge/Sweep.hpp"

namespace El {

namespace hess_schur {

template<typename Field>
HessenbergSchurInfo
MultiBulge
( Matrix<Field>& H,
  Matrix<Complex<Base<Field>>>& w,
  Matrix<Field>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
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
    Matrix<Field> U, W, WAccum;

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
        auto winInd = IR(winBeg,winEnd);

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        const Int iterOffset = DetectSmallSubdiagonal( H(winInd,winInd) );
        const Int iterBeg = winBeg + iterOffset;
        const Int iterWinSize = winEnd-iterBeg;
        if( iterOffset > 0 )
            H(iterBeg,iterBeg-1) = Field(0);
        if( iterWinSize == 1 )
        {
            w(iterBeg) = H(iterBeg,iterBeg);

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterWinSize == 2 )
        {
            multibulge::TwoByTwo( H, w, Z, iterBeg, ctrl );

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterWinSize < minMultiBulgeSize )
        {
            // The window is small enough to switch to the simple scheme
            auto ctrlSub( ctrl );
            ctrlSub.winBeg = iterBeg;
            ctrlSub.winEnd = winEnd;
            Simple( H, w, Z, ctrlSub );

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
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

template<typename Field>
HessenbergSchurInfo
MultiBulge
( DistMatrix<Field,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Base<Field>>,STAR,STAR>& w,
  DistMatrix<Field,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    const Grid& grid = H.Grid();

    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const Int blockSize = H.BlockHeight();

    // TODO(poulson): Implement a more reasonable/configurable means of deciding
    // when to call the sequential implementation
    Int minMultiBulgeSize = Max( ctrl.minMultiBulgeSize, 2*blockSize );
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
    DistMatrix<Field,STAR,STAR> hMainWin(grid), hSubWin(grid), hSuperWin(grid);
    while( winBeg < winEnd )
    {
        if( info.numIterations >= maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("MultiBulge QR iteration did not converge");
            else
                break;
        }
        auto winInd = IR(winBeg,winEnd);
        if( ctrl.progress && grid.Rank() == 0 )
            Output("winBeg=",winBeg,", winEnd=",winEnd);

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        // TODO(poulson): Have the interblock chase from the previous sweep
        // collect the main and sub diagonal of H along the diagonal workers
        // and then broadcast across the "cross" communicator.
        util::GatherTridiagonal( H, winInd, hMainWin, hSubWin, hSuperWin );

        const Int iterOffset =
          DetectSmallSubdiagonal
          ( hMainWin.Matrix(), hSubWin.Matrix(), hSuperWin.Matrix() );
        const Int iterBeg = winBeg + iterOffset;
        const Int iterWinSize = winEnd-iterBeg;
        if( iterOffset > 0 )
        {
            if( ctrl.progress && grid.Rank() == 0 )
                Output("iterOffset was ",iterOffset);
            H.Set( iterBeg, iterBeg-1, Field(0) );
            hSubWin.Set( iterOffset-1, 0, Field(0) );
        }
        if( iterWinSize == 1 )
        {
            if( ctrl.progress && grid.Rank() == 0 )
                Output("One-by-one window at ",iterBeg);
            w.Set( iterBeg, 0, hMainWin.GetLocal(iterOffset,0) );

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterWinSize == 2 )
        {
            if( ctrl.progress && grid.Rank() == 0 )
                Output("Two-by-two window at ",iterBeg);
            const Field eta00 = hMainWin.GetLocal(iterOffset,0);
            const Field eta01 = hSuperWin.GetLocal(iterOffset,0);
            const Field eta10 = hSubWin.GetLocal(iterOffset,0);
            const Field eta11 = hMainWin.GetLocal(iterOffset+1,0);
            multibulge::TwoByTwo
            ( H, eta00, eta01, eta10, eta11, w, Z, iterBeg, ctrl );

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterWinSize < minMultiBulgeSize )
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

        const Int numShiftsRec = ctrl.numShifts( n, iterWinSize );
        if( ctrl.progress && grid.Rank() == 0 )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
            Output("  recommending ",numShiftsRec," shifts");
        }

        // NOTE(poulson): In the case where exceptional shifts are used, the
        // main and subdiagonals of H in the window are currently redundantly
        // gathered. It could be worthwhile to pass in hMainWin and hSubWin.
        const Int shiftBeg = multibulge::ComputeShifts
        ( H, w, iterBeg, winBeg, winEnd, numShiftsRec, numIterSinceDeflation,
          numStaleIterBeforeExceptional, ctrlShifts );
        auto shiftInd = IR(shiftBeg,winEnd);
        auto wShifts = w(shiftInd,ALL);

        // Perform a small-bulge sweep
        auto ctrlSweep( ctrl );
        ctrlSweep.winBeg = iterBeg;
        ctrlSweep.winEnd = winEnd;
        multibulge::Sweep( H, wShifts, Z, ctrlSweep );

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
