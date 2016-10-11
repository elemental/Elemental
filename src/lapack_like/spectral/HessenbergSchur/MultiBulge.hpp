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
#include "./MultiBulge/Sweep.hpp"

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
        auto winInd = IR(winBeg,winEnd);

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        const Int iterOffset = DetectSmallSubdiagonal( H(winInd,winInd) );
        const Int iterBeg = winBeg + iterOffset;
        const Int iterWinSize = winEnd-iterBeg;
        if( iterOffset > 0 )
            H(iterBeg,iterBeg-1) = zero;
        if( iterWinSize == 1 )
        {
            w(iterBeg) = H(iterBeg,iterBeg);
            --winEnd;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterWinSize == 2 )
        {
            multibulge::TwoByTwo( H, w, Z, iterBeg, ctrl );
            winEnd -= 2;
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

namespace multibulge {

// TODO(poulson): Move this into MultiBulge/RedundantlyHandleWindow.hpp
template<typename F>
HessenbergSchurInfo
RedundantlyHandleWindow
( DistMatrix<F,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Base<F>>,STAR,STAR>& w,
  DistMatrix<F,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const Int blockSize = H.BlockHeight();

    auto winInd = IR(winBeg,winEnd);
    auto HWin = H(winInd,winInd);
    auto& wLoc = w.Matrix();

    // Compute the Schur decomposition HWin = ZWin TWin ZWin',
    // where HWin is overwritten by TWin, and wWin by diag(TWin).
    DistMatrix<F,STAR,STAR> HWinFull( HWin );
    auto wWin = wLoc(winInd,ALL);
    Matrix<F> ZWin;
    Identity( ZWin, winSize, winSize );
    auto info = HessenbergSchur( HWinFull.Matrix(), wWin, ZWin );
    HWin = HWinFull;

    if( ctrl.fullTriangle )
    {
        if( n > winEnd )
        {
            // Overwrite H(winInd,winEnd:n) *= ZWin'
            // (applied from the left)
            auto HRight = H( winInd, IR(winEnd,n) );
            // Since we only need to overwrite HRight, it is overkill to
            // fully collect the columns over the columns of the process
            // grid.
            const Int firstBlockHeight = blockSize - HRight.ColCut();
            if( winSize <= firstBlockHeight )
            {
                Matrix<F> HRightLocCopy( HRight.Matrix() );
                Gemm
                ( NORMAL, ADJOINT, F(1), ZWin, HRightLocCopy, HRight.Matrix() );
            }
            else
            {
                DistMatrix<F,STAR,MR,BLOCK> HRight_STAR_MR( HRight );
                Matrix<F> HRightLocCopy( HRight_STAR_MR.Matrix() );
                Gemm
                ( NORMAL, ADJOINT, F(1), ZWin, HRightLocCopy,
                  HRight_STAR_MR.Matrix() );
                HRight = HRight_STAR_MR;
            }
        }
        // Overwrite H(0:winBeg,winInd) *= ZWin
        auto HTop = H( IR(0,winBeg), winInd );
        // Again, we only need to overwrite HTop, so it is overkill to
        // always fully collect the rows over the rows of the grid.
        const Int firstBlockWidth = blockSize - HTop.RowCut();
        if( winSize <= firstBlockWidth )
        {
            Matrix<F> HTopLocCopy( HTop.Matrix() );
            Gemm( NORMAL, NORMAL, F(1), HTopLocCopy, ZWin, HTop.Matrix() );
        }
        else
        {
            DistMatrix<F,MC,STAR,BLOCK> HTop_MC_STAR( HTop );
            Matrix<F> HTopLocCopy( HTop_MC_STAR.Matrix() );
            Gemm
            ( NORMAL, NORMAL, F(1), HTopLocCopy, ZWin, HTop_MC_STAR.Matrix() );
            HTop = HTop_MC_STAR;
        }
    }
    if( ctrl.wantSchurVecs )
    {
        // Overwrite Z(:,winInd) *= ZWin
        auto ZBlock = Z( ALL, winInd );
        // Again, we only need to overwrite ZBlock, so it is overkill to
        // always fully collect the rows over the rows of the grid.
        const Int firstBlockWidth = blockSize - ZBlock.RowCut();
        if( winSize <= firstBlockWidth )
        {
            Matrix<F> ZBlockLocCopy( ZBlock.Matrix() );
            Gemm
            ( NORMAL, NORMAL, F(1), ZBlockLocCopy, ZWin, ZBlock.Matrix() );
        }
        else
        {
            DistMatrix<F,MC,STAR,BLOCK> ZBlock_MC_STAR( ZBlock );
            Matrix<F> ZBlockLocCopy( ZBlock_MC_STAR.Matrix() );
            Gemm
            ( NORMAL, NORMAL, F(1), ZBlockLocCopy, ZWin,
              ZBlock_MC_STAR.Matrix() );
            ZBlock = ZBlock_MC_STAR;
        }
    }

    return info;
}

} // namespace multibulge

template<typename F>
HessenbergSchurInfo
MultiBulge
( DistMatrix<F,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Base<F>>,STAR,STAR>& w,
  DistMatrix<F,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE 
    typedef Base<F> Real;
    const Real zero(0);
    const Grid& grid = H.Grid();

    const Int n = H.Height();
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
    // TODO(poulson): Re-enable this
    //minMultiBulgeSize = Max( minMultiBulgeSize, 500 );

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
    DistMatrix<F,STAR,STAR> hMainWin(grid), hSuperWin(grid);
    DistMatrix<Real,STAR,STAR> hSubWin(grid);
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
            H.Set( iterBeg, iterBeg-1, zero );
            hSubWin.Set( iterOffset-1, 0, zero );
        }
        if( iterWinSize == 1 )
        {
            w.Set( iterBeg, 0, hMainWin.GetLocal(iterOffset,0) );

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterWinSize == 2 )
        {
            const F eta00 = hMainWin.GetLocal(iterOffset,0);
            const F eta01 = hSuperWin.GetLocal(iterOffset,0);
            const Real eta10 = hSubWin.GetLocal(iterOffset,0);
            const F eta11 = hMainWin.GetLocal(iterOffset+1,0);
            multibulge::TwoByTwo
            ( H, eta00, eta01, eta10, eta11, w, Z, iterBeg, ctrl );

            winEnd = iterBeg;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterWinSize < minMultiBulgeSize )
        {
            // The window is small enough to switch to the simple scheme
            auto ctrlIter( ctrl );
            ctrlIter.winBeg = iterBeg;
            ctrlIter.winEnd = winEnd;
            auto iterInfo =
              multibulge::RedundantlyHandleWindow( H, w, Z, ctrlIter );
            info.numIterations += iterInfo.numIterations;
             
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
