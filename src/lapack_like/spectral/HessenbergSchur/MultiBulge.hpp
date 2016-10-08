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
            multibulge::TwoByTwo( H, w, Z, iterBeg, ctrl );
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

    // TODO(poulson): Implement a more reasonable/configurable means of deciding
    // when to call the sequential implementation
    Int minMultiBulgeSize = Max( ctrl.minMultiBulgeSize, 2*H.BlockSize() );
    // This maximum is meant to account for parallel overheads and needs to be
    // more principled (and perhaps based upon the number of workers and the 
    // cluster characteristics)
    minMultiBulgeSize = Max( minMultiBulgeSize, 500 );

    HessenbergSchurInfo info;

    w.Resize( n, 1 );
    if( winSize < minMultiBulgeSize )
    {
        DistMatrix<F,STAR,STAR> HFull( H );
        DistMatrix<F,STAR,STAR> ZFull( H.Grid() );
        if( ctrl.wantSchurVecs )
        {
            ZFull = Z;
            // In case Z was not n x n, ensure that ZFull is
            ZFull.Resize( n, n );
        }

        auto info = Simple( HFull.Matrix(), w.Matrix(), ZFull.Matrix(), ctrl );

        H = HFull;
        if( ctrl.wantSchurVecs )
            Z = ZFull;
        return info;
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
        multibulge::GatherTridiagonal
        ( H, winInd, hMainWin, hSubWin, hSuperWin );

        const Int iterOffset =
          DetectSmallSubdiagonal( hMainWin, hSubWin, hSuperWin );
        const Int iterBeg = winEnd + iterOffset;
        if( iterOffset > 0 )
        {
            H.Set( iterBeg, iterBeg-1, zero );
            hSubWin.Set( iterOffset, 0, zero );
        }
        if( iterBeg == winEnd-1 )
        {
            w.Set( iterBeg, 0, hMainWin.GetLocal(iterOffset,0) );
            --winEnd;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterBeg == winEnd-2 )
        {
            const Real eta00 = hMainWin.GetLocal(iterOffset,0);
            const Real eta01 = hSuperWin.GetLocal(iterOffset,0);
            const Real eta10 = hSubWin.GetLocal(iterOffset,0);
            const Real eta11 = hMainWin.GetLocal(iterOffset+1,0);
            multibulge::TwoByTwo
            ( H, eta00, eta01, eta10, eta11, Z, iterBeg, ctrl );
            winEnd -= 2;
            numIterSinceDeflation = 0;
            continue;
        }

        const Int iterWinSize = winEnd-iterBeg;
        if( iterWinSize < minMultiBulgeSize )
        {
            // The window is small enough to switch to the simple scheme
            //
            // TODO(poulson): Have all workers independently process the
            // diagonal block and then apply the transformation into the
            // remainder of the matrix.
            LogicError("This section is not yet finished");
            /*
            auto ctrlSub( ctrl );
            ctrlSub.winBeg = iterBeg;
            ctrlSub.winEnd = winEnd;
            Simple( H, w, Z, ctrlSub );
            */
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
