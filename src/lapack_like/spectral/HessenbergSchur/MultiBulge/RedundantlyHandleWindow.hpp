/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_REDUNDANTLY_HANDLE_WINDOW_HPP
#define EL_HESS_SCHUR_MULTIBULGE_REDUNDANTLY_HANDLE_WINDOW_HPP

#include "./Transform.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename F>
void ConsistentlyComputeDecomposition
(       DistMatrix<F,MC,MR,BLOCK>& H,
        DistMatrix<Complex<Base<F>>,STAR,STAR>& w,
        Matrix<F>& Z,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() )
{
    DEBUG_CSE
    // Because double-precision floating-point computation is often
    // non-deterministic due to extra-precision computation being frequent but
    // not guaranteed, we must be careful to not allow this non-determinism to
    // be amplified by the forward instability of Francis sweeps.
    const Grid& grid = H.Grid();
    const int owner = H.Owner(0,0);

    // Test round-trip first
    {
        DistMatrix<F,CIRC,CIRC> H_CIRC_CIRC( grid, owner );
        H_CIRC_CIRC = H;
        DistMatrix<F,MC,MR,BLOCK> HRound( grid );
        HRound = H_CIRC_CIRC;
        Output(grid.Rank(),": H ~ ",H.Height()," x ",H.Width()," ",H.ColAlign(),";",H.ColCut(),";",H.LocalHeight()," ",H.RowAlign(),";",H.RowCut(),";",H.LocalWidth(),", H_CIRC_CIRC ~ ",H_CIRC_CIRC.Height()," x ",H_CIRC_CIRC.Width()," ",H_CIRC_CIRC.Root()," ",H_CIRC_CIRC.LocalHeight()," ",H_CIRC_CIRC.LocalWidth(),", HRound ~ ",HRound.Height()," x ",HRound.Width()," ",HRound.ColAlign(),";",HRound.ColCut(),";",HRound.LocalHeight()," ",HRound.RowAlign(),";",HRound.RowCut(),";",HRound.LocalWidth());
        //Axpy( F(-1), H, HRound );
        {
            DEBUG_ONLY(AssertSameGrids( H, HRound ))
            const DistData HDistData = H.DistData();
            const DistData HRoundDistData = HRound.DistData();

            if( HDistData == HRoundDistData )
            {
                Axpy( F(-1), H.LockedMatrix(), HRound.Matrix() );
            }
            else
            {
                DistMatrix<F,MC,MR,BLOCK> HCopy( HRound.Grid(), HRound.Root() );
                //unique_ptr<BlockMatrix<T>>
                //  HCopy( HRound.Construct(HRound.Grid(),HRound.Root()) );
                HCopy.AlignWith( HRoundDistData );
                Copy( H, HCopy );
                Output(grid.Rank(),": HCopy ~ ",HCopy.Height()," x ",HCopy.Width()," ",HCopy.ColAlign(),";",HCopy.ColCut(),";",HCopy.LocalHeight()," ",HCopy.RowAlign(),";",HCopy.RowCut(),";",HCopy.LocalWidth(),", HRound ~ ",HRound.Height()," x ",HRound.Width()," ",HRound.ColAlign(),";",HRound.ColCut(),";",HRound.LocalHeight()," ",HRound.RowAlign(),";",HRound.RowCut(),";",HRound.LocalWidth());
                Axpy( F(-1), HCopy.LockedMatrix(), HRound.Matrix() );
            }
        }
        const Base<F> errNorm = FrobeniusNorm( HRound );
        if( errNorm != Base<F>(0) )
        {
            Print( H, "H" );
            Print( H_CIRC_CIRC, "H_CIRC_CIRC" );
            Print( HRound, "E" );
            LogicError("Round-trip error norm was ",errNorm);
        }
    }

    DistMatrix<F,CIRC,CIRC> H_CIRC_CIRC( grid, owner );
    H_CIRC_CIRC = H;
    w.Resize( H.Height(), 1 );
    if( H_CIRC_CIRC.CrossRank() == H_CIRC_CIRC.Root() )
        HessenbergSchur( H_CIRC_CIRC.Matrix(), w.Matrix(), Z, ctrl );
    else
        Z.Resize( H.Height(), H.Height() );
    H = H_CIRC_CIRC;
    El::Broadcast( w.Matrix(), H_CIRC_CIRC.CrossComm(), H_CIRC_CIRC.Root() );
    El::Broadcast( Z, H_CIRC_CIRC.CrossComm(), H_CIRC_CIRC.Root() );
}

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
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );

    // TODO(poulson): Fill this information structure
    HessenbergSchurInfo info;

    const auto winInd = IR(winBeg,winEnd);
    auto HWin = H(winInd,winInd);
    auto wWin = w(winInd,ALL);

    // Compute the Schur decomposition HWin = ZWin TWin ZWin',
    // where HWin is overwritten by TWin, and wWin by diag(TWin).
    Matrix<F> ZWin;
    ConsistentlyComputeDecomposition( HWin, wWin, ZWin );

    if( ctrl.fullTriangle )
    {
        if( n > winEnd )
        {
            // Overwrite H(winInd,winEnd:n) *= ZWin'
            // (applied from the left)
            auto HRight = H( winInd, IR(winEnd,n) );
            TransformRows( ZWin, HRight );
        }

        // Overwrite H(0:winBeg,winInd) *= ZWin
        auto HTop = H( IR(0,winBeg), winInd );
        TransformColumns( ZWin, HTop );
    }
    if( ctrl.wantSchurVecs )
    {
        // Overwrite Z(:,winInd) *= ZWin
        auto ZBlock = Z( ALL, winInd );
        TransformColumns( ZWin, ZBlock );
    }

    return info;
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_REDUNDANTLY_HANDLE_WINDOW_HPP
