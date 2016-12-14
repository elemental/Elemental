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

template<typename Field>
void ConsistentlyComputeDecomposition
(       DistMatrix<Field,MC,MR,BLOCK>& H,
        DistMatrix<Complex<Base<Field>>,STAR,STAR>& w,
        Matrix<Field>& Z,
  const HessenbergSchurCtrl& ctrl=HessenbergSchurCtrl() )
{
    EL_DEBUG_CSE
    // Because double-precision floating-point computation is often
    // non-deterministic due to extra-precision computation being frequent but
    // not guaranteed, we must be careful to not allow this non-determinism to
    // be amplified by the forward instability of Francis sweeps.
    const Grid& grid = H.Grid();
    const int owner = H.Owner(0,0);

    DistMatrix<Field,CIRC,CIRC> H_CIRC_CIRC( grid, owner );
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

template<typename Field>
HessenbergSchurInfo
RedundantlyHandleWindow
( DistMatrix<Field,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Base<Field>>,STAR,STAR>& w,
  DistMatrix<Field,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
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
    Matrix<Field> ZWin;
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
