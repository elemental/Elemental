/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_INTRA_BLOCK_CHASE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_INTRA_BLOCK_CHASE_HPP

#include "./PairShifts.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

// We keep track of the number of diagonal blocks of tightly-coupled shifts
// that have been chased down the diagonal so far.
//
// If numChasedBlocks is zero, then the Min(numBulgesPerBlock,numShifts/2)
// bulges have already been introduced in the top-left block.
// 
template<typename F>
void IntraBlockChase
(       Int numBulgesPerBlock,
        Int numChasedBlocks,
        DistMatrix<F,MC,MR,BLOCK>& H,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
        Matrix<F>& Z,
        Matrix<F>& U,
        Matrix<F>& W,
        Matrix<F>& WAccum,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int numShifts = shifts.Height();
    const Int blockHeight = H.BlockHeight();
    const Int blockWidth = H.BlockWidth();
    if( blockHeight != blockWidth )
        LogicError("IntraBlockChase assumes square distribution blocks"); 
    const Int colCut = H.ColCut();
    const Int rowCut = H.RowCut();
    if( colCut != 0 || rowCut != 0 )
        LogicError("IntraBlockChase assumes zero row and column cuts");

    const Int colShift = H.ColShift();
    const Int rowShift = H.RowShift();
    const Int colStride = H.ColStride();
    const Int rowStride = H.RowStride();
    // HERE
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTRA_BLOCK_CHASE_HPP
