/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_INTRA_BLOCK_CHASE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_INTRA_BLOCK_CHASE_HPP

#include "./ComputeIntraBlockReflectors.hpp"
#include "./ApplyIntraBlockReflectors.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

// Chase the separated packets of tightly-packed 4x4 bulges from the top-left
// corners of the diagonal blocks down to the bottom-right corners. This is
// accomplished by locally accumulating the reflections into a dense matrix
// and then broadcasting/allgathering said matrix within the rows and columns
// of the process grid. See Fig. 3 of 
//
//   R. Granat, Bo Kagstrom, and D. Kressner, "LAPACK Working Note #216:
//   A novel parallel QR algorithm for hybrid distributed memory HPC systems",
//
// [CITATION] for a diagram.
//

struct DistChaseState
{
  Int activeBlockBeg; // The block index of the first active block
  Int activeBlockEnd; // The block index of the end of the active blocks

  // The activeBlock{Beg,End} indices are relative to the active window implied
  // by the coupled HessenbergSchurCtrl.win{Beg,End} parameters.

  Int shiftBeg; // The first active shift index
  Int shiftEnd; // The end of the active shift indices
};

// The following extends the discussion in Granat et al. to handle active
// windows. For example, ctrl.winBeg and ctrl.winEnd define index sets
//
//   ind0 = [0,ctrl.winBeg),
//   ind1 = [ctrl.winBeg,ctrl.winEnd),
//   ind2 = [ctrl.winEnd,n),
//
// and a partitioning
//
//    H = | H00 H01 H02 |
//        | H10 H11 H12 |
//        | 0   H21 H22 |
//
// such that H11 is the active submatrix, and H10 and H21 contain a single
// nonzero entry in their upper right corners (if they are non-empty). If
// a full Schur decomposition was requested, then the appropriate pieces of H01
// and H12 must be updated by the accumulated Householder transformations
// rather than just the off-diagonal blocks of H11. Note that the fact that
// [winBeg,winEnd) may take on arbitrary values implies that we must handle
// windows which begin in the middle of distribution blocks, but thankfully this
// only effects the *inter*-block chases and not the *intra*-block chases
// (Cf. the diagrams in the *inter*-block chase source for how these
// complications are handled).
//
// Each of the intra-block multibulge chases takes a single form (though the
// last diagonal block may have a different number of bulges): the bulges are
// locally chased from the top-left to the bottom-right of the diagonal block,
// e.g., if the distribution block size was 12 and there were two bulges in the
// diagonal block, we would have the transformation
//
//         ~ ~ ~ ~ ~ ~ ~ ~ ~ ~                  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  
//      ------------------------             -------------------------
//     | B B B B x x x x x x x x |          | x x x x x x x x x x x x |
//   ~ | B B B B x x x x x x x x |        ~ | x x x x x x x x x x x x |
//   ~ | B B B B x x x x x x x x |        ~ |   x x x x x x x x x x x |
//   ~ | B B B B B B B x x x x x |        ~ |     x x x x x x x x x x |
//   ~ |       B B B B x x x x x |        ~ |       x x x x x x x x x |
//   ~ |       B B B B x x x x x |  |->   ~ |         x B B B B x x x |
//   ~ |       B B B B x x x x x |        ~ |           B B B B x x x |
//   ~ |             x x x x x x |        ~ |           B B B B x x x |
//   ~ |               x x x x x |        ~ |           B B B B B B B |
//   ~ |                 x x x x |        ~ |                 B B B B |
//   ~ |                   x x x |        ~ |                 B B B B |
//     |                     x x |          |                 B B B B |
//      -------------------------            -------------------------
//
// It is worth noting that the accumulation of the ten 3x3 Householder
// reflections for this diagram effect all but the first and last rows when
// applied from the left, and all but the first and last columns when applied
// from the right. 
//
// It is also worth noting that none of the intra-block chases involve the last
// diagonal block, as inter-block chases that introduce bulges into the last
// diagonal block are immediately chased out of the window.
//

template<typename F>
void IntraBlockChase
(       DistMatrix<F,MC,MR,BLOCK>& H,
        DistMatrix<F,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
        Matrix<F>& W,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int numShifts = shifts.Height();
    const Int numBulges = numShifts / 2;
    const Int winBeg = ( ctrl.winBeg == END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd == END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    auto& HLoc = H.Matrix();
    auto& ZLoc = Z.Matrix();
    const auto& shiftsLoc = shifts.LockedMatrix();

    // We will immediately apply the accumulated Householder transformations
    // after receiving them
    const bool immediatelyApply = true;

    // Store the local indices for bounding the off-diagonal application of the
    // accumulated Householder transformations for each diagonal block
    const Int localTransformRowBeg =
      ( ctrl.fullTriangle ? 0 : H.LocalRowOffset(winBeg) );
    const Int localTransformColEnd =
      ( ctrl.fullTriangle ? H.LocalColOffset(n) : H.LocalColOffset(winEnd) );

    const Grid& grid = H.Grid();
    const int gridHeight = grid.Height();
    const int gridWidth = grid.Width();
    const int gridRow = grid.Row();
    const int gridCol = grid.Col();

    const Int blockSize = H.BlockHeight();
    if( blockSize != H.BlockWidth() )
        LogicError("IntraBlockChase assumes square distribution blocks"); 
    if( winSize < 2 * blockSize )
        LogicError("The window size must be at least twice the block size");
    const Int blockCut = H.ColCut();
    if( blockCut != H.RowCut() )
        LogicError("IntraBlockChase assumes that the cuts are equal");
    if( ctrl.wantSchurVecs )
    {
        // Ensure that H and Z have the same distributions
        if( Z.DistData() != H.DistData() )
            LogicError("The distributions of H and Z should match");
    }

    // Apply the black-box function for determining the number of bulges
    const Int numBulgesPerBlock = ctrl.numBulgesPerBlock(blockSize);
    if( numBulgesPerBlock*6 + 1 > blockSize )
        LogicError
        ("Cannot incorporate ",numBulgesPerBlock,
         " bulges into a distributed block size of ",blockSize);

    // We push the leftover shifts through first, so the last block only has
    // the leftover number of shifts if shiftEnd is equal to the number of
    // shifts.
    const Int numBulgesInLastBlock =
      ( state.shiftEnd == numShifts ?
        Mod(numBulges,numBulgesPerBlock) :
        numBulgesPerBlock );

    // Compute the offset into the distribution block that the window begins in
    const Int windowCut = Mod( blockCut+winBeg, blockSize );

    // Compute the remaining size of the first block in the window
    const Int firstBlockSize = blockSize - windowCut;

    // Compute the row and column owners of the beginning of the window.
    const int winColAlign = H.RowOwner( winBeg );
    const int winRowAlign = H.ColOwner( winBeg );

    // Compute the row and column owners of the beginning of the active bulge
    // region and the first block assigned to us.
    const Int activeBeg = winBeg +
      ( state.activeBlockBeg == 0 ?
        0 :
        firstBlockSize + (state.activeBlockBeg-1)*blockSize );
    const int activeColAlign = H.RowOwner( activeBeg );
    const int activeRowAlign = H.ColOwner( activeBeg );
    const int activeColShift = Shift( gridRow, activeColAlign, gridHeight );
    const int activeRowShift = Shift( gridCol, activeRowAlign, gridWidth );
    const Int activeRowBlockBeg = state.activeBlockBeg + activeColShift;
    const Int activeColBlockBeg = state.activeBlockBeg + activeRowShift;

    // Count the number of diagonal blocks assigned to this process
    // TODO(poulson): Move this into a subroutine?
    Int numLocalDiagBlocks;
    {
        // Only loop over the row blocks that are assigned to our process row
        // and occur within the active window.
        Int diagBlock = activeRowBlockBeg;
        Int localDiagBlock = 0;
        while( diagBlock < state.activeBlockEnd )
        {
            const int ownerCol = (winRowAlign + diagBlock) % gridWidth;
            if( ownerCol == gridCol )
                ++localDiagBlock;
            diagBlock += gridHeight;
        }
        numLocalDiagBlocks = localDiagBlock;
    }
    
    // Chase bulges down the local diagonal blocks and store the accumulations
    // of the Householder reflections.
    vector<Matrix<F>> UList(numLocalDiagBlocks);
    {
        // Only loop over the row blocks that are assigned to our process row
        Int diagBlock = activeRowBlockBeg;
        Int localDiagBlock = 0;
        Zeros( W, 3, numBulgesPerBlock );
        while( diagBlock < state.activeBlockEnd )
        {
            const Int thisBlockHeight = 
              ( diagBlock == 0 ? firstBlockSize : blockSize );
            const Int numBlockBulges =
              ( diagBlock==state.activeBlockEnd-1 ?
                numBulgesInLastBlock :
                numBulgesPerBlock );

            const int ownerCol = (winRowAlign + diagBlock) % gridWidth;
            if( ownerCol == gridCol )
            {
                const Int diagOffset = winBeg +
                  (diagBlock == 0 ?
                   0 :
                   firstBlockSize + (diagBlock-1)*blockSize );

                // View the local diagonal block of H
                const Int localRowOffset = H.LocalRowOffset( diagOffset );
                const Int localColOffset = H.LocalColOffset( diagOffset );
                auto HBlockLoc = 
                  HLoc
                  ( IR(0,thisBlockHeight)+localRowOffset,
                    IR(0,thisBlockHeight)+localColOffset );

                // View the local shifts for this diagonal block
                const Int shiftOffset = state.shiftBeg +
                  (2*numBulgesPerBlock)*(diagBlock-state.activeBlockBeg);
                auto shiftsBlockLoc =
                  shiftsLoc( IR(0,2*numBlockBulges)+shiftOffset, ALL );

                // Initialize the accumulated reflection matrix
                auto& UBlock = UList[localDiagBlock];
                Identity( UBlock, thisBlockHeight-1, thisBlockHeight-1 );

                // Perform the diagonal block sweep and accumulate the
                // reflections in UBlock
                const Int numSteps = (thisBlockHeight-1) - 3*numBlockBulges;
                for( Int step=0; step<numSteps; ++step )
                {
                    ComputeIntraBlockReflectors
                    ( step, numBlockBulges, HBlockLoc, shiftsLoc, W,
                      ctrl.progress );
                    ApplyIntraBlockReflectorsOpt
                    ( step, numBlockBulges, HBlockLoc, UBlock, W,
                      ctrl.progress );
                }
                ++localDiagBlock;
            }
            diagBlock += gridHeight;
        }
    }

    // Broadcast/Allgather the accumulated reflections within rows
    // TODO(poulson): Move this into a subroutine
    if( immediatelyApply )
    {
        Matrix<F> UBlock;
        Matrix<F> tempMatrix;

        // Only loop over the row blocks assigned to this grid row
        Int diagBlock = activeRowBlockBeg;
        Int localDiagBlock = 0;
        while( diagBlock < state.activeBlockEnd )
        {
            const Int thisBlockHeight = 
              ( diagBlock == 0 ? firstBlockSize : blockSize );
            const Int numBlockBulges =
              ( diagBlock==state.activeBlockEnd-1 ?
                numBulgesInLastBlock :
                numBulgesPerBlock );

            const int ownerCol = (winRowAlign + diagBlock) % gridWidth;
            if( ownerCol == gridCol )
            {
                UBlock = UList[localDiagBlock];
                ++localDiagBlock;
            }

            const Int diagOffset = winBeg +
              (diagBlock == 0 ?
               0 :
               firstBlockSize + (diagBlock-1)*blockSize );
            const Int localRowOffset = H.LocalRowOffset( diagOffset );
            const Int localColOffset = H.LocalColOffset( diagOffset );
            El::Broadcast( UBlock, H.RowComm(), ownerCol );
            const auto applyRowInd = IR(1,thisBlockHeight-1) + localRowOffset;
            const auto applyColInd =
              IR(localColOffset+thisBlockHeight,localTransformColEnd);

            auto HLocRight = HLoc( applyRowInd, applyColInd );
            tempMatrix = HLocRight;
            Gemm( ADJOINT, NORMAL, F(1), UBlock, tempMatrix, HLocRight );

            diagBlock += gridHeight;
        }
    }
    else
    {
        // TODO(poulson): Add support for AllGather variant
        LogicError("This option is not yet supported");
    }

    // Broadcast/Allgather the accumulated reflections within columns
    // TODO(poulson): Move this into a subroutine
    if( immediatelyApply )
    {
        Matrix<F> UBlock;
        Matrix<F> tempMatrix;

        // Only loop over the row blocks assigned to this grid row
        Int diagBlock = activeColBlockBeg;
        Int localDiagBlock = 0;
        while( diagBlock < state.activeBlockEnd )
        {
            const Int thisBlockHeight = 
              ( diagBlock == 0 ? firstBlockSize : blockSize );
            const Int numBlockBulges =
              ( diagBlock==state.activeBlockEnd-1 ?
                numBulgesInLastBlock :
                numBulgesPerBlock );

            const int ownerRow = (winColAlign + diagBlock) % gridHeight;
            if( ownerRow == gridRow )
            {
                UBlock = UList[localDiagBlock];
                ++localDiagBlock;
            }

            const Int diagOffset = winBeg +
              ( diagBlock == 0 ?
                0 :
                firstBlockSize + (diagBlock-1)*blockSize );
            const Int localRowOffset = H.LocalRowOffset( diagOffset );
            const Int localColOffset = H.LocalColOffset( diagOffset );
            El::Broadcast( UBlock, H.ColComm(), ownerRow );
            const auto applyRowInd = IR(localTransformRowBeg,localRowOffset);
            const auto applyColInd = IR(1,thisBlockHeight-1) + localColOffset;

            auto HLocAbove = HLoc( applyRowInd, applyColInd );
            tempMatrix = HLocAbove;
            Gemm( NORMAL, NORMAL, F(1), tempMatrix, UBlock, HLocAbove );
            if( ctrl.wantSchurVecs )
            {
                auto ZLocBlock = ZLoc( ALL, applyColInd ); 
                tempMatrix = ZLocBlock;
                Gemm( NORMAL, NORMAL, F(1), tempMatrix, UBlock, ZLocBlock );
            }

            diagBlock += gridWidth;
        }
    }
    else
    {
        // TODO(poulson): Add support for AllGather variant
        LogicError("This option is not yet supported");
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTRA_BLOCK_CHASE_HPP
