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

#include "IntraBlockChase/ComputeReflectors.hpp"

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
//
// We keep track of the number of diagonal blocks of tightly-coupled shifts
// that have been chased down the diagonal so far.
//
// If numChasedBlocks is zero, then the Min(numBulgesPerBlock,numShifts/2)
// bulges have already been introduced in the top-left block.
// 
template<typename F>
void IntraBlockChase
(       Int numBulgesPerBlock,
        Int numIntroducedBulges,
        DistMatrix<F,MC,MR,BLOCK>& H,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
        Matrix<F>& Z,
        Matrix<F>& U,
        Matrix<F>& W,
        Matrix<F>& WAccum,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int numShifts = shifts.Height();

    const Grid& grid = H.Grid();
    const int gridHeight = grid.Height();
    const int gridWidth = grid.Width();
    const int gridSize = gridHeight * gridWidth;
    const int gridRow = grid.Row();
    const int gridCol = grid.Col();

    const Int blockHeight = H.BlockHeight();
    const Int blockWidth = H.BlockWidth();
    if( blockHeight != blockWidth )
        LogicError("IntraBlockChase assumes square distribution blocks"); 
    const Int colCut = H.ColCut();
    const Int rowCut = H.RowCut();
    if( colCut != 0 || rowCut != 0 )
        LogicError("IntraBlockChase assumes zero row and column cuts");

    const int colAlign = H.ColAlign();
    const int rowAlign = H.RowAlign();
    const int colShift = H.ColShift();
    const int rowShift = H.RowShift();

    // Count the number of diagonal blocks assigned to this process
    // TODO(poulson): Move this into a subroutine
    Int numLocalDiagBlocks;
    {
        // Only loop over the row blocks that are assigned to our process row
        Int diagBlock = colShift;
        Int diagOffset = colShift * blockHeight;
        Int bulgeOffset = colShift * numBulgesPerBlock;
        Int localDiagBlock = 0;
        while( diagOffset < n && bulgeOffset < numIntroducedBulges )
        {
            const int ownerCol = (rowAlign + diagBlock) % gridWidth;
            if( ownerCol == gridCol )
                ++localDiagBlock;

            diagBlock += gridHeight;
            diagOffset += gridHeight * blockHeight;
            bulgeOffset += gridHeight * numBulgesPerBlock;
        }
        numLocalDiagBlocks = localDiagBlock;
    }
    
    // Chase bulges down the local diagonal blocks and store the accumulations
    // of the Householder reflections.
    vector<Matrix<F>> UList(numLocalDiagBlocks);
    {
        // Only loop over the row blocks that are assigned to our process row
        Int diagBlock = colShift;
        Int diagOffset = colShift * blockHeight;
        Int bulgeOffset = colShift * numBulgesPerBlock;
        Int localDiagBlock = 0;
        auto& HLoc = H.Matrix();
        const auto& shiftsLoc = shifts.LockedMatrix();
        Matrix<F> W( 3, numBulgesPerBlock );
        while( diagOffset < n && bulgeOffset < numIntroducedBulges )
        {
            const Int thisBlockHeight = Min( blockHeight, n-diagOffset );
            const Int numBlockBulges =
              Min( numBulgesPerBlock, numIntroducedBulges-bulgeOffset );

            const int ownerCol = (rowAlign + diagBlock) % gridWidth;
            if( ownerCol == gridCol )
            {
                const Int localRowOffset = H.LocalRowOffset( diagOffset );
                const Int localColOffset = H.LocalColOffset( diagOffset );
                auto HBlockLoc = 
                  HLoc
                  ( IR(0,thisBlockHeight)+localRowOffset,
                    IR(0,thisBlockHeight)+localColOffset );
                auto shiftsBlockLoc =
                  shiftsLoc( IR(0,2*numBlockBulges)+2*bulgeOffset, ALL );
                auto& UBlock = UList[localDiagBlock];
                Identity( UBlock, thisBlockHeight-1, thisBlockHeight-1 );
                const Int numSteps = (thisBlockHeight-1) - 3*numBlockBulges;
                for( Int step=0; step<numSteps; ++step )
                {
                    ComputeDiagonalBlockReflectors
                    ( step, numBlockBulges, HBlockLoc, shiftsLoc, W,
                      ctrl.progress );
                    ApplyReflectorsToDiagonalBlockOpt
                    ( step, numBlockBulges, HBlockLoc, UBlock, W,
                      ctrl.progress );
                }
                ++localDiagBlock;
            }

            diagBlock += gridHeight;
            diagOffset += gridHeight * blockHeight;
            bulgeOffset += gridHeight * numBulgesPerBlock;
        }
    }

    // Broadcast/Allgather the accumulated reflections within rows
    // TODO(poulson): Move this into a subroutine
    const bool immediatelyApply = true;
    if( immediatelyApply )
    {
        Matrix<F> UBlock;
        Matrix<F> HLocRightCopy;

        // Only loop over the row blocks assigned to this grid row
        Int diagBlock = colShift;
        Int diagOffset = colShift * blockHeight;
        Int bulgeOffset = colShift * numBulgesPerBlock;
        Int localDiagBlock = 0;
        auto& HLoc = H.Matrix();
        while( diagOffset < n && bulgeOffset < numIntroducedBulges )
        {
            const Int thisBlockHeight = Min( blockHeight, n-diagOffset );
            const Int numBlockBulges =
              Min( numBulgesPerBlock, numIntroducedBulges-bulgeOffset );

            const int ownerCol = (rowAlign + diagBlock) % gridWidth;
            if( ownerCol == gridCol )
            {
                UBlock = UList[localDiagBlock];
                ++localDiagBlock;
            }

            const Int localRowOffset = H.LocalRowOffset( diagOffset );
            const Int localColOffset = H.LocalColOffset( diagOffset );
            El::Broadcast( UBlock, H.RowComm(), ownerCol );
            const auto applyRowInd = IR(1,thisBlockHeight-1) + localRowOffset;
            const auto applyColInd = IR(localColOffset+thisBlockHeight,END);
            auto HLocRight = HLoc( applyRowInd, applyColInd );
            HLocRightCopy = HLocRight;
            Gemm( NORMAL, NORMAL, F(1), UBlock, HLocRightCopy, HLocRight );

            diagBlock += gridHeight;
            diagOffset += gridHeight * blockHeight;
            bulgeOffset += gridHeight * numBulgesPerBlock;
        }
    }
    else
    {
        // TODO(poulson): Add support for AllGather variant
    }

    // Broadcast/Allgather the accumulated reflections within columns
    // TODO(poulson): Move this into a subroutine
    if( immediatelyApply )
    {
        Matrix<F> UBlock;
        Matrix<F> HLocAboveCopy;

        // Only loop over the column blocks assigned to this grid column
        Int diagBlock = rowShift;
        Int diagOffset = rowShift * blockHeight;
        Int bulgeOffset = rowShift * numBulgesPerBlock;
        Int localDiagBlock = 0;
        auto& HLoc = H.Matrix();
        while( diagOffset < n && bulgeOffset < numIntroducedBulges )
        {
            const Int thisBlockHeight = Min( blockHeight, n-diagOffset );
            const Int numBlockBulges =
              Min( numBulgesPerBlock, numIntroducedBulges-bulgeOffset );

            const int ownerRow = (colAlign + diagBlock) % gridHeight;
            if( ownerRow == gridRow )
            {
                UBlock = UList[localDiagBlock];
                ++localDiagBlock;
            }

            const Int localRowOffset = H.LocalRowOffset( diagOffset );
            const Int localColOffset = H.LocalColOffset( diagOffset );
            El::Broadcast( UBlock, H.ColComm(), ownerRow );
            const auto applyRowInd = IR(0,localRowOffset);
            const auto applyColInd = IR(1,thisBlockHeight-1) + localColOffset;
            auto HLocAbove = HLoc( applyRowInd, applyColInd );
            HLocAboveCopy = HLocAbove;
            Gemm( NORMAL, ADJOINT, F(1), HLocAboveCopy, UBlock, HLocAbove );

            diagBlock += gridWidth;
            diagOffset += gridWidth * blockHeight;
            bulgeOffset += gridWidth * numBulgesPerBlock;
        }
    }
    else
    {
        // TODO(poulson): Add support for AllGather variant
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTRA_BLOCK_CHASE_HPP
