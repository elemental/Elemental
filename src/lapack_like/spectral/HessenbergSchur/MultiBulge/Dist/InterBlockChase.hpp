/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_INTER_BLOCK_CHASE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_INTER_BLOCK_CHASE_HPP

#include "ComputeIntraBlockReflectors.hpp"
//#include "ComputeIntroductoryReflectors.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename F>
void InterBlockChase
(       Int numBulgesPerBlock,
        Int numIntroducedBulges,
        DistMatrix<F,MC,MR,BLOCK>& H,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
        Matrix<F>& W,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int numBulges = shifts.Height() / 2;
    const Int numRemainingBulges = numBulges - numIntroducedBulges;
    const Int numNewBulges = Min( numRemainingBulges, numBulgesPerBlock ); 

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

    vector<Matrix<F>> UList(numLocalDiagBlocks);
    if( colShift == 0 && rowShift == 0 )
    {
        // TODO(poulson)
        // Introduce the new bulges into the upper-left corner
    }
    // TODO(poulson): Perform the remaining inter-block chases
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTER_BLOCK_CHASE_HPP
