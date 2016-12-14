/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTRA_BLOCK_CHASE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTRA_BLOCK_CHASE_HPP

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

namespace intrablock {

// Form the list of accumulated Householder transformations, which should be
// applied as
//
//     \hat{U}_i' H_i \hat{U}_i,
//
// where \hat{U}_i = | 1, 0,   0 |
//                   | 0, U_i, 0 |
//                   | 0, 0,   1 |
//
// is the extension of U_i to the entire diagonal block (as the transformation
// leaves the first and last rows unchanged when applied from the left), and
// H_i is the i'th locally-owned diagonal block of H.
//
template<typename Field>
void LocalChase
(       DistMatrix<Field,MC,MR,BLOCK>& H,
  const DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl,
        vector<Matrix<Field>>& UList )
{
    EL_DEBUG_CSE
    const Grid& grid = H.Grid();

    Matrix<Field> W;
    auto& HLoc = H.Matrix();
    const auto& shiftsLoc = shifts.LockedMatrix();

    // Packets are originally pushed directly into the second block if the first
    // is not full.
    const Int intraBlockStart =
      ( state.firstBlockSize == state.blockSize ?
        state.introBlock+1 : Max(state.introBlock+1,1) );

    // Count the number of diagonal blocks assigned to this process
    {
        // Only loop over the row blocks that are assigned to our process row
        // and occur within the active window.
        Int diagBlock = state.activeRowBlockBeg;
        while( diagBlock < intraBlockStart )
            diagBlock += grid.Height();

        Int numLocalBlocks = 0;
        // Recall that packets are never left in the last block of the window
        while( diagBlock < Min(state.endBlock,state.numWinBlocks-1) )
        {
            const int ownerCol =
              Mod( state.winRowAlign+diagBlock, grid.Width() );
            if( ownerCol == grid.Col() )
                ++numLocalBlocks;
            diagBlock += grid.Height();
        }
        UList.resize(numLocalBlocks);
    }

    // Chase bulges down the local diagonal blocks and store the accumulations
    // of the Householder reflections. We only loop over the row blocks that
    // are assigned to our process row and filter based upon whether or not
    // we are in the correct process column.
    Int localDiagBlock = 0;
    Int diagBlock = state.activeRowBlockBeg;
    while( diagBlock < intraBlockStart )
        diagBlock += grid.Height();
    Zeros( W, 3, state.numBulgesPerBlock );
    Matrix<Field> ZDummy;
    while( diagBlock < Min(state.endBlock,state.numWinBlocks-1) )
    {
        const Int numBlockBulges =
          ( diagBlock==state.endBlock-1 ?
            state.numBulgesInLastBlock :
            state.numBulgesPerBlock );

        const int ownerCol = Mod( state.winRowAlign+diagBlock, grid.Width() );
        if( ownerCol == grid.Col() )
        {
            const Int diagOffset = state.winBeg +
              ( diagBlock == 0 ?
                0 :
                state.firstBlockSize + (diagBlock-1)*state.blockSize );

            // View the local diagonal block of H
            const Int localRowOffset = H.LocalRowOffset( diagOffset );
            const Int localColOffset = H.LocalColOffset( diagOffset );
            auto HBlockLoc =
              HLoc
              ( IR(0,state.blockSize)+localRowOffset,
                IR(0,state.blockSize)+localColOffset );

            // View the local shifts for this diagonal block
            const Int bulgeOffset = state.bulgeBeg +
              state.numBulgesPerBlock*(diagBlock-intraBlockStart);
            auto packetShifts =
              shiftsLoc( IR(0,2*numBlockBulges)+(2*bulgeOffset), ALL );

            // Initialize the accumulated reflection matrix; recall that it
            // does not effect the first or last index of the block. For
            // example, consider the effects of a single 3x3 Householder
            // similarity bulge chase step
            //
            //        ~ ~ ~                 ~ ~ ~
            //     -----------           -----------
            //    | B B B B x |  |->    | x x x x x |
            //  ~ | B B B B x |       ~ | x B B B B |
            //  ~ | B B B B x |       ~ |   B B B B |.
            //  ~ | B B B B x |       ~ |   B B B B |
            //    |       x x |         |   B B B B |
            //     -----------           -----------
            //
            auto& UBlock = UList[localDiagBlock];
            Identity( UBlock, state.blockSize-2, state.blockSize-2 );

            // Perform the diagonal block sweep and accumulate the
            // reflections in UBlock. The number of diagonal entries spanned
            // by numBlockBulges bulges is 1 + 3 numBlockBulges, so the number
            // of steps is blockSize - (1 + 3*numBlockBulges).
            const Int numSteps = state.blockSize - (1 + 3*numBlockBulges);
            const Int blockWinBeg = 0;
            const Int blockWinEnd = state.blockSize;
            const Int chaseBeg = 0;
            const Int transformRowBeg = 0;
            const Int transformColEnd = state.blockSize;
            const bool wantSchurVecsSub = false;
            const bool accumulateSub = true;
            const Int firstBulge = 0;
            for( Int step=0; step<numSteps; ++step )
            {
                const Int packetBeg = step;
                ComputeReflectors
                ( HBlockLoc, blockWinBeg, blockWinEnd, packetShifts, W,
                  packetBeg, firstBulge, numBlockBulges, ctrl.progress );
                ApplyReflectorsOpt
                ( HBlockLoc, blockWinBeg, blockWinEnd, chaseBeg, packetBeg,
                  transformRowBeg, transformColEnd, ZDummy, wantSchurVecsSub,
                  UBlock, W, firstBulge, numBlockBulges, accumulateSub,
                  ctrl.progress );
            }
            ++localDiagBlock;
        }
        diagBlock += grid.Height();
    }
}

template<typename Field>
void ApplyAccumulatedReflections
(       DistMatrix<Field,MC,MR,BLOCK>& H,
        DistMatrix<Field,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl,
  const vector<Matrix<Field>>& UList )
{
    EL_DEBUG_CSE
    const Grid& grid = H.Grid();

    // We will immediately apply the accumulated Householder transformations
    // after receiving them
    const bool immediatelyApply = true;

    auto& HLoc = H.Matrix();
    auto& ZLoc = Z.Matrix();

    // Packets are originally pushed directly into the second block if the first
    // is not full.
    const Int intraBlockStart =
      ( state.firstBlockSize == state.blockSize ?
        state.introBlock+1 : Max(state.introBlock+1,1) );

    // Broadcast/Allgather the accumulated reflections within rows
    if( immediatelyApply )
    {
        Matrix<Field> UBlock;
        Matrix<Field> tempMatrix;

        // Only loop over the row blocks assigned to this grid row
        Int diagBlockRow = state.activeRowBlockBeg;
        while( diagBlockRow < intraBlockStart )
            diagBlockRow += grid.Height();

        Int localDiagBlock = 0;
        while( diagBlockRow < Min(state.endBlock,state.numWinBlocks-1) )
        {
            const int ownerCol =
              Mod( state.winRowAlign+diagBlockRow, grid.Width() );
            if( ownerCol == grid.Col() )
            {
                UBlock = UList[localDiagBlock++];
            }
            else
                Zeros( UBlock, state.blockSize-2, state.blockSize-2 );
            if( UBlock.Height() != state.blockSize-2 ||
                UBlock.Width() != state.blockSize-2 )
                LogicError
                ("UBlock was ",UBlock.Height()," x ",UBlock.Width(),
                 " instead of ",state.blockSize-2," x ",state.blockSize-2);
            El::Broadcast( UBlock, H.RowComm(), ownerCol );

            // TODO(poulson): Move into subroutine
            const Int diagOffset = state.winBeg +
              ( diagBlockRow == 0 ?
                0 :
                state.firstBlockSize + (diagBlockRow-1)*state.blockSize );
            const Int localRowOffset = H.LocalRowOffset( diagOffset );
            const Int localColOffset =
              H.LocalColOffset( diagOffset+state.blockSize );
            const auto applyRowInd = IR(1,state.blockSize-1) + localRowOffset;
            const auto applyColInd =
              IR(localColOffset,state.localTransformColEnd);

            auto HLocRight = HLoc( applyRowInd, applyColInd );
            tempMatrix = HLocRight;
            Gemm( ADJOINT, NORMAL, Field(1), UBlock, tempMatrix, HLocRight );

            diagBlockRow += grid.Height();
        }
    }
    else
    {
        // TODO(poulson): Add support for AllGather variant
        LogicError("This option is not yet supported");
    }

    // Broadcast/Allgather the accumulated reflections within columns
    if( immediatelyApply )
    {
        Matrix<Field> UBlock;
        Matrix<Field> tempMatrix;

        // Only loop over the row blocks assigned to this grid row
        Int diagBlockCol = state.activeColBlockBeg;
        while( diagBlockCol < intraBlockStart )
            diagBlockCol += grid.Width();

        Int localDiagBlock = 0;
        while( diagBlockCol < Min(state.endBlock,state.numWinBlocks-1) )
        {
            const int ownerRow =
              Mod( state.winColAlign+diagBlockCol, grid.Height() );
            if( ownerRow == grid.Row() )
            {
                UBlock = UList[localDiagBlock++];
            }
            else
                Zeros( UBlock, state.blockSize-2, state.blockSize-2 );
            if( UBlock.Height() != state.blockSize-2 ||
                UBlock.Width() != state.blockSize-2 )
                LogicError
                ("UBlock was ",UBlock.Height()," x ",UBlock.Width(),
                 " instead of ",state.blockSize-2," x ",state.blockSize-2);
            El::Broadcast( UBlock, H.ColComm(), ownerRow );

            // TODO(poulson): Move into subroutine
            const Int diagOffset = state.winBeg +
              ( diagBlockCol == 0 ?
                0 :
                state.firstBlockSize + (diagBlockCol-1)*state.blockSize );
            const Int localRowOffset = H.LocalRowOffset( diagOffset );
            const Int localColOffset = H.LocalColOffset( diagOffset );
            const auto applyRowInd =
              IR(state.localTransformRowBeg,localRowOffset);
            const auto applyColInd = IR(1,state.blockSize-1) + localColOffset;

            auto HLocAbove = HLoc( applyRowInd, applyColInd );
            tempMatrix = HLocAbove;
            Gemm( NORMAL, NORMAL, Field(1), tempMatrix, UBlock, HLocAbove );
            if( ctrl.wantSchurVecs )
            {
                auto ZLocBlock = ZLoc( ALL, applyColInd );
                tempMatrix = ZLocBlock;
                Gemm( NORMAL, NORMAL, Field(1), tempMatrix, UBlock, ZLocBlock );
            }

            diagBlockCol += grid.Width();
        }
    }
    else
    {
        // TODO(poulson): Add support for AllGather variant
        LogicError("This option is not yet supported");
    }
}

} // namespace intrablock

template<typename Field>
void IntraBlockChase
(       DistMatrix<Field,MC,MR,BLOCK>& H,
        DistMatrix<Field,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE

    // Form the list of accumulated Householder transformations, which should be
    // applied as
    //
    //     \hat{U}_i' H_i \hat{U}_i,
    //
    // where \hat{U}_i = | 1, 0,   0 |
    //                   | 0, U_i, 0 |
    //                   | 0, 0,   1 |
    //
    // is the extension of U_i to the entire diagonal block (as the
    // transformation leaves the first and last rows unchanged when applied from
    // the left), and H_i is the i'th locally-owned diagonal block of H.
    //
    vector<Matrix<Field>> UList;
    intrablock::LocalChase( H, shifts, state, ctrl, UList );

    // Broadcast the accumulated transformations from the owning diagonal block
    // over the entire process row/column teams and then apply them to Z from
    // the right (if the Schur vectors are desired), to the above-diagonal
    // portion of H from the right, and their adjoints to the relevant
    // right-of-diagonal portions of H from the left.
    intrablock::ApplyAccumulatedReflections( H, Z, shifts, state, ctrl, UList );
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTRA_BLOCK_CHASE_HPP
