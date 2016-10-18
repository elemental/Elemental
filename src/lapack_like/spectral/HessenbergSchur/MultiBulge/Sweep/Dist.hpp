/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_DIST_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_DIST_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

// A structure for storing auxiliary metadata that is somewhat tedious to 
// repeatedly recompute
struct DistChaseState
{
  // The following is fixed throughout a sweep
  // -----------------------------------------

  // These are the same as ctrl.win{Beg,End} but with END properly replaced
  Int winBeg;
  Int winEnd;

  // The process rows and columns owning the block at the top-left of the window
  int winRowAlign;
  int winColAlign;

  Int blockSize;

  // The size of the piece of the bottom-right quadrant of the distribution
  // block that the window starts at. For example, if the distribution block
  // size is 16 but the window starts two diagonals down the diagonal of such 
  // a block, then 'firstBlockSize' is 14.
  Int firstBlockSize;
  // The size of the last block in the window.
  Int lastBlockSize;
  // The number of (full or partial) diagonal blocks in the window.
  Int numWinBlocks;

  // The first local row index of H that should be updated
  Int localTransformRowBeg;
  // The last local column index of H that should be updated
  Int localTransformColEnd;

  // The number of bulges in the vast majority of the packets
  Int numBulgesPerBlock;

  // The following changes throughout a sweep as the active window changes
  // ---------------------------------------------------------------------
  Int bulgeBeg; // The first active bulge index
  Int bulgeEnd; // The (non-inclusive) end of the active bulge indices

  // The number of bulges currently in the very last active block
  Int numBulgesInLastBlock;

  // The index of the diagonal entry that the first active shift is in. If 
  // activeBlockBeg is zero, then this is the top-left entry of the bottom-right
  // quadrant of the distribution block lying within the window; otherwise, it
  // is the index of the top-left entry of the full distribution block.
  Int activeBeg;
  Int activeEnd;

  // The activeBlock{Beg,End} indices are relative to the active window implied
  // by win{Beg,End}.
  Int activeBlockBeg; // The block index of the first active block
  Int activeBlockEnd; // The block index of the end of the active blocks

  // The first diagonal block assigned to our process row that occurs on or
  // after activeBlockBeg
  Int activeRowBlockBeg;
  // The first diagonal block assigned to our process column that occurs on or
  // after activeBlockBeg
  Int activeColBlockBeg; 
};

template<typename F>
DistChaseState BuildDistChaseState
( const DistMatrix<F,MC,MR,BLOCK>& H,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int numBulges = shifts.Height() / 2;
    const Grid& grid = H.Grid();

    DistChaseState state;

    state.winBeg = ( ctrl.winBeg == END ? n : ctrl.winBeg );
    state.winEnd = ( ctrl.winEnd == END ? n : ctrl.winEnd );
    state.winRowAlign = H.ColOwner( state.winBeg );
    state.winColAlign = H.RowOwner( state.winBeg );

    state.blockSize = H.BlockHeight();

    // Compute the remaining size of the first block in the window
    const Int winCut = Mod( state.winBeg + H.ColCut(), state.blockSize );
    state.firstBlockSize = state.blockSize - winCut;

    // Compute the number of (full or partial) diagonal blocks in the window 
    const Int winSize = state.winEnd - state.winBeg;
    const Int winSizeAfterFirst = winSize - state.firstBlockSize;
    if( winSizeAfterFirst == 0 )
    {
        LogicError("Single-block windows are not allowed.");
    }
    else
    {
        state.numWinBlocks = 2 + (winSizeAfterFirst-1)/state.blockSize;
        state.lastBlockSize = winSizeAfterFirst % state.blockSize;
    }

    // Apply the black-box function for determining the number of bulges
    state.numBulgesPerBlock = ctrl.numBulgesPerBlock(state.blockSize);

    // Store the local indices for bounding the off-diagonal application of the
    // accumulated Householder transformations for each diagonal block
    state.localTransformRowBeg =
      ( ctrl.fullTriangle ?
        0 :
        H.LocalRowOffset(state.winBeg) );
    state.localTransformColEnd =
      ( ctrl.fullTriangle ?
        H.LocalColOffset(n) :
        H.LocalColOffset(state.winEnd) );

    // Do some preliminary sanity checks
    if( winSize < 2*state.blockSize )
        LogicError("The window size must be at least twice the block size");
    if( state.blockSize != H.BlockWidth() )
        LogicError("IntraBlockChase assumes square distribution blocks"); 
    if( H.ColCut() != H.RowCut() )
        LogicError("IntraBlockChase assumes that the cuts are equal");
    if( state.numBulgesPerBlock*6 + 1 > state.blockSize )
        LogicError
        ("Cannot incorporate ",state.numBulgesPerBlock,
         " bulges into a distributed block size of ",state.blockSize);

    // The following will change as the active window is updated
    // ---------------------------------------------------------

    // Bulges are introduced starting from the bottom, and we introduce the
    // leftover number of shifts first.
    const Int bulgeLeftover = Mod(numBulges,state.numBulgesPerBlock);
    state.numBulgesInLastBlock =
      ( bulgeLeftover==0 ? state.numBulgesPerBlock : bulgeLeftover );
    state.bulgeBeg = numBulges - state.numBulgesInLastBlock;
    state.bulgeEnd = numBulges;

    // We initialize before any shifts have been introduced. If the first block
    // is of full size, then the initial inter-block chase will simply involve
    // introducing a packet into block 0; otherwise, a packet will be 
    // immediately chased through to block 1.
    const bool fullFirstBlock = ( state.blockSize == state.firstBlockSize );
    state.activeBlockBeg = 0;
    state.activeBeg = state.winBeg;
    if( fullFirstBlock )
    {
        state.activeBlockEnd = 1;
        state.activeEnd = state.winBeg + state.blockSize;
    }
    else
    {
        state.activeBlockEnd = 2;
        state.activeEnd = state.winBeg + state.firstBlockSize + state.blockSize;
    }

    // We now compute the first row/column blocks of the active window that are
    // assigned to this process
    const int activeColAlign = H.RowOwner( state.activeBeg );
    const int activeRowAlign = H.ColOwner( state.activeBeg );
    const int activeColShift =
      Shift( grid.Row(), activeColAlign, grid.Height() );
    const int activeRowShift =
      Shift( grid.Col(), activeRowAlign, grid.Width() );
    state.activeRowBlockBeg = state.activeBlockBeg + activeColShift;
    state.activeColBlockBeg = state.activeBlockBeg + activeRowShift;

    return state;
}

// Update the chase state after performing an interblock chase followed by an
// intrablock chase.
template<typename F>
void AdvanceChaseState
( const DistMatrix<F,MC,MR,BLOCK>& H,
        DistChaseState& state )
{
    DEBUG_CSE
    const Grid& grid = H.Grid();

    if( state.activeBlockEnd == state.numWinBlocks )
    {
        // The last packet was chased off of the matrix
        state.bulgeEnd -= state.numBulgesInLastBlock;
        state.numBulgesInLastBlock = state.numBulgesPerBlock;
    }
    else
    {
        // The last packet was chased down a block of the matrix
        if( state.activeBlockEnd == state.numWinBlocks-1 )
            state.activeEnd += state.lastBlockSize;
        else
            state.activeEnd += state.blockSize;
        ++state.activeBlockEnd;
    }

    if( state.bulgeBeg == 0 )
    {
        // All of the packets are in play, so the sweep is winding down

        // Push the active block indices up one block
        if( state.activeBlockBeg == 0 &&
            state.firstBlockSize != state.blockSize )
        {
            state.activeBeg += state.firstBlockSize + state.blockSize;
            state.activeBlockBeg += 2; 
        }
        else
        {
            state.activeBeg += state.blockSize;
            state.activeBlockBeg += 1;
        }

        // Compute the active region alignments and first assigned rows/col's
        const int activeColAlign = H.RowOwner( state.activeBeg );
        const int activeRowAlign = H.ColOwner( state.activeBeg );
        const int activeColShift =
          Shift( grid.Row(), activeColAlign, grid.Height() );
        const int activeRowShift =
          Shift( grid.Col(), activeRowAlign, grid.Width() );
        state.activeRowBlockBeg = state.activeBlockBeg + activeColShift;
        state.activeColBlockBeg = state.activeBlockBeg + activeRowShift;
    }
    else
    {
        // Introduce another packet at the beginning of the window
        DEBUG_ONLY(
          if( state.activeBlockBeg != 0 )
              LogicError("Inconsistent chase state");
        )
        state.bulgeBeg -= state.numBulgesPerBlock;
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#include "./Dist/IntraBlockChase.hpp"
#include "./Dist/InterBlockChase.hpp"

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_DIST_HPP
