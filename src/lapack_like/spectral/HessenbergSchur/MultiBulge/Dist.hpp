/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_DIST_HPP
#define EL_SCHUR_HESS_MULTIBULGE_DIST_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

struct DistChaseState
{
  Int activeBlockBeg; // The block index of the first active block
  Int activeBlockEnd; // The block index of the end of the active blocks

  // The activeBlock{Beg,End} indices are relative to the active window implied
  // by the coupled HessenbergSchurCtrl.win{Beg,End} parameters.

  Int shiftBeg; // The first active shift index
  Int shiftEnd; // The end of the active shift indices
};

// A structure for storing auxiliary metadata that is somewhat tedious to 
// repeatedly recompute
struct DistChaseContext
{
  // These are the same as ctrl.win{Beg,End} but with END replaced with 'n'
  Int winBeg, winEnd;
  // The process rows and columns owning the block at the top-left of the window
  Int winRowAlign, winColAlign;
  
  // The full block size
  Int blockSize;
  // The size of the piece of the bottom-right quadrant of the distribution
  // block that the window starts at. For example, if the distribution block
  // size is 16 but the window starts two diagonals down the diagonal of such 
  // a block, then 'firstBlockSize' is 14.
  Int firstBlockSize;

  // The number of (full or partial) diagonal blocks in the window.
  Int numWinBlocks; 

  // The index of the diagonal entry that the first active shift is in. If 
  // activeBlockBeg is zero, then this is the top-left entry of the bottom-right
  // quadrant of the distribution block lying within the window; otherwise, it
  // is the index of the top-left entry of the full distribution block.
  Int activeBeg;
  // The first diagonal block assigned to our process row that occurs on or
  // after activeBlockBeg
  Int activeRowBlockBeg;
  // The first diagonal block assigned to our process column that occurs on or
  // after activeBlockBeg
  Int activeColBlockBeg; 

  // The number of bulges in the vast majority of the packets
  Int numBulgesPerBlock;
  // The number of bulges currently in the very last active block
  Int numBulgesInLastBlock;

  // The first local row index of H that should be updated
  Int localTransformRowBeg;
  // The last local column index of H that should be updated
  Int localTransformColEnd;
};

template<typename F>
DistChaseContext BuildDistChaseContext
( const DistMatrix<F,MC,MR,BLOCK>& H,
  const DistMatrix<F,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int numShifts = shifts.Height();
    const Int numBulges = numShifts / 2;
    const Grid& grid = H.Grid();

    DistChaseContext context;
    context.winBeg = ( ctrl.winBeg == END ? n : ctrl.winBeg );
    context.winEnd = ( ctrl.winEnd == END ? n : ctrl.winEnd );
    context.winRowAlign = H.ColOwner( context.winBeg );
    context.winColAlign = H.RowOwner( context.winBeg );

    context.blockSize = H.BlockHeight();

    // Compute the offset into the distribution block that the window begins in
    const Int winCut = Mod( context.winBeg + H.ColCut(), context.blockSize );

    // Compute the remaining size of the first block in the window
    context.firstBlockSize = context.blockSize - winCut;

    // Compute the number of (full or partial) diagonal blocks in the window 
    const Int winSize = context.winEnd - context.winBeg;
    const Int winSizeAfterFirst = winSize - context.firstBlockSize;
    if( winSizeAfterFirst == 0 )
    {
        // This should never be allowed to happen, but we handle it anyway.
        context.numWinBlocks = 1;
    }
    else
    {
        context.numWinBlocks = 2 + (winSizeAfterFirst-1)/context.blockSize;
    }

    // Apply the black-box function for determining the number of bulges
    context.numBulgesPerBlock = ctrl.numBulgesPerBlock(context.blockSize);

    // We push the leftover shifts through first, so the last block only has
    // the leftover number of shifts if shiftEnd is equal to the number of
    // shifts.
    context.numBulgesInLastBlock =
      ( state.shiftEnd == numShifts ?
        Mod(numBulges,context.numBulgesPerBlock) :
        context.numBulgesPerBlock );

    // Compute the row and column owners of the beginning of the active bulge
    // region and the first block assigned to us.
    context.activeBeg = context.winBeg +
      ( state.activeBlockBeg == 0 ?
        0 :
        context.firstBlockSize + (state.activeBlockBeg-1)*context.blockSize );
    const int activeColAlign = H.RowOwner( context.activeBeg );
    const int activeRowAlign = H.ColOwner( context.activeBeg );
    const int activeColShift =
      Shift( grid.Row(), activeColAlign, grid.Height() );
    const int activeRowShift =
      Shift( grid.Col(), activeRowAlign, grid.Width() );
    context.activeRowBlockBeg = state.activeBlockBeg + activeColShift;
    context.activeColBlockBeg = state.activeBlockBeg + activeRowShift;

    // Store the local indices for bounding the off-diagonal application of the
    // accumulated Householder transformations for each diagonal block
    context.localTransformRowBeg =
      ( ctrl.fullTriangle ?
        0 :
        H.LocalRowOffset(context.winBeg) );
    context.localTransformColEnd =
      ( ctrl.fullTriangle ?
        H.LocalColOffset(n) :
        H.LocalColOffset(context.winEnd) );

    // Do some preliminary sanity checks
    if( winSize < 2*context.blockSize )
        LogicError("The window size must be at least twice the block size");
    if( context.blockSize != H.BlockWidth() )
        LogicError("IntraBlockChase assumes square distribution blocks"); 
    if( H.ColCut() != H.RowCut() )
        LogicError("IntraBlockChase assumes that the cuts are equal");
    if( ctrl.wantSchurVecs )
    {
        // Ensure that H and Z have the same distributions
        if( Z.DistData() != H.DistData() )
            LogicError("The distributions of H and Z should match");
    }
    if( context.numBulgesPerBlock*6 + 1 > context.blockSize )
        LogicError
        ("Cannot incorporate ",context.numBulgesPerBlock,
         " bulges into a distributed block size of ",context.blockSize);

    return context;
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#include "./Dist/IntraBlockChase.hpp"
#include "./Dist/InterBlockChase.hpp"

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_DIST_HPP
