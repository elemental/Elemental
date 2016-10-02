/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_INTER_BLOCK_CHASE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_INTER_BLOCK_CHASE_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

// For the vast majority of the inter-block packet chases, the transformations
// should start in the bottom-right of one diagonal block and end in the
// top-left of the subsequent diagonal block. And, to support non-overlapping
// inter-block chases of "even" and "odd" diagonal blocks (with the parity
// determined by number of blocks of separation from the last block), we require
// that the distribution block size be at least one plus six times the number
// of bulges.
//
// If we had two bulges per block and a distribution block size of 16, then 
// most of the 'even' intra-block chases would begin in the form
//
//                          ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//     -------------------------------------------------------------------
//    | x x x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    | x x x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |   x x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |     x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |       x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |         x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |           x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |             x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |               x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |                 x E E E E x x x | x x x x x x x x x x x x x x x x |
//  ~ |                   E E E E x x x | x x x x x x x x x x x x x x x x |
//  ~ |                   E E E E x x x | x x x x x x x x x x x x x x x x |
//  ~ |                   E E E E E E E | x x x x x x x x x x x x x x x x |
//  ~ |                         E E E E | x x x x x x x x x x x x x x x x |
//  ~ |                         E E E E | x x x x x x x x x x x x x x x x |
//  ~ |                         E E E E | x x x x x x x x x x x x x x x x |
//    |---------------------------------|---------------------------------|,
//  ~ |                               x | x x x x x x x x x x x x x x x x |
//  ~ |                                 | x x x x x x x x x x x x x x x x |
//  ~ |                                 |   x x x x x x x x x x x x x x x |
//  ~ |                                 |     x x x x x x x x x x x x x x |
//  ~ |                                 |       x x x x x x x x x x x x x |
//  ~ |                                 |         x x x x x x x x x x x x |
//    |                                 |           x x x x x x x x x x x |
//    |                                 |             x x x x x x x x x x |
//    |                                 |               x x x x x x x x x |
//    |                                 |                 x O O O O x x x |
//    |                                 |                   O O O O x x x |
//    |                                 |                   O O O O x x x |
//    |                                 |                   O O O O O O O |
//    |                                 |                         O O O O |
//    |                                 |                         O O O O |
//    |                                 |                         O O O O |
//     -------------------------------------------------------------------
//
// with the even and odd bulges respectively marked with an 'E' and an 'O' and
// with the partitioning denoting process boundaries. The even chase would then
// end with the bulges in the positions shown by
//
//                          ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//     -------------------------------------------------------------------
//    | x x x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    | x x x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |   x x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |     x x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |       x x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |         x x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |           x x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |             x x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |               x x x x x x x x x | x x x x x x x x x x x x x x x x |
//    |                 x x x x x x x x | x x x x x x x x x x x x x x x x |
//  ~ |                   x x x x x x x | x x x x x x x x x x x x x x x x |
//  ~ |                     x x x x x x | x x x x x x x x x x x x x x x x |
//  ~ |                       x x x x x | x x x x x x x x x x x x x x x x |
//  ~ |                         x x x x | x x x x x x x x x x x x x x x x |
//  ~ |                           x x x | x x x x x x x x x x x x x x x x |
//  ~ |                             x x | x x x x x x x x x x x x x x x x |
//    |---------------------------------|---------------------------------|.
//  ~ |                               x | E E E E x x x x x x x x x x x x |
//  ~ |                                 | E E E E x x x x x x x x x x x x |
//  ~ |                                 | E E E E x x x x x x x x x x x x |
//  ~ |                                 | E E E E E E E x x x x x x x x x |
//  ~ |                                 |       E E E E x x x x x x x x x |
//  ~ |                                 |       E E E E x x x x x x x x x |
//    |                                 |       E E E E x x x x x x x x x |
//    |                                 |             x x x x x x x x x x |
//    |                                 |               x x x x x x x x x |
//    |                                 |                 x O O O O x x x |
//    |                                 |                   O O O O x x x |
//    |                                 |                   O O O O x x x |
//    |                                 |                   O O O O O O O |
//    |                                 |                         O O O O |
//    |                                 |                         O O O O |
//    |                                 |                         O O O O |
//     -------------------------------------------------------------------
//
// The smallest admissible block size (if we demand that the bulges are entirely
// within the distribution block; Cf. ScaLAPACK's minimum block size of 6) with
// two bulges per block would be 1 + 6*2 = 13, and the equivalent 'even' chase
// would begin with the bulges in the following positions:
//
//                    ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//     -------------------------------------------------------
//    | x x x x x x x x x x x x x | x x x x x x x x x x x x x |
//    | x x x x x x x x x x x x x | x x x x x x x x x x x x x |
//    |   x x x x x x x x x x x x | x x x x x x x x x x x x x |
//    |     x x x x x x x x x x x | x x x x x x x x x x x x x |
//    |       x x x x x x x x x x | x x x x x x x x x x x x x |
//    |         x x x x x x x x x | x x x x x x x x x x x x x |
//    |           x E E E E x x x | x x x x x x x x x x x x x |
//  ~ |             E E E E x x x | x x x x x x x x x x x x x |
//  ~ |             E E E E x x x | x x x x x x x x x x x x x |
//  ~ |             E E E E E E E | x x x x x x x x x x x x x |
//  ~ |                   E E E E | x x x x x x x x x x x x x |
//  ~ |                   E E E E | x x x x x x x x x x x x x |
//  ~ |                   E E E E | x x x x x x x x x x x x x |
//    |---------------------------|---------------------------|
//  ~ |                         x | x x x x x x x x x x x x x |
//  ~ |                           | x x x x x x x x x x x x x |
//  ~ |                           |   x x x x x x x x x x x x |
//  ~ |                           |     x x x x x x x x x x x |
//  ~ |                           |       x x x x x x x x x x |
//  ~ |                           |         x x x x x x x x x |
//    |                           |           x O O O O x x x |
//    |                           |             O O O O x x x |
//    |                           |             O O O O x x x |
//    |                           |             O O O O O O O |
//    |                           |                   O O O O |
//    |                           |                   O O O O |
//    |                           |                   O O O O |
//     -------------------------------------------------------
//
// and end with the bulges in the following positions:
//
//                    ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//     -------------------------------------------------------
//    | x x x x x x x x x x x x x | x x x x x x x x x x x x x |
//    | x x x x x x x x x x x x x | x x x x x x x x x x x x x |
//    |   x x x x x x x x x x x x | x x x x x x x x x x x x x |
//    |     x x x x x x x x x x x | x x x x x x x x x x x x x |
//    |       x x x x x x x x x x | x x x x x x x x x x x x x |
//    |         x x x x x x x x x | x x x x x x x x x x x x x |
//    |           x x x x x x x x | x x x x x x x x x x x x x |
//  ~ |             x x x x x x x | x x x x x x x x x x x x x |
//  ~ |               x x x x x x | x x x x x x x x x x x x x |
//  ~ |                 x x x x x | x x x x x x x x x x x x x |
//  ~ |                   x x x x | x x x x x x x x x x x x x |
//  ~ |                     x x x | x x x x x x x x x x x x x |
//  ~ |                       x x | x x x x x x x x x x x x x |
//    |---------------------------|---------------------------|.
//  ~ |                         x | E E E E x x x x x x x x x |
//  ~ |                           | E E E E x x x x x x x x x |
//  ~ |                           | E E E E x x x x x x x x x |
//  ~ |                           | E E E E E E E x x x x x x |
//  ~ |                           |       E E E E x x x x x x |
//  ~ |                           |       E E E E x x x x x x |
//    |                           |       E E E O O O O x x x |
//    |                           |             O O O O x x x |
//    |                           |             O O O O x x x |
//    |                           |             O O O O O O O |
//    |                           |                   O O O O |
//    |                           |                   O O O O |
//    |                           |                   O O O O |
//     -------------------------------------------------------
//
// We note that the overlap of the second even bulge with the first odd bulge
// is acceptable for the same reason that the two even bulges are allowed to
// overlap: the Householder transformations that pushed the second even bulge
// down the diagonal to overlap with the first odd bulge did not disturb the 
// odd bulge.
//
// Keeping with the minimal block size of 13, the subsequent 'odd' chase would
// typically begin in the form
//
//                    ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//     -------------------------------------------------------
//    | E E E E x x x x x x x x x | x x x x x x x x x x x x x |
//    | E E E E x x x x x x x x x | x x x x x x x x x x x x x |
//    | E E E E x x x x x x x x x | x x x x x x x x x x x x x |
//    | E E E E E E E x x x x x x | x x x x x x x x x x x x x |
//    |       E E E E x x x x x x | x x x x x x x x x x x x x |
//    |       E E E E x x x x x x | x x x x x x x x x x x x x |
//    |       E E E O O O O x x x | x x x x x x x x x x x x x |
//  ~ |             O O O O x x x | x x x x x x x x x x x x x |
//  ~ |             O O O O x x x | x x x x x x x x x x x x x |
//  ~ |             O O O O O O O | x x x x x x x x x x x x x |
//  ~ |                   O O O O | x x x x x x x x x x x x x |
//  ~ |                   O O O O | x x x x x x x x x x x x x |
//  ~ |                   O O O O | x x x x x x x x x x x x x |
//    |---------------------------|---------------------------|,
//  ~ |                         x | x x x x x x x x x x x x x |
//  ~ |                           | x x x x x x x x x x x x x |
//  ~ |                           |   x x x x x x x x x x x x |
//  ~ |                           |     x x x x x x x x x x x |
//  ~ |                           |       x x x x x x x x x x |
//  ~ |                           |         x x x x x x x x x |
//    |                           |           x x x x x x x x |
//    |                           |             x x x x x x x |
//    |                           |               x x x x x x |
//    |                           |                 x x x x x |
//    |                           |                   x x x x |
//    |                           |                     x x x |
//    |                           |                       x x |
//     -------------------------------------------------------
//
// and end in the form
//                    ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~ 
//     -------------------------------------------------------
//    | E E E E x x x x x x x x x | x x x x x x x x x x x x x |
//    | E E E E x x x x x x x x x | x x x x x x x x x x x x x |
//    | E E E E x x x x x x x x x | x x x x x x x x x x x x x |
//    | E E E E E E E x x x x x x | x x x x x x x x x x x x x |
//    |       E E E E x x x x x x | x x x x x x x x x x x x x |
//    |       E E E E x x x x x x | x x x x x x x x x x x x x |
//    |       E E E E x x x x x x | x x x x x x x x x x x x x |
//  ~ |             x x x x x x x | x x x x x x x x x x x x x |
//  ~ |               x x x x x x | x x x x x x x x x x x x x |
//  ~ |                 x x x x x | x x x x x x x x x x x x x |
//  ~ |                   x x x x | x x x x x x x x x x x x x |
//  ~ |                     x x x | x x x x x x x x x x x x x |
//  ~ |                       x x | x x x x x x x x x x x x x |
//    |---------------------------|---------------------------|.
//  ~ |                         x | O O O O x x x x x x x x x |
//  ~ |                           | O O O O x x x x x x x x x |
//  ~ |                           | O O O O x x x x x x x x x |
//  ~ |                           | O O O O O O O x x x x x x |
//  ~ |                           |       O O O O x x x x x x |
//  ~ |                           |       O O O O x x x x x x |
//    |                           |       O O O O x x x x x x |
//    |                           |             x x x x x x x |
//    |                           |               x x x x x x |
//    |                           |                 x x x x x |
//    |                           |                   x x x x |
//    |                           |                     x x x |
//    |                           |                       x x |
//     -------------------------------------------------------
//
// As a consequence of the windowing scheme, it is possible for the upper-left
// corner of H11 to begin in the middle of a block of the distributed matrix.
// In these cases, we modify the inter-block chase that would have introduced
// the bulges into said partial block to directly put them into the upper-left
// corner of the *next* diagonal block (and call it an "introductory"
// inter-block chase). For example, if the distribution block size was 13, but
// the window begins in the tenth position of a distribution block, then the
// upper-left inter-block chase of two bulges would begin in the
// form
//
//        ~ ~ ~   ~ ~ ~ ~ ~ ~ 
//       -----------------------------------
//   ~  | x x x | x x x x x x x x x x x x x |
//   ~  | x x x | x x x x x x x x x x x x x |
//   ~  |   x x | x x x x x x x x x x x x x |
//      |-------|---------------------------|
//   ~  |     x | x x x x x x x x x x x x x |
//   ~  |       | x x x x x x x x x x x x x |
//   ~  |       |   x x x x x x x x x x x x |
//   ~  |       |     x x x x x x x x x x x |
//   ~  |       |       x x x x x x x x x x |
//   ~  |       |         x x x x x x x x x |
//      |       |           x x x x x x x x |
//      |       |             x x x x x x x |
//      |       |               x x x x x x |
//      |       |                 x x x x x |
//      |       |                   x x x x |
//      |       |                     x x x |
//      |       |                       x x |
//       -----------------------------------
//
// and end in the form
//
//        ~ ~ ~   ~ ~ ~ ~ ~ ~ 
//       -----------------------------------
//   ~  | x x x | x x x x x x x x x x x x x |
//   ~  | x x x | x x x x x x x x x x x x x |
//   ~  |   x x | x x x x x x x x x x x x x |
//      |-------|---------------------------|
//   ~  |     x | B B B B x x x x x x x x x |
//   ~  |       | B B B B x x x x x x x x x |
//   ~  |       | B B B B x x x x x x x x x |
//   ~  |       | B B B B B B B x x x x x x |
//   ~  |       |       B B B B x x x x x x |
//   ~  |       |       B B B B x x x x x x |,
//      |       |       B B B B x x x x x x |
//      |       |             x x x x x x x |
//      |       |               x x x x x x |
//      |       |                 x x x x x |
//      |       |                   x x x x |
//      |       |                     x x x |
//      |       |                       x x |
//       -----------------------------------
//
// usually with two bulges of odd parity packed into the bottom-right corner
// throughout the process. In cases where the upper-left corner of the window
// begins at the beginning of a distribution block, the above would simplify to
//
//       ~ ~ ~ ~ ~ ~                           ~ ~ ~ ~ ~ ~
//      ---------------------------           ---------------------------
//   ~ | x x x x x x x x x x x x x |       ~ | B B B B x x x x x x x x x |
//   ~ | x x x x x x x x x x x x x |       ~ | B B B B x x x x x x x x x |
//   ~ |   x x x x x x x x x x x x |       ~ | B B B B x x x x x x x x x |
//   ~ |     x x x x x x x x x x x |       ~ | B B B B B B B x x x x x x |
//   ~ |       x x x x x x x x x x |  |->  ~ |       B B B B x x x x x x |
//   ~ |         x x x x x x x x x |       ~ |       B B B B x x x x x x |
//     |           x x x x x x x x |         |       B B B B x x x x x x |,
//     |             x x x x x x x |         |             x x x x x x x |
//     |               x x x x x x |         |               x x x x x x |
//     |                 x x x x x |         |                 x x x x x |
//     |                   x x x x |         |                   x x x x |
//     |                     x x x |         |                     x x x |
//     |                       x x |         |                       x x |
//      ---------------------------           ---------------------------
//
// again, usually with two bulges of opposite parity packed into the
// bottom-right corner. We emphasize that, in this simple case where only one
// diagonal block is involved in the introductory inter-block chase, the parity
// of the single diagonal block should be opposite to the parity of the chase,
// as an even first block should have bulges created during the odd chase due to
// the interpretation of an implicit zero-by-zero even parity block living in
// the top-left corner. The reason for this interpretation is that it cleanly
// generalizes to the partial block case.
//
// Lastly, if we chase a packet of bulges into the last diagonal block of the
// window, regardless of whether or not it is a full diagonal block, we
// immediately chase them out of the window. For example, we might begin in the
// form
//
//                      ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//       -----------------------------------------
//      | x x x x x x x x x x x x x | x x x x x x |
//      | x x x x x x x x x x x x x | x x x x x x |
//      |   x x x x x x x x x x x x | x x x x x x |
//      |     x x x x x x x x x x x | x x x x x x |
//      |       x x x x x x x x x x | x x x x x x |
//      |         x x x x x x x x x | x x x x x x |
//      |           x O O O O x x x | x x x x x x |
//    ~ |             O O O O x x x | x x x x x x |
//    ~ |             O O O O x x x | x x x x x x |
//    ~ |             O O O O O O O | x x x x x x |,
//    ~ |                   O O O O | x x x x x x |
//    ~ |                   O O O O | x x x x x x |
//    ~ |                   O O O O | x x x x x x |
//      |---------------------------|-------------|
//    ~ |                         x | x x x x x x |
//    ~ |                           | x x x x x x |
//    ~ |                           |   x x x x x |
//    ~ |                           |     x x x x |
//    ~ |                           |       x x x |
//    ~ |                           |         x x |
//       -----------------------------------------
//
// and end in the form
//
//                      ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//       -----------------------------------------
//      | x x x x x x x x x x x x x | x x x x x x |
//      | x x x x x x x x x x x x x | x x x x x x |
//      |   x x x x x x x x x x x x | x x x x x x |
//      |     x x x x x x x x x x x | x x x x x x |
//      |       x x x x x x x x x x | x x x x x x |
//      |         x x x x x x x x x | x x x x x x |
//      |           x x x x x x x x | x x x x x x |
//    ~ |             x x x x x x x | x x x x x x |
//    ~ |               x x x x x x | x x x x x x |
//    ~ |                 x x x x x | x x x x x x |.
//    ~ |                   x x x x | x x x x x x |
//    ~ |                     x x x | x x x x x x |
//    ~ |                       x x | x x x x x x |
//      |---------------------------|-------------|
//    ~ |                         x | x x x x x x |
//    ~ |                           | x x x x x x |
//    ~ |                           |   x x x x x |
//    ~ |                           |     x x x x |
//    ~ |                           |       x x x |
//    ~ |                           |         x x |
//       -----------------------------------------
//
// As in the case of the "introductory" inter-block chase, we perform such an
// "exit" inter-block chase during the stage implied by the parity of the 
// top-left block.
//
// In order to prevent an inter-block chase that is both an "introduction" and
// an "exit", we require that the distributed algorithm is only run on windows
// of size at least 2 * blockHeight.

enum DistChaseType {
  SIMPLE_INTRO_CHASE,
  COUPLED_INTRO_CHASE,
  STANDARD_CHASE,
  EXIT_CHASE,
  NO_CHASE
};

struct InterBlockInteraction
{
  DistChaseType chaseType;
  Int block0;
  // We always have block1 = block0 + 1 for non-trivial interactions, but we
  // keep it here for the sake of simplicity.
  Int block1;
  Int blockSize0;
  Int blockSize1;
  Int numBulges;

  // The (entry-wise) indices of the beginning and end of the interaction
  Int beg;
  Int end;
};

namespace interblock {

inline InterBlockInteraction
DetermineInteraction
( bool evenToOdd,
  Int diagBlockRow,
  const Grid& grid,
  const DistChaseState& state,
  const DistChaseContext& context )
{
    DEBUG_CSE
    if( diagBlockRow < state.activeBlockBeg ||
        diagBlockRow >= state.activeBlockEnd )
        LogicError("Diagonal block row was not in the active range");

    const bool fullFirstBlock = ( context.firstBlockSize == context.blockSize );

    const Int distFromEnd = state.activeBlockEnd - diagBlockRow;
    const bool evenBlock = ( Mod( distFromEnd, 2 ) == 0 );
    const bool sameParity = ( evenBlock == evenToOdd );
    if( diagBlockRow+1 == state.activeBlockEnd && sameParity )
        LogicError("This step would push a packet out of the active window");

    const int gridCol = grid.Col();
    const int ownerCol = Mod( context.winRowAlign+diagBlockRow, grid.Width() );
    const int prevCol = Mod( context.winRowAlign-1+diagBlockRow, grid.Width() );
    const int nextCol = Mod( context.winRowAlign+1+diagBlockRow, grid.Width() );

    InterBlockInteraction interaction;
    interaction.chaseType = NO_CHASE;
    interaction.block0 = -1;
    interaction.block1 = -1;
    interaction.blockSize0 = 0;
    interaction.blockSize1 = 0;
    if( diagBlockRow == 0 )
    {
        if( sameParity )
        {
            // We are pushing the packet (possibly after introducing it)
            if( gridCol == ownerCol || gridCol == nextCol )
            { 
                if( fullFirstBlock )
                {
                    // A packet is already lying the (0,0) diagonal block.
                    interaction.chaseType = STANDARD_CHASE;
                    interaction.block0 = 0;
                    interaction.block1 = 1;
                    interaction.blockSize0 = context.blockSize;
                    interaction.blockSize1 = context.blockSize;
                }
                else
                {
                    // We need to introduce and chase a packet.
                    interaction.chaseType = COUPLED_INTRO_CHASE;
                    interaction.block0 = 0;
                    interaction.block1 = 1;
                    interaction.blockSize0 = context.firstBlockSize;
                    interaction.blockSize1 = context.blockSize; 
                }
            }
        }
        else
        {
            if( gridCol == ownerCol )
            {
                if( fullFirstBlock )
                {
                    // We will introduce a packet into this block
                    interaction.chaseType = SIMPLE_INTRO_CHASE;
                    interaction.block0 = -1;
                    interaction.block1 = 0;
                    interaction.blockSize0 = 0;
                    interaction.blockSize1 = context.blockSize;
                }
                else
                {
                    // No packets are left in partial blocks; so our process
                    // does not interact via this row.
                }
            }
        }
    }
    else if( diagBlockRow == 1 )
    {
        if( sameParity )
        {
            // We are pushing the packet
            if( gridCol == ownerCol || gridCol == nextCol )
            {
                // A packet was put in the (1,1) diagonal block in the previous
                // chase (of opposite parity); this chase will push it to the
                // (2,2) block.
                if( context.numWinBlocks == 3 )
                {
                    interaction.chaseType = EXIT_CHASE;
                    interaction.block0 = 1;
                    interaction.block1 = 2;
                    interaction.blockSize0 = context.blockSize;
                    interaction.blockSize1 = context.lastBlockSize;
                }
                else
                {
                    interaction.chaseType = STANDARD_CHASE;
                    interaction.block0 = 1;
                    interaction.block1 = 2;
                    interaction.blockSize0 = context.blockSize;
                    interaction.blockSize1 = context.blockSize;
                }
            }
        }
        else
        {
            // We are pulling the packet
            if( gridCol == ownerCol || gridCol == prevCol )
            {
                if( fullFirstBlock )
                {
                    // This chase will pull a packet from (0,0) to (1,1)
                    interaction.chaseType = STANDARD_CHASE;
                    interaction.block0 = 0;
                    interaction.block1 = 1;
                    interaction.blockSize0 = context.blockSize;
                    interaction.blockSize1 = context.blockSize;
                }
                else
                {
                    // This chase will introduce a packet and immediately chase 
                    // it into diagonal block (1,1)
                    interaction.chaseType = COUPLED_INTRO_CHASE;
                    interaction.block0 = 0;
                    interaction.block1 = 1;
                    interaction.blockSize0 = context.firstBlockSize;
                    interaction.blockSize1 = context.blockSize;
                }
            }
        }
    }
    else if( diagBlockRow == context.numWinBlocks-2 )
    {
        if( sameParity )
        {
            // We will push a packet into the last block (and then exit it)
            if( gridCol == ownerCol || gridCol == nextCol )
            {
                interaction.chaseType = EXIT_CHASE;
                interaction.block0 = context.numWinBlocks-2;
                interaction.block1 = context.numWinBlocks-1;
                interaction.blockSize0 = context.blockSize;
                interaction.blockSize1 = context.lastBlockSize; 
            }
        }
        else
        {
            // We will pull a packet into the next-to-last diagonal block
            if( gridCol == ownerCol || gridCol == prevCol )
            {
                interaction.chaseType = STANDARD_CHASE;
                interaction.block0 = context.numWinBlocks-3;
                interaction.block1 = context.numWinBlocks-2;
                interaction.blockSize0 = context.blockSize;
                interaction.blockSize1 = context.blockSize;
            }
        }
    }
    else if( diagBlockRow == context.numWinBlocks-1 )
    {
        if( sameParity )
        {
            // Packets are never left in the last diagonal block
        }
        else
        {
            // We will pull a packet into this last diagonal block and exit it
            if( gridCol == ownerCol || gridCol == prevCol )
            {
                interaction.chaseType = EXIT_CHASE;
                interaction.block0 = context.numWinBlocks-2;
                interaction.block1 = context.numWinBlocks-1;
                interaction.blockSize0 = context.blockSize;
                interaction.blockSize1 = context.lastBlockSize;
            }
        }
    }
    else
    {
        if( sameParity )
        {
            // We will push a packet into the next diagonal block
            if( gridCol == ownerCol || gridCol == nextCol )
            {
                interaction.chaseType = STANDARD_CHASE;
                interaction.block0 = diagBlockRow;
                interaction.block1 = diagBlockRow+1;
                interaction.blockSize0 = context.blockSize;
                interaction.blockSize1 = context.blockSize;
            }
        }
        else
        {
            // We will pull a packet into this diagonal block.
            if( gridCol == ownerCol || gridCol == prevCol )
            {
                interaction.chaseType = STANDARD_CHASE;
                interaction.block0 = diagBlockRow-1;
                interaction.block1 = diagBlockRow;
                interaction.blockSize0 = context.blockSize;
                interaction.blockSize1 = context.blockSize;
            }
        }
    }

    // The number of bulges is guaranteed to be equal to
    // context.numBulgesPerBlock except (possibly) in the last active
    // interaction, which must have its second block at position
    // state.activeBlockEnd-1.
    interaction.numBulges =
      ( interaction.block1 < state.activeBlockEnd-1 ?
        context.numBulgesPerBlock :
        context.numBulgesInLastBlock );

    interaction.beg = context.winBeg +
      ( interaction.block0 <= 0 ?
        0 :
        context.firstBlockSize + (interaction.block0-1)*context.blockSize );
    interaction.end = interaction.beg + interaction.blockSize0 +
      interaction.blockSize1;

    return interaction;
}

inline vector<InterBlockInteraction>
FormInteractionList
( bool evenToOdd,
  const Grid& grid,
  const DistChaseState& state,
  const DistChaseContext& context )
{
    vector<InterBlockInteraction> interactionList;
    // Only loop over the row blocks that are assigned to our process row
    // and occur within the active window.
    Int diagBlock = context.activeRowBlockBeg;
    while( diagBlock < state.activeBlockEnd )
    {
        auto interaction =
          interblock::DetermineInteraction
          ( evenToOdd, diagBlock, grid, state, context );
        if( interaction.chaseType == NO_CHASE )
        {
            diagBlock += grid.Height();
        }
        else
        {
            interactionList.push_back(interaction);
            // We must take care to not participate twice in one chase
            if( grid.Height() == 1 )
                diagBlock = interaction.block1 + 1;
            else
                diagBlock += grid.Height();
        }
    }
    return interactionList;
}

// TODO(poulson): Shave down the communication volume by only sending the
// necessary pieces of the interaction window. For example, if one bulge was
// being chased from the end of a block of size 7 into another block of size 7,
// we would have the transition
//
//              ~ ~ ~   ~ ~ ~                           ~ ~ ~   ~ ~ ~
//     -------------------------------         -------------------------------
//    | x x x x x x x | x x x x x x x |       | x x x x x x x | x x x x x x x |
//    | x x x x x x x | x x x x x x x |       | x x x x x x x | x x x x x x x |
//    |   x x x x x x | x x x x x x x |       |   x x x x x x | x x x x x x x |
//    |     x B B B B | x x x x x x x |       |     x x x x x | x x x x x x x |
//  ~ |       B B B B | x x x x x x x |     ~ |       x x x x | x x x x x x x |
//  ~ |       B B B B | x x x x x x x |     ~ |         x x x | x x x x x x x |
//  ~ |       B B B B | x x x x x x x |     ~ |           x x | x x x x x x x |
//    |---------------|---------------| |->   |---------------|---------------|,
//  ~ |             x | x x x x x x x |     ~ |             x | B B B B x x x |
//  ~ |               | x x x x x x x |     ~ |               | B B B B x x x |
//  ~ |               |   x x x x x x |     ~ |               | B B B B x x x |
//    |               |     x x x x x |       |               | B B B B x x x |
//    |               |       x x x x |       |               |       x x x x |
//    |               |         x x x |       |               |         x x x |
//    |               |           x x |       |               |           x x |
//     -------------------------------         -------------------------------
//
// of which only six rows and six columns are involved in the inter-block packet
// chase. We can also obviously only transmit a single entry of the bottom-left
// block (and the nonzero portions of the other quadrants).
//
template<typename F>
void CollectInterBlock
(       InterBlockInteraction interaction,
  const DistMatrix<F,MC,MR,BLOCK>& H,
  const Grid& grid,
  const DistChaseState& state,
  const DistChaseContext& context, 
        Matrix<F>& HBlock )
{
    DEBUG_CSE
    const Int blockSize = context.blockSize;
    const Int firstBlockSize = context.firstBlockSize;
    const auto& HLoc = H.LockedMatrix();

    if( interaction.chaseType == SIMPLE_INTRO_CHASE )
    {
        // Only a single process participates in introductory chases, and they
        // occur over the entire top-left block (which must have been full).
        const Int indexBeg = context.winBeg;
        const Int indexEnd = context.winBeg + blockSize;
        const Int localRowBeg = H.LocalRowOffset( indexBeg );
        const Int localRowEnd = H.LocalRowOffset( indexEnd );
        const Int localColBeg = H.LocalColOffset( indexBeg );
        const Int localColEnd = H.LocalColOffset( indexEnd );
        auto HInteractLoc =
          HLoc( IR(localRowBeg,localRowEnd), IR(localColBeg,localColEnd) );
        HBlock = HInteractLoc;
        return;
    }
    else if( interaction.chaseType == NO_CHASE )
        LogicError("Invalid request to collect an inter-block window");

    // Only "simple" introductory chases label their first block as "-1".
    const Int block0 = interaction.block0;

    // (Up to) four processes may participate in introducing and chasing
    // a packet into and out of a partial block.
    const int firstRow = Mod( context.winColAlign+block0, grid.Height() );
    const int secondRow = Mod( firstRow+1, grid.Height() );
    const int firstCol = Mod( context.winRowAlign+block0, grid.Width() );
    const int secondCol = Mod( firstCol+1, grid.Width() );

    // We can grab the indices of our local portion of the 2x2 interaction
    // window in a black-box manner.
    const Int localRowBeg = H.LocalRowOffset( interaction.beg );
    const Int localRowEnd = H.LocalRowOffset( interaction.end );
    const Int localColBeg = H.LocalColOffset( interaction.beg );
    const Int localColEnd = H.LocalColOffset( interaction.end );
    auto HInteractLoc =
      HLoc( IR(localRowBeg,localRowEnd), IR(localColBeg,localColEnd) );

    // The interior blocks are all full, and we know the first block size and
    // the inter-block interaction size, so we can easily compute the two
    // interaction block sizes.
    const Int interactionSize = interaction.blockSize0 + interaction.blockSize1;
    const auto ind0 = IR(0,interaction.blockSize0);
    const auto ind1 = IR(interaction.blockSize0,interactionSize);
    Zeros( HBlock, interactionSize, interactionSize );

    if( grid.Height() == 1 && grid.Width() == 1 )
    {
        // Only our process participates
        HBlock = HInteractLoc;
    }
    else if( grid.Height() == 1 )
    {
        // Two processes in the same row participate
        auto HBlockLeft = HBlock( ALL, ind0 );
        auto HBlockRight = HBlock( ALL, ind1 );
        if( grid.Col() == firstCol )
        {
            HBlockLeft = HInteractLoc;
            El::Send( HBlockLeft, secondCol, grid.RowComm() );
            El::Recv( HBlockRight, secondCol, grid.RowComm() );
        }
        else
        {
            HBlockRight = HInteractLoc;
            El::Recv( HBlockLeft, firstCol, grid.RowComm() );
            El::Send( HBlockRight, firstCol, grid.RowComm() );
        }
    }
    else if( grid.Width() == 1 )
    {
        // Two processes in the same column participate
        auto HBlockTop = HBlock( ind0, ALL );
        auto HBlockBottom = HBlock( ind1, ALL );
        if( grid.Row() == firstRow )
        {
            HBlockTop = HInteractLoc; 
            El::Send( HBlockTop, secondRow, grid.ColComm() );
            El::Recv( HBlockBottom, secondRow, grid.ColComm() );
        }
        else
        {
            HBlockBottom = HInteractLoc;
            El::Recv( HBlockTop, secondRow, grid.ColComm() );
            El::Send( HBlockBottom, secondRow, grid.ColComm() );
        }
    }
    else
    {
        // Four processes participate, though only the upper-left and
        // bottom-right ones will chase the packet, so only they receive
        // any data.
        auto HBlock00 = HBlock( ind0, ind0 );
        auto HBlock01 = HBlock( ind0, ind1 );
        auto HBlock10 = HBlock( ind1, ind0 );
        auto HBlock11 = HBlock( ind1, ind1 );
        // We will use the column-major ordering (which is the VC comm.)
        const int proc00 = firstRow + firstCol*grid.Height();
        const int proc01 = firstRow + secondCol*grid.Height();
        const int proc10 = secondRow + firstCol*grid.Height();
        const int proc11 = secondRow + secondCol*grid.Height();
        if( grid.Row() == firstRow && grid.Col() == firstCol )
        {
            HBlock00 = HInteractLoc;
            // Receive the off-diagonal blocks
            El::Recv( HBlock01, grid.VCComm(), proc01 ); 
            El::Recv( HBlock10, grid.VCComm(), proc10 );
            // Exchange diagonal blocks with proc11
            El::Recv( HBlock11, grid.VCComm(), proc11 );
            El::Send( HBlock00, grid.VCComm(), proc11 );
        }
        else if( grid.Row() == firstRow && grid.Col() == secondCol )
        {
            HBlock01 = HInteractLoc;
            // Send our off-diagonal block to the two diagonal processes
            El::Send( HBlock01, grid.VCComm(), proc00 );
            El::Send( HBlock01, grid.VCComm(), proc11 );
        }
        else if( grid.Row() == secondRow && grid.Col() == firstCol )
        {
            HBlock10 = HInteractLoc;
            // Send our off-diagonal block to the two diagonal processes
            El::Send( HBlock10, grid.VCComm(), proc00 );
            El::Send( HBlock10, grid.VCComm(), proc11 );
        }
        else if( grid.Row() == secondRow && grid.Col() == secondCol )
        {
            HBlock11 = HInteractLoc;
            // Receive the off-diagonal blocks
            El::Recv( HBlock01, grid.VCComm(), proc01 );
            El::Recv( HBlock10, grid.VCComm(), proc10 );
            // Exchange diagonal blocks with proc00
            El::Send( HBlock11, grid.VCComm(), proc00 );
            El::Recv( HBlock00, grid.VCComm(), proc00 );
        }
    }
}

} // namespace interblock

template<typename F>
void InterBlockChase
(       DistMatrix<F,MC,MR,BLOCK>& H,
        DistMatrix<F,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
        bool evenToOdd,
  const DistChaseState& state,
  const DistChaseContext& context,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Grid& grid = H.Grid();
    const bool fullFirstBlock = ( context.firstBlockSize == context.blockSize );
    
    if( state.activeBlockBeg < 0 )
        LogicError("state.activeBlockBeg was negative");
    if( state.activeBlockBeg > state.activeBlockEnd )
        LogicError("state.activeBlockBeg > state.activeBlockEnd");
    if( state.activeBlockEnd == 1 && !fullFirstBlock )
        LogicError("Cannot introduce any bulges");

    auto interactionList =
      interblock::FormInteractionList( evenToOdd, grid, state, context );

    const Int numInteractions = interactionList.size();
    vector<Matrix<F>> UList(numInteractions);

    // For now, the following is just meant to exercise CollectInterBlock
    for( Int whichInteraction=0; whichInteraction<numInteractions;
         ++whichInteraction )
    {
        Matrix<F> HBlock;
        auto interaction = interactionList[whichInteraction];
        CollectInterBlock( interaction, H, state, context, HBlock );
    }

    // TODO(poulson): Finish implementing this routine
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTER_BLOCK_CHASE_HPP
