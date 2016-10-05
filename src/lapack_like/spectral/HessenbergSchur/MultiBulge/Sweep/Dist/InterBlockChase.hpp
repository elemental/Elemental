/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTER_BLOCK_CHASE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTER_BLOCK_CHASE_HPP

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

  // The (entry-wise) indices of the beginning and end of the range effected
  // by the Householder transformations
  Int householderBeg;
  Int householderEnd;

  bool participating;
};

namespace interblock {

inline InterBlockInteraction
DetermineInteraction
( bool evenToOdd,
  Int diagBlockRow,
  const Grid& grid,
  const DistChaseState& state )
{
    DEBUG_CSE
    if( diagBlockRow < state.activeBlockBeg ||
        diagBlockRow >= state.activeBlockEnd )
        LogicError("Diagonal block row was not in the active range");

    const bool fullFirstBlock = ( state.firstBlockSize == state.blockSize );

    const Int distFromEnd = state.activeBlockEnd - diagBlockRow;
    const bool evenBlock = ( Mod( distFromEnd, 2 ) == 0 );
    const bool sameParity = ( evenBlock == evenToOdd );
    if( diagBlockRow+1 == state.activeBlockEnd && sameParity )
        LogicError("This step would push a packet out of the active window");

    InterBlockInteraction interaction;
    interaction.chaseType = NO_CHASE;
    interaction.block0 = ( sameParity ? diagBlockRow : diagBlockRow-1 );
    interaction.block1 = ( sameParity ? diagBlockRow+1 : diagBlockRow );
    interaction.blockSize0 = 0;
    interaction.blockSize1 = 0;
    if( diagBlockRow == 0 )
    {
        if( sameParity )
        {
            // This block pushes the packet (possibly after introducing it)
            if( fullFirstBlock )
            {
                // A packet is already lying the (0,0) diagonal block.
                interaction.chaseType = STANDARD_CHASE;
                interaction.blockSize0 = state.blockSize;
                interaction.blockSize1 = state.blockSize;
            }
            else
            {
                // We need to introduce and chase a packet.
                interaction.chaseType = COUPLED_INTRO_CHASE;
                interaction.blockSize0 = state.firstBlockSize;
                interaction.blockSize1 = state.blockSize; 
            }
        }
        else
        {
            if( fullFirstBlock )
            {
                // We will introduce a packet into this block
                interaction.chaseType = SIMPLE_INTRO_CHASE;
                interaction.blockSize0 = 0;
                interaction.blockSize1 = state.blockSize;
            }
            else
            {
                // No packets are left in partial blocks; so our process
                // does not interact via this row.
                interaction.chaseType = NO_CHASE;
            }
        }
    }
    else if( diagBlockRow == 1 )
    {
        if( sameParity )
        {
            // A packet was put in the (1,1) diagonal block in the previous
            // chase (of opposite parity); this chase will push it to the
            // (2,2) block.
            if( state.numWinBlocks == 3 )
            {
                interaction.chaseType = EXIT_CHASE;
                interaction.blockSize0 = state.blockSize;
                interaction.blockSize1 = state.lastBlockSize;
            }
            else
            {
                interaction.chaseType = STANDARD_CHASE;
                interaction.blockSize0 = state.blockSize;
                interaction.blockSize1 = state.blockSize;
            }
        }
        else
        {
            if( fullFirstBlock )
            {
                // This chase will pull a packet from (0,0) to (1,1)
                interaction.chaseType = STANDARD_CHASE;
                interaction.blockSize0 = state.blockSize;
                interaction.blockSize1 = state.blockSize;
            }
            else
            {
                // This chase will introduce a packet and immediately chase 
                // it into diagonal block (1,1)
                interaction.chaseType = COUPLED_INTRO_CHASE;
                interaction.blockSize0 = state.firstBlockSize;
                interaction.blockSize1 = state.blockSize;
            }
        }
    }
    else if( diagBlockRow == state.numWinBlocks-2 )
    {
        if( sameParity )
        {
            // We will push a packet into the last block (and then exit it)
            interaction.chaseType = EXIT_CHASE;
            interaction.blockSize0 = state.blockSize;
            interaction.blockSize1 = state.lastBlockSize; 
        }
        else
        {
            // We will pull a packet into the next-to-last diagonal block
            interaction.chaseType = STANDARD_CHASE;
            interaction.blockSize0 = state.blockSize;
            interaction.blockSize1 = state.blockSize;
        }
    }
    else if( diagBlockRow == state.numWinBlocks-1 )
    {
        if( sameParity )
        {
            // Packets are never left in the last diagonal block
            interaction.chaseType = NO_CHASE;
        }
        else
        {
            // We will pull a packet into this last diagonal block and exit it
            interaction.chaseType = EXIT_CHASE;
            interaction.blockSize0 = state.blockSize;
            interaction.blockSize1 = state.lastBlockSize;
        }
    }
    else
    {
        if( sameParity )
        {
            // We will push a packet into the next diagonal block
            interaction.chaseType = STANDARD_CHASE;
            interaction.blockSize0 = state.blockSize;
            interaction.blockSize1 = state.blockSize;
        }
        else
        {
            // We will pull a packet into this diagonal block.
            interaction.chaseType = STANDARD_CHASE;
            interaction.blockSize0 = state.blockSize;
            interaction.blockSize1 = state.blockSize;
        }
    }

    // The number of bulges is guaranteed to be equal to
    // state.numBulgesPerBlock except (possibly) in the last active
    // interaction, which must have its second block at position
    // state.activeBlockEnd-1.
    interaction.numBulges =
      ( interaction.block1 < state.activeBlockEnd-1 ?
        state.numBulgesPerBlock :
        state.numBulgesInLastBlock );

    interaction.beg = state.winBeg +
      ( interaction.block0 <= 0 ?
        0 :
        state.firstBlockSize + (interaction.block0-1)*state.blockSize );
    interaction.end = interaction.beg + interaction.blockSize0 +
      interaction.blockSize1;

    interaction.householderBeg =
      ( (interaction.chaseType==SIMPLE_INTRO_CHASE || 
         interaction.chaseType==COUPLED_INTRO_CHASE) ?
        interaction.beg :
        interaction.beg + interaction.blockSize0 - 3*interaction.numBulges );
    interaction.householderEnd =
      ( interaction.chaseType==EXIT_CHASE ?
        interaction.end :
        interaction.beg + interaction.blockSize0 + 3*interaction.numBulges );

    // (Up to) four processes may participate.
    const int firstRow =
      Mod( state.winColAlign+interaction.block0, grid.Height() );
    const int firstCol =
      Mod( state.winRowAlign+interaction.block0, grid.Width() );
    const int secondRow = Mod( firstRow+1, grid.Height() );
    const int secondCol = Mod( firstCol+1, grid.Width() );
    if( (grid.Row() == firstRow || grid.Row() == secondRow) &&
        (grid.Col() == firstCol || grid.Col() == secondCol) )
        interaction.participating = true;
    else
        interaction.participating = false;

    return interaction;
}

inline vector<InterBlockInteraction>
FormRowInteractionList
( bool evenToOdd,
  const Grid& grid,
  const DistChaseState& state )
{
    DEBUG_CSE
    vector<InterBlockInteraction> interactionList;
    // Only loop over the blocks that are assigned to our process row
    // and occur within the active window.
    Int diagBlock = state.activeRowBlockBeg;
    while( diagBlock < state.activeBlockEnd )
    {
        auto interaction =
          interblock::DetermineInteraction( evenToOdd, diagBlock, grid, state );
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

inline vector<InterBlockInteraction>
FormColumnInteractionList
( bool evenToOdd,
  const Grid& grid,
  const DistChaseState& state )
{
    DEBUG_CSE
    vector<InterBlockInteraction> interactionList;
    // Only loop over the blocks that are assigned to our process column
    // and occur within the active window.
    Int diagBlock = state.activeColBlockBeg;
    while( diagBlock < state.activeBlockEnd )
    {
        auto interaction =
          interblock::DetermineInteraction( evenToOdd, diagBlock, grid, state );
        if( interaction.chaseType == NO_CHASE )
        {
            diagBlock += grid.Width();
        }
        else
        {
            interactionList.push_back(interaction);
            // We must take care to not participate twice in one chase
            if( grid.Width() == 1 )
                diagBlock = interaction.block1 + 1;
            else
                diagBlock += grid.Width();
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
void CollectBlock
( const InterBlockInteraction& interaction,
  const DistMatrix<F,MC,MR,BLOCK>& H,
        Matrix<F>& HBlock,
  const DistChaseState& state )
{
    DEBUG_CSE
    const auto& HLoc = H.LockedMatrix();
    const Grid& grid = H.Grid();

    // (Up to) four processes may participate in introducing and chasing
    // a packet into and out of a partial block.
    const int firstRow =
      Mod( state.winColAlign+interaction.block0, grid.Height() );
    const int firstCol =
      Mod( state.winRowAlign+interaction.block0, grid.Width() );
    const int secondRow = Mod( firstRow+1, grid.Height() );
    const int secondCol = Mod( firstCol+1, grid.Width() );
    DEBUG_ONLY(
      if( (grid.Row() != firstRow && grid.Row() != secondRow) ||
          (grid.Col() != firstCol && grid.Col() != secondCol) )
          LogicError("This process does not participate in this interaction");
    )

    if( interaction.chaseType == SIMPLE_INTRO_CHASE )
    {
        // Only a single process participates in introductory chases, and they
        // occur over the entire top-left block (which must have been full).
        if( grid.Row() == secondRow && grid.Col() == secondCol )
        {
            const Int indexBeg = state.winBeg;
            const Int indexEnd = state.winBeg + interaction.blockSize1;
            const Int localRowBeg = H.LocalRowOffset( indexBeg );
            const Int localRowEnd = H.LocalRowOffset( indexEnd );
            const Int localColBeg = H.LocalColOffset( indexBeg );
            const Int localColEnd = H.LocalColOffset( indexEnd );
            auto HInteractLoc =
              HLoc( IR(localRowBeg,localRowEnd), IR(localColBeg,localColEnd) );
            HBlock = HInteractLoc;
        }
        else
        {
            LogicError("Invalid SIMPLE_INTRO_CHASE StoreBlock");
        }
        return;
    }
    else if( interaction.chaseType == NO_CHASE )
        LogicError("Invalid request to collect an inter-block window");

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

template<typename F>
void StoreBlock
( const InterBlockInteraction& interaction,
        DistMatrix<F,MC,MR,BLOCK>& H,
  const Matrix<F>& HBlock,
  const DistChaseState& state )
{
    DEBUG_CSE
    const auto& HLoc = H.LockedMatrix();
    const Grid& grid = H.Grid();

    // (Up to) four processes may participate in introducing and chasing
    // a packet into and out of a partial block.
    const int firstRow =
      Mod( state.winColAlign+interaction.block0, grid.Height() );
    const int firstCol =
      Mod( state.winRowAlign+interaction.block0, grid.Width() );
    const int secondRow = Mod( firstRow+1, grid.Height() );
    const int secondCol = Mod( firstCol+1, grid.Width() );
    DEBUG_ONLY(
      if( (grid.Row() != firstRow && grid.Row() != secondRow) ||
          (grid.Col() != firstCol && grid.Col() != secondCol) )
          LogicError("This process does not participate in this interaction");
    )

    if( interaction.chaseType == SIMPLE_INTRO_CHASE )
    {
        // Only a single process participates in introductory chases, and they
        // occur over the entire top-left block (which must have been full).
        if( grid.Row() == secondRow && grid.Col() == secondCol )
        { 
            const Int indexBeg = state.winBeg;
            const Int indexEnd = state.winBeg + interaction.blockSize1;
            const Int localRowBeg = H.LocalRowOffset( indexBeg );
            const Int localRowEnd = H.LocalRowOffset( indexEnd );
            const Int localColBeg = H.LocalColOffset( indexBeg );
            const Int localColEnd = H.LocalColOffset( indexEnd );
            auto HInteractLoc =
              HLoc( IR(localRowBeg,localRowEnd), IR(localColBeg,localColEnd) );
            HInteractLoc = HBlock;
        }
        else
        {
            LogicError("Invalid SIMPLE_INTRO_CHASE StoreBlock");
        }
        return;
    }
    else if( interaction.chaseType == NO_CHASE )
        LogicError("Invalid request to collect an inter-block window");

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

    if( grid.Height() == 1 && grid.Width() == 1 )
    {
        // Only our process participates
        HInteractLoc = HBlock;
    }
    else if( grid.Height() == 1 )
    {
        // Two processes in the same row participate
        if( grid.Col() == firstCol )
        {
            auto HBlockLeft = HBlock( ALL, ind0 );
            HInteractLoc = HBlockLeft;
        }
        else
        {
            auto HBlockRight = HBlock( ALL, ind1 );
            HInteractLoc = HBlockRight;
        }
    }
    else if( grid.Width() == 1 )
    {
        // Two processes in the same column participate
        auto HBlockTop = HBlock( ind0, ALL );
        auto HBlockBottom = HBlock( ind1, ALL );
        if( grid.Row() == firstRow )
        {
            HInteractLoc = HBlockTop;
        }
        else
        {
            HInteractLoc = HBlockBottom;
        }
    }
    else
    {
        // Four processes participate, though only the upper-left and
        // bottom-right ones will chase the packet, so only they receive
        // any data.
        if( grid.Row() == firstRow && grid.Col() == firstCol )
        {
            auto HBlock00 = HBlock( ind0, ind0 );
            HInteractLoc = HBlock00;
        }
        else if( grid.Row() == firstRow && grid.Col() == secondCol )
        {
            auto HBlock01 = HBlock( ind0, ind1 );
            HInteractLoc = HBlock01;
        }
        else if( grid.Row() == secondRow && grid.Col() == firstCol )
        {
            auto HBlock10 = HBlock( ind1, ind0 );
            HInteractLoc = HBlock10;
        }
        else if( grid.Row() == secondRow && grid.Col() == secondCol )
        {
            auto HBlock11 = HBlock( ind1, ind1 );
            HInteractLoc = HBlock11;
        }
    }
}

template<typename F>
void LocalChase
( bool evenToOdd,
  const InterBlockInteraction& interaction,
        Matrix<F>& HBlock,
        Matrix<F>& UBlock,
        Matrix<F>& W,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
  const DistChaseState& state,
        bool progress )
{
    DEBUG_CSE
    const auto& shiftsLoc = shifts.LockedMatrix();
    const Int blockWinBeg = 0;
    const Int blockWinEnd = interaction.blockSize0 + interaction.blockSize1;
    Identity( UBlock, blockWinEnd, blockWinEnd );

    const Int householderSize =
      interaction.householderEnd - interaction.householderBeg;

    // If this is not an exit chase, each of the Householder transformations is
    // 3x3 and overlaps in one entry with the previous transformation(s).
    // If this *is* an exit chase, the last transformation is 2x2 and again
    // overlaps in one entry.
    const Int numSteps = 
      ( interaction.chaseType == EXIT_CHASE ?
        householderSize - 1 :
        householderSize - 2 );

    // All non-exit blocks can carry a full load of shifts, with the exception
    // of non-full first diagonal blocks. Further, the block indices of a
    // non-full introductory chase are (0,1), whereas they are (-1,0) for a full
    // introductory chase.
    //
    // Let us consider the four scenarios: the first block is either full or
    // non-full, and the chase is either of the same or different parity. The
    // following diagrams mark the sequences of interactions with the (maximum) 
    // number of packets that will live in each at the end of the chase.
    //
    // Full, Same parity:
    //
    //  (0,1), (2,3), (4,5), ...
    //    2      2      2
    //
    // Non-full, Same parity:
    //
    //  (0,1), (2,3), (4,5), ...
    //    2      2      2
    //
    // Full, Different parity:
    //
    //  (-1,0), (1,2), (3,4), ...
    //     2      2      2 
    //
    // Non-full, Different parity:
    //
    //  (1,2), (3,4), (5,6), ...
    //    2      2      2
    //
    const bool fullFirstBlock = ( state.firstBlockSize == state.blockSize );
    const bool evenFirst = ( Mod( state.activeBlockEnd, 2 ) == 0 );
    const bool sameParity = ( evenFirst == evenToOdd );
    Int packetOffset;
    if( sameParity )
        packetOffset = interaction.block0;
    else if( fullFirstBlock )
        packetOffset = interaction.block0 + 1;
    else
        packetOffset = interaction.block0 - 1;

    const Int firstBulge = state.bulgeBeg + packetOffset;

    Matrix<F> ZDummy;
    const Int chaseBeg = interaction.householderBeg-1;
    const Int transformRowBeg = blockWinBeg;
    const Int transformColEnd = blockWinEnd;
    const bool wantSchurVecsSub = false;
    const bool accumulateSub = true;
    for( Int step=0; step<numSteps; ++step )
    {
        const Int packetBeg = chaseBeg + step;
        ComputeReflectors
        ( HBlock, blockWinBeg, blockWinEnd, shiftsLoc, W, packetBeg,
          firstBulge, interaction.numBulges, progress );
        ApplyReflectorsOpt
        ( HBlock, blockWinBeg, blockWinEnd, chaseBeg, packetBeg,
          transformRowBeg, transformColEnd, ZDummy, wantSchurVecsSub,
          UBlock, W, firstBulge, accumulateSub, progress );
    }
}

template<typename F>
void ApplyAccumulatedFromLeft
( const InterBlockInteraction& interaction,
        DistMatrix<F,MC,MR,BLOCK>& H,
  const Matrix<F>& U,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    auto& HLoc = H.Matrix();

    // HRight := U' HRight
    {
        const Int colBeg = interaction.end; 
        const Int colEnd = ( ctrl.fullTriangle ? n : state.winEnd );
        const Int localRowBeg = H.LocalRowOffset( interaction.householderBeg ); 
        const Int localRowEnd = H.LocalRowOffset( interaction.householderEnd );
        const Int localColBeg = H.LocalColOffset( colBeg );
        const Int localColEnd = H.LocalColOffset( colEnd );
        auto HRightLoc =
          HLoc( IR(localRowBeg,localRowEnd), IR(localColBeg,localColEnd) );
        auto HRightCopy( HRightLoc );
        Gemm( ADJOINT, NORMAL, F(1), U, HRightCopy, HRightLoc );
    }
}

template<typename F>
void ApplyAccumulatedFromRight
( const InterBlockInteraction& interaction,
        DistMatrix<F,MC,MR,BLOCK>& H,
        DistMatrix<F,MC,MR,BLOCK>& Z,
  const Matrix<F>& U,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    auto& HLoc = H.Matrix();
    auto& ZLoc = Z.Matrix();

    // HTop := HTop U
    {
        const Int rowBeg = ( ctrl.fullTriangle ? 0 : state.winBeg );
        const Int rowEnd = interaction.beg;
        const Int localRowBeg = H.LocalRowOffset( rowBeg ); 
        const Int localRowEnd = H.LocalRowOffset( rowEnd );
        const Int localColBeg = H.LocalColOffset( interaction.householderBeg );
        const Int localColEnd = H.LocalColOffset( interaction.householderEnd );
        auto HTopLoc =
          HLoc( IR(localRowBeg,localRowEnd), IR(localColBeg,localColEnd) );
        auto HTopCopy( HTopLoc );
        Gemm( NORMAL, NORMAL, F(1), HTopCopy, U, HTopLoc );
    }

    // ZInteract := ZInteract U
    if( ctrl.wantSchurVecs )
    {
        const Int localColBeg = Z.LocalColOffset( interaction.householderBeg );
        const Int localColEnd = Z.LocalColOffset( interaction.householderEnd );
        auto ZInteractLoc = ZLoc( ALL, IR(localColBeg,localColEnd) );
        auto ZInteractCopy( ZInteractLoc );
        Gemm( NORMAL, NORMAL, F(1), ZInteractCopy, U, ZInteractLoc );
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
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const auto& shiftsLoc = shifts.LockedMatrix();
    const Grid& grid = H.Grid();
    const bool fullFirstBlock = ( state.firstBlockSize == state.blockSize );
    // If fullFirstBlock is false, then we need to subtract one from the block
    // index when computing the beginning shift.
    
    if( state.activeBlockBeg < 0 )
        LogicError("state.activeBlockBeg was negative");
    if( state.activeBlockBeg > state.activeBlockEnd )
        LogicError("state.activeBlockBeg > state.activeBlockEnd");
    if( state.activeBlockEnd == 1 && !fullFirstBlock )
        LogicError("Cannot introduce any bulges");

    auto rowInteractionList =
      interblock::FormRowInteractionList( evenToOdd, grid, state );
    auto colInteractionList =
      interblock::FormColumnInteractionList( evenToOdd, grid, state );

    const Int numRowInteractions = rowInteractionList.size();
    const Int numColInteractions = colInteractionList.size();

    // Count the number of interactions our process participates in
    Int numLocalInteractions = 0;
    for( const auto& interaction : rowInteractionList )
        if( interaction.participating ) 
            ++numLocalInteractions;
    vector<Matrix<F>> UList(numLocalInteractions);

    // Chase the packets that we interact with in this step and store the
    // accumulated Householder reflections
    Matrix<F> W;
    Matrix<F> HBlock;
    Int localInteraction = 0;
    for( Int rowInteraction=0; rowInteraction<numRowInteractions;
         ++rowInteraction )
    {
        auto interaction = rowInteractionList[rowInteraction];
        if( interaction.participating )
        {
            auto& UBlock = UList[localInteraction];
            interblock::CollectBlock( interaction, H, HBlock, state );
            interblock::LocalChase
            ( evenToOdd, interaction, HBlock, UBlock, W, shifts, state,
              ctrl.progress );
            interblock::StoreBlock( interaction, H, HBlock, state );
            ++localInteraction;
        }
    }

    localInteraction = 0;
    Matrix<F> U;
    for( Int rowInteraction=0; rowInteraction<numRowInteractions;
         ++rowInteraction )
    {
        auto interaction = rowInteractionList[rowInteraction];
        const Int interactionSize =
          interaction.blockSize0 + interaction.blockSize1;
        if( interaction.participating )
            U = UList[localInteraction++];
        else
            Zeros( U, interactionSize, interactionSize );

        const int firstRow =
          Mod( state.winColAlign+interaction.block0, grid.Height() );
        const int firstCol =
          Mod( state.winRowAlign+interaction.block0, grid.Width() );
        const int secondRow = Mod( firstRow+1, grid.Height() );
        const int secondCol = Mod( firstCol+1, grid.Width() );

        int ownerCol;
        if( interaction.chaseType == SIMPLE_INTRO_CHASE )
            ownerCol = secondCol; 
        else if( firstRow == grid.Row() )
            ownerCol = firstCol;
        else
            ownerCol = secondCol;

        El::Broadcast( U, grid.RowComm(), ownerCol );

        interblock::ApplyAccumulatedFromLeft( interaction, H, U, state, ctrl );
    }

    localInteraction = 0;
    for( Int colInteraction=0; colInteraction<numColInteractions;
         ++colInteraction )
    {
        auto interaction = colInteractionList[colInteraction];
        const Int interactionSize =
          interaction.blockSize0 + interaction.blockSize1;
        if( interaction.participating )
            U = UList[localInteraction++];
        else
            Zeros( U, interactionSize, interactionSize );

        const int firstRow =
          Mod( state.winColAlign+interaction.block0, grid.Height() );
        const int firstCol =
          Mod( state.winRowAlign+interaction.block0, grid.Width() );
        const int secondRow = Mod( firstRow+1, grid.Height() );
        const int secondCol = Mod( firstCol+1, grid.Width() );

        int ownerRow;
        if( interaction.chaseType == SIMPLE_INTRO_CHASE )
            ownerRow = secondRow; 
        else if( firstCol == grid.Col() )
            ownerRow = firstRow;
        else
            ownerRow = secondRow;

        El::Broadcast( U, grid.ColComm(), ownerRow );

        interblock::ApplyAccumulatedFromRight
        ( interaction, H, Z, U, state, ctrl );
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTER_BLOCK_CHASE_HPP
