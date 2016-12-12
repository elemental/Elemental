/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTER_BLOCK_CHASE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTER_BLOCK_CHASE_HPP

#include "../../Transform.hpp"

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
  bool onDiagonal;
  bool chasing;
};

namespace interblock {

inline InterBlockInteraction
DetermineInteraction
( bool evenToOdd,
  Int diagBlockRow,
  const Grid& grid,
  const DistChaseState& state )
{
    EL_DEBUG_CSE
    if( diagBlockRow < Max(state.introBlock,0) ||
        diagBlockRow >= state.endBlock )
        LogicError("Diagonal block row was not in the active range");

    const bool fullFirstBlock = ( state.firstBlockSize == state.blockSize );

    const Int distFromEnd = state.endBlock - diagBlockRow;
    const bool evenBlock = ( Mod( distFromEnd, 2 ) == 0 );
    const bool sameParity = ( evenBlock == evenToOdd );

    InterBlockInteraction interaction;
    interaction.chaseType = NO_CHASE;
    interaction.block0 = ( sameParity ? diagBlockRow : diagBlockRow-1 );
    interaction.block1 = ( sameParity ? diagBlockRow+1 : diagBlockRow );
    interaction.blockSize0 = 0;
    interaction.blockSize1 = 0;
    interaction.participating = false;
    interaction.onDiagonal = false;
    interaction.chasing = false;
    if( diagBlockRow == state.introBlock && !sameParity )
    {
        // We cannot pull into the intro block
        return interaction;
    }
    if( diagBlockRow == state.endBlock-1 && sameParity )
    {
        // Packets are never left in the next-to-last block, so we cannot push
        // a packet from said position
        return interaction;
    }

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
                return interaction;
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
    // state.endBlock-1.
    interaction.numBulges =
      ( interaction.block1 < state.endBlock-1 ?
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
    if( interaction.chaseType == SIMPLE_INTRO_CHASE )
    {
        interaction.participating =
          (grid.Row() == secondRow && grid.Col() == secondCol);
        interaction.onDiagonal = interaction.participating;
        interaction.chasing = interaction.participating;
    }
    else
    {
        const bool inTwoByTwo =
          (grid.Row() == firstRow || grid.Row() == secondRow) &&
          (grid.Col() == firstCol || grid.Col() == secondCol);
        interaction.participating = inTwoByTwo;
        interaction.onDiagonal =
          (grid.Row() == firstRow && grid.Col() == firstCol) ||
          (grid.Row() == secondRow && grid.Col() == secondCol);
        interaction.chasing =
          (grid.Row() == firstRow && grid.Col() == firstCol);
    }

    return interaction;
}

inline vector<InterBlockInteraction>
FormRowInteractionList
( bool evenToOdd,
  const Grid& grid,
  const DistChaseState& state )
{
    EL_DEBUG_CSE
    vector<InterBlockInteraction> interactionList;
    // Only loop over the blocks that are assigned to our process row
    // and occur within the active window.
    Int diagBlock = state.activeRowBlockBeg;
    while( diagBlock < state.endBlock )
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
    EL_DEBUG_CSE
    vector<InterBlockInteraction> interactionList;
    // Only loop over the blocks that are assigned to our process column
    // and occur within the active window.
    Int diagBlock = state.activeColBlockBeg;
    while( diagBlock < state.endBlock )
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
template<typename Field>
void CollectBlock
( const InterBlockInteraction& interaction,
  const DistMatrix<Field,MC,MR,BLOCK>& H,
        Matrix<Field>& HBlock,
  const DistChaseState& state )
{
    EL_DEBUG_CSE
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
    EL_DEBUG_ONLY(
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
            LogicError("Invalid SIMPLE_INTRO_CHASE in CollectBlock");
        }
        return;
    }
    if( interaction.chaseType == NO_CHASE )
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
        HBlock = HInteractLoc;
    }
    else if( grid.Height() == 1 )
    {
        // Only the first of two processes will chase
        if( grid.Col() == firstCol )
        {
            Zeros( HBlock, interactionSize, interactionSize );
            auto HBlockLeft = HBlock( ALL, ind0 );
            auto HBlockRight = HBlock( ALL, ind1 );

            HBlockLeft = HInteractLoc;
            El::Recv( HBlockRight, grid.RowComm(), secondCol );
        }
        else
        {
            El::Send( HInteractLoc, grid.RowComm(), firstCol );
        }
    }
    else if( grid.Width() == 1 )
    {
        // Only the first of two processes will chase
        if( grid.Row() == firstRow )
        {
            Zeros( HBlock, interactionSize, interactionSize );
            auto HBlockTop = HBlock( ind0, ALL );
            auto HBlockBottom = HBlock( ind1, ALL );

            HBlockTop = HInteractLoc;
            El::Recv( HBlockBottom, grid.ColComm(), secondRow );
        }
        else
        {
            El::Send( HInteractLoc, grid.ColComm(), firstRow );
        }
    }
    else
    {
        // Four processes participate, though only the upper-left will chase
        // We will use the column-major ordering (which is the VC comm.)
        const int proc00 = firstRow + firstCol*grid.Height();
        const int proc01 = firstRow + secondCol*grid.Height();
        const int proc10 = secondRow + firstCol*grid.Height();
        const int proc11 = secondRow + secondCol*grid.Height();
        if( grid.Row() == firstRow && grid.Col() == firstCol )
        {
            Zeros( HBlock, interactionSize, interactionSize );
            auto HBlock00 = HBlock( ind0, ind0 );
            auto HBlock01 = HBlock( ind0, ind1 );
            auto HBlock10 = HBlock( ind1, ind0 );
            auto HBlock11 = HBlock( ind1, ind1 );

            HBlock00 = HInteractLoc;
            El::Recv( HBlock01, grid.VCComm(), proc01 );
            El::Recv( HBlock10, grid.VCComm(), proc10 );
            El::Recv( HBlock11, grid.VCComm(), proc11 );
        }
        else if( grid.Row() == firstRow && grid.Col() == secondCol )
        {
            El::Send( HInteractLoc, grid.VCComm(), proc00 );
        }
        else if( grid.Row() == secondRow && grid.Col() == firstCol )
        {
            El::Send( HInteractLoc, grid.VCComm(), proc00 );
        }
        else if( grid.Row() == secondRow && grid.Col() == secondCol )
        {
            El::Send( HInteractLoc, grid.VCComm(), proc00 );
        }
    }
}

template<typename Field>
void StoreBlock
( const InterBlockInteraction& interaction,
        DistMatrix<Field,MC,MR,BLOCK>& H,
        Matrix<Field>& HBlock,
  const DistChaseState& state )
{
    EL_DEBUG_CSE
    auto& HLoc = H.Matrix();
    const Grid& grid = H.Grid();

    // (Up to) four processes may participate in introducing and chasing
    // a packet into and out of a partial block.
    const int firstRow =
      Mod( state.winColAlign+interaction.block0, grid.Height() );
    const int firstCol =
      Mod( state.winRowAlign+interaction.block0, grid.Width() );
    const int secondRow = Mod( firstRow+1, grid.Height() );
    const int secondCol = Mod( firstCol+1, grid.Width() );
    EL_DEBUG_ONLY(
      if( (grid.Row() != firstRow && grid.Row() != secondRow) ||
          (grid.Col() != firstCol && grid.Col() != secondCol) )
          LogicError("This process does not participate in this interaction");
    )

    // We can grab the indices of our local portion of the interaction window
    // in a black-box manner
    const Int localRowBeg = H.LocalRowOffset( interaction.beg );
    const Int localRowEnd = H.LocalRowOffset( interaction.end );
    const Int localColBeg = H.LocalColOffset( interaction.beg );
    const Int localColEnd = H.LocalColOffset( interaction.end );
    auto HInteractLoc =
      HLoc( IR(localRowBeg,localRowEnd), IR(localColBeg,localColEnd) );

    if( interaction.chaseType == SIMPLE_INTRO_CHASE )
    {
        // Only a single process participates in introductory chases, and they
        // occur over the entire top-left block (which must have been full).
        if( grid.Row() == secondRow && grid.Col() == secondCol )
        {
            HInteractLoc = HBlock;
        }
        else
        {
            LogicError("Invalid SIMPLE_INTRO_CHASE StoreBlock");
        }
        return;
    }
    if( interaction.chaseType == NO_CHASE )
        LogicError("Invalid request to collect an inter-block window");

    const Int interactionSize = interaction.blockSize0 + interaction.blockSize1;
    const auto ind0 = IR(0,interaction.blockSize0);
    const auto ind1 = IR(interaction.blockSize0,interactionSize);

    if( grid.Height() == 1 && grid.Width() == 1 )
    {
        HInteractLoc = HBlock;
    }
    else if( grid.Height() == 1 )
    {
        if( grid.Col() == firstCol )
        {
            auto HBlockLeft = HBlock( ALL, ind0 );
            auto HBlockRight = HBlock( ALL, ind1 );
            HInteractLoc = HBlockLeft;
            El::Send( HBlockRight, grid.RowComm(), secondCol );
        }
        else
        {
            El::Recv( HInteractLoc, grid.RowComm(), firstCol );
        }
    }
    else if( grid.Width() == 1 )
    {
        if( grid.Row() == firstRow )
        {
            auto HBlockTop = HBlock( ind0, ALL );
            auto HBlockBottom = HBlock( ind1, ALL );
            HInteractLoc = HBlockTop;
            El::Send( HBlockBottom, grid.ColComm(), secondRow );
        }
        else
        {
            El::Recv( HInteractLoc, grid.ColComm(), firstRow );
        }
    }
    else
    {
        // Four processes participate, though only the upper-left and
        // bottom-right ones will chase the packet, so only they receive
        // any data.
        const int proc00 = firstRow + firstCol*grid.Height();
        const int proc01 = firstRow + secondCol*grid.Height();
        const int proc10 = secondRow + firstCol*grid.Height();
        const int proc11 = secondRow + secondCol*grid.Height();
        if( grid.Row() == firstRow && grid.Col() == firstCol )
        {
            auto HBlock00 = HBlock( ind0, ind0 );
            auto HBlock01 = HBlock( ind0, ind1 );
            auto HBlock10 = HBlock( ind1, ind0 );
            auto HBlock11 = HBlock( ind1, ind1 );
            HInteractLoc = HBlock00;
            El::Send( HBlock01, grid.VCComm(), proc01 );
            El::Send( HBlock10, grid.VCComm(), proc10 );
            El::Send( HBlock11, grid.VCComm(), proc11 );
        }
        else if( grid.Row() == firstRow && grid.Col() == secondCol )
        {
            El::Recv( HInteractLoc, grid.VCComm(), proc00 );
        }
        else if( grid.Row() == secondRow && grid.Col() == firstCol )
        {
            El::Recv( HInteractLoc, grid.VCComm(), proc00 );
        }
        else if( grid.Row() == secondRow && grid.Col() == secondCol )
        {
            El::Recv( HInteractLoc, grid.VCComm(), proc00 );
        }
    }
}

template<typename Field>
void LocalChase
( bool evenToOdd,
  const InterBlockInteraction& interaction,
        Matrix<Field>& HBlock,
        Matrix<Field>& UBlock,
        Matrix<Field>& W,
  const DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
  const DistChaseState& state,
        bool progress )
{
    EL_DEBUG_CSE
    const auto& shiftsLoc = shifts.LockedMatrix();
    const Int blockWinBeg = 0;
    const Int blockWinEnd = interaction.blockSize0 + interaction.blockSize1;

    const Int householderSize =
      interaction.householderEnd - interaction.householderBeg;
    Identity( UBlock, householderSize, householderSize );
    Zeros( W, 3, interaction.numBulges );

    Int numSteps;
    if( interaction.chaseType == STANDARD_CHASE )
    {
        // Standard chases involve stepHouseholderSize x stepHouseholderSize
        // transformations; the effected index range expands by one in each step
        const Int stepHouseholderSize = 3*interaction.numBulges;
        numSteps = householderSize - stepHouseholderSize + 1;
    }
    else if( interaction.chaseType == EXIT_CHASE )
    {
        // Exit chases involve a 2x2 in the last step and expand by one with
        // each previous transformation
        numSteps = householderSize - 1;
    }
    else
    {
        // Introductory chases involve 3x3 rotations in the first step and
        // expand by one in each subsequent step
        numSteps = householderSize - 2;
    }

    // TODO(poulson): Description of the following
    const bool fullFirstBlock = ( state.firstBlockSize == state.blockSize );
    Int packetOffset;
    if( state.introBlock == - 1 && !fullFirstBlock )
    {
        packetOffset =
          (interaction.block0-(state.introBlock+1))*state.numBulgesPerBlock;
    }
    else
    {
        packetOffset =
          (interaction.block0-state.introBlock)*state.numBulgesPerBlock;
    }
    const Int bulgeOffset = state.bulgeBeg + packetOffset;

    Matrix<Field> ZDummy;
    const Int chaseBeg = (interaction.householderBeg-interaction.beg)-1;
    const Int transformRowBeg = blockWinBeg;
    const Int transformColEnd = blockWinEnd;
    const bool wantSchurVecsSub = false;
    const bool accumulateSub = true;
    for( Int step=0; step<numSteps; ++step )
    {
        const Int firstActiveBulgePosition = 0;
        Int packetBeg, numActiveBulges, firstActiveBulge;
        if( interaction.chaseType == SIMPLE_INTRO_CHASE )
        {
            // One bulge is introduced every three steps
            packetBeg = chaseBeg + Mod(step,3);
            numActiveBulges = (step/3) + 1;
            firstActiveBulge = interaction.numBulges - numActiveBulges;
        }
        else if( interaction.chaseType == COUPLED_INTRO_CHASE )
        {
            // At most one bulge is introduced every three steps;
            // recall that step-1 is the starting position of the right-most
            // bulge
            packetBeg =
              Max(chaseBeg+Mod(step,3),(step-1)-3*(interaction.numBulges-1));
            numActiveBulges = Min( (step/3)+1, interaction.numBulges );
            firstActiveBulge = interaction.numBulges - numActiveBulges;
        }
        else if( interaction.chaseType == STANDARD_CHASE )
        {
            // All bulges are active
            packetBeg = chaseBeg + step;
            numActiveBulges = interaction.numBulges;
            firstActiveBulge = 0;
        }
        else
        {
            // At most one bulge is removed every three steps
            packetBeg = chaseBeg + step;
            numActiveBulges = Min( (numSteps-step+2)/3, interaction.numBulges );
            firstActiveBulge = 0;
        }

        const IR activeInd(firstActiveBulge,firstActiveBulge+numActiveBulges);
        const auto& activeShifts =
          shiftsLoc( IR(2*activeInd.beg,2*activeInd.end)+(2*bulgeOffset), ALL );

        ComputeReflectors
        ( HBlock, blockWinBeg, blockWinEnd, activeShifts, W, packetBeg,
          firstActiveBulgePosition, numActiveBulges, progress );
        ApplyReflectorsOpt
        ( HBlock, blockWinBeg, blockWinEnd, chaseBeg, packetBeg,
          transformRowBeg, transformColEnd, ZDummy, wantSchurVecsSub,
          UBlock, W, firstActiveBulgePosition, numActiveBulges, accumulateSub,
          progress );
    }
}

template<typename Field>
void ApplyAccumulatedFromLeft
( const InterBlockInteraction& interaction,
        DistMatrix<Field,MC,MR,BLOCK>& H,
  const Matrix<Field>& U,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    const Int colBeg = interaction.end;
    const Int colEnd = ( ctrl.fullTriangle ? n : state.winEnd );
    const Int houseBeg = interaction.householderBeg;
    const Int houseEnd = interaction.householderEnd;
    EL_DEBUG_ONLY(
      if( houseEnd-houseBeg != U.Height() )
          LogicError
          ("U was of size ",U.Height()," but householder indices are [",
           houseBeg,",",houseEnd,")");
    )

    // HRight := U' HRight
    auto HRight = H( IR(houseBeg,houseEnd), IR(colBeg,colEnd) );
    TransformRows( U, HRight );
}

template<typename Field>
void ApplyAccumulatedFromRight
( const InterBlockInteraction& interaction,
        DistMatrix<Field,MC,MR,BLOCK>& H,
        DistMatrix<Field,MC,MR,BLOCK>& Z,
  const Matrix<Field>& U,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int rowBeg = ( ctrl.fullTriangle ? 0 : state.winBeg );
    const Int rowEnd = interaction.beg;
    const Int houseBeg = interaction.householderBeg;
    const Int houseEnd = interaction.householderEnd;
    EL_DEBUG_ONLY(
      if( houseEnd-houseBeg != U.Height() )
          LogicError
          ("U was of size ",U.Height()," but householder indices are [",
           houseBeg,",",houseEnd,")");
      if( Z.RowAlign() != H.RowAlign() )
          LogicError("Z and H had different row alignments");
      if( Z.RowCut() != H.RowCut() )
          LogicError("Z and H had different row cuts");
    )

    auto HTop = H( IR(rowBeg,rowEnd), IR(houseBeg,houseEnd) );
    TransformColumns( U, HTop );

    // ZBlock := ZBlock U
    auto ZBlock = Z( ALL, IR(houseBeg,houseEnd) );
    TransformColumns( U, ZBlock );
}

} // namespace interblock

template<typename Field>
void InterBlockChase
(       DistMatrix<Field,MC,MR,BLOCK>& H,
        DistMatrix<Field,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
        bool evenToOdd,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = H.Grid();
    const bool fullFirstBlock = ( state.firstBlockSize == state.blockSize );
    // If fullFirstBlock is false, then we need to subtract one from the block
    // index when computing the beginning shift.

    if( state.introBlock < -1 )
        LogicError("state.introBlock was less than -1");
    if( state.introBlock >= state.endBlock )
        LogicError("state.introBlock >= state.endBlock");
    if( state.endBlock == 1 && !fullFirstBlock )
        LogicError("Cannot introduce any bulges");

    auto rowInteractionList =
      interblock::FormRowInteractionList( evenToOdd, grid, state );
    auto colInteractionList =
      interblock::FormColumnInteractionList( evenToOdd, grid, state );

    const Int numRowInteractions = rowInteractionList.size();
    const Int numColInteractions = colInteractionList.size();

    // Count the number of chases handled by this process
    Int numDiagInteractions = 0;
    for( const auto& interaction : rowInteractionList )
        if( interaction.onDiagonal )
            ++numDiagInteractions;
    vector<Matrix<Field>> UList(numDiagInteractions);

    const int prevGridRow = Mod( grid.Row()-1, grid.Height() );
    const int nextGridRow = Mod( grid.Row()+1, grid.Height() );

    const int prevGridCol = Mod( grid.Col()-1, grid.Width() );
    const int nextGridCol = Mod( grid.Col()+1, grid.Width() );

    const int prevDiagProc = prevGridRow + prevGridCol*grid.Height();
    const int nextDiagProc = nextGridRow + nextGridCol*grid.Height();

    // Chase the packets that we interact with in this step and store the
    // accumulated Householder reflections
    Matrix<Field> W;
    Matrix<Field> HBlock;
    Int diagInteraction = 0;
    for( Int rowInteraction=0; rowInteraction<numRowInteractions;
         ++rowInteraction )
    {
        auto interaction = rowInteractionList[rowInteraction];
        if( interaction.participating )
        {
            interblock::CollectBlock( interaction, H, HBlock, state );
            if( interaction.chasing )
            {
                auto& UBlock = UList[diagInteraction++];
                interblock::LocalChase
                ( evenToOdd, interaction, HBlock, UBlock, W, shifts, state,
                  ctrl.progress );
                if( interaction.chaseType != SIMPLE_INTRO_CHASE &&
                    (grid.Height() != 1 || grid.Width() != 1) )
                    El::Send( UBlock, grid.VCComm(), nextDiagProc );
            }
            else if( interaction.onDiagonal )
            {
                auto& UBlock = UList[diagInteraction++];
                const Int householderSize =
                  interaction.householderEnd - interaction.householderBeg;
                UBlock.Resize( householderSize, householderSize );
                El::Recv( UBlock, grid.VCComm(), prevDiagProc );
            }
            interblock::StoreBlock( interaction, H, HBlock, state );
        }
    }

    diagInteraction = 0;
    Matrix<Field> U;
    for( Int rowInteraction=0; rowInteraction<numRowInteractions;
         ++rowInteraction )
    {
        auto interaction = rowInteractionList[rowInteraction];
        const Int householderSize =
          interaction.householderEnd - interaction.householderBeg;
        if( interaction.onDiagonal )
            U = UList[diagInteraction++];
        else
            Zeros( U, householderSize, householderSize );
        EL_DEBUG_ONLY(
          if( U.Height() != householderSize || U.Width() != householderSize )
              LogicError
              ("U was ",U.Height()," x ",U.Width()," instead of ",
               householderSize," for row interaction ",rowInteraction);
        )

        const int firstRow =
          Mod( state.winColAlign+interaction.block0, grid.Height() );
        const int firstCol =
          Mod( state.winRowAlign+interaction.block0, grid.Width() );
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

    diagInteraction = 0;
    for( Int colInteraction=0; colInteraction<numColInteractions;
         ++colInteraction )
    {
        auto interaction = colInteractionList[colInteraction];
        const Int householderSize =
          interaction.householderEnd - interaction.householderBeg;
        if( interaction.onDiagonal )
            U = UList[diagInteraction++];
        else
            Zeros( U, householderSize, householderSize );
        EL_DEBUG_ONLY(
          if( U.Height() != householderSize || U.Width() != householderSize )
              LogicError
              ("U was ",U.Height()," x ",U.Width()," instead of ",
               householderSize," for column interaction ",colInteraction);
        )

        const int firstRow =
          Mod( state.winColAlign+interaction.block0, grid.Height() );
        const int firstCol =
          Mod( state.winRowAlign+interaction.block0, grid.Width() );
        const int secondRow = Mod( firstRow+1, grid.Height() );

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
