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

// For the vast majority of the inter-block packet chases, the transformations
// should start in the bottom-right of one diagonal block and end in the
// top-left of the subsequent diagonal block. And, to support non-overlapping
// inter-block chases of "even" and "odd" diagonal blocks (with the parity
// determined by counting up from zero starting at the first active block), we
// require that the distribution block size be at least one plus six times the
// number of bulges.
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
//   ~  |     x | E E E E x x x x x x x x x |
//   ~  |       | E E E E x x x x x x x x x |
//   ~  |       | E E E E x x x x x x x x x |
//   ~  |       | E E E E E E E x x x x x x |
//   ~  |       |       E E E E x x x x x x |
//   ~  |       |       E E E E x x x x x x |,
//      |       |       E E E E x x x x x x |
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
//   ~ | x x x x x x x x x x x x x |       ~ | E E E E x x x x x x x x x |
//   ~ | x x x x x x x x x x x x x |       ~ | E E E E x x x x x x x x x |
//   ~ |   x x x x x x x x x x x x |       ~ | E E E E x x x x x x x x x |
//   ~ |     x x x x x x x x x x x |       ~ | E E E E E E E x x x x x x |
//   ~ |       x x x x x x x x x x |  |->  ~ |       E E E E x x x x x x |
//   ~ |         x x x x x x x x x |       ~ |       E E E E x x x x x x |
//     |           x x x x x x x x |         |       E E E E x x x x x x |,
//     |             x x x x x x x |         |             x x x x x x x |
//     |               x x x x x x |         |               x x x x x x |
//     |                 x x x x x |         |                 x x x x x |
//     |                   x x x x |         |                   x x x x |
//     |                     x x x |         |                     x x x |
//     |                       x x |         |                       x x |
//      ---------------------------           ---------------------------
//
// again, usually with two bulges of odd parity packed into the bottom-right
// corner.
//
// Lastly, if we chase a packet of bulges into the last diagonal block of the
// window, we go ahead and fully chase them out of the window. For example,
// we might begin in the form
//
//                      ~ ~ ~ ~ ~ ~   ~ ~ ~ ~ ~ ~
//       -----------------------------------------
//      | x x x x x x x x x x x x x | x x x x x x |
//      | x x x x x x x x x x x x x | x x x x x x |
//      |   x x x x x x x x x x x x | x x x x x x |
//      |     x x x x x x x x x x x | x x x x x x |
//      |       x x x x x x x x x x | x x x x x x |
//      |         x x x x x x x x x | x x x x x x |
//      |           x B B B B x x x | x x x x x x |
//    ~ |             B B B B x x x | x x x x x x |
//    ~ |             B B B B x x x | x x x x x x |
//    ~ |             B B B B B B B | x x x x x x |,
//    ~ |                   B B B B | x x x x x x |
//    ~ |                   B B B B | x x x x x x |
//    ~ |                   B B B B | x x x x x x |
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

template<typename F>
void InterBlockEvenChase
(       DistMatrix<F,MC,MR,BLOCK>& H,
        DistMatrix<F,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
        Matrix<F>& W,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE

    const Int blockSize = H.BlockHeight();
    if( blockSize != H.BlockWidth() )
        LogicError("InterBlockEvenChase assumes square distribution blocks");
    const Int blockCut = H.ColCut();
    if( blockCut != H.RowCut() )
        LogicError("InterBlockEvenChase assumes that the cuts are equal");

    // TODO(poulson): Handle 'introductory', standard, and 'exit' chases here
}

template<typename F>
void InterBlockOddChase
(       DistMatrix<F,MC,MR,BLOCK>& H,
        DistMatrix<F,MC,MR,BLOCK>& Z,
  const DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
        Matrix<F>& W,
  const DistChaseState& state,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE

    const Int blockSize = H.BlockHeight();
    if( blockSize != H.BlockWidth() )
        LogicError("InterBlockOddChase assumes square distribution blocks");
    const Int blockCut = H.ColCut();
    if( blockCut != H.RowCut() )
        LogicError("InterBlockOddChase assumes that the cuts are equal");

    // TODO(poulson): Handle standard and 'exit' chases here
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTER_BLOCK_CHASE_HPP
