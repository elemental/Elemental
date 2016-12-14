/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_HPP

#include "./PairShifts.hpp"
#include "./Sweep/ComputeReflectors.hpp"
#include "./Sweep/ApplyReflectors.hpp"
#include "./Sweep/Dist.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename Field>
void SweepHelper
(       Matrix<Field>& H,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Field>& Z,
        Matrix<Field>& U,
        Matrix<Field>& W,
        Matrix<Field>& WAccum,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    EL_DEBUG_ONLY(
      if( winEnd-winBeg < 4 )
          LogicError
          ("multibulge::Sweep shouldn't be called for window sizes < 4");
    )

    const Int numShifts = shifts.Height();
    EL_DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts...");
      if( numShifts % 2 != 0 )
          LogicError("Expected an even number of shifts");
    )
    const Int numBulges = numShifts / 2;

    // Set aside space for storing either the three nonzero entries of the first
    // column of a quadratic polynomial for introducing each bulge or the scalar
    // 'tau' and non-unit subvector v such that the 2x2 or 3x3 Householder
    // reflection I - tau*[1;v] [1;v]' can safely be stored within a vector of
    // length 3 (following the LAPACK convention that 'tau' will lie in the
    // first entry and 'v' will begin in the second entry).
    W.Resize( 3, numBulges );

    // Each time a new 4x4 bulge is introduced into the upper-left corner
    // of the matrix, we think of the bulge as having started at index
    // (winBeg-1), as, after the introduction of the 3x3 Householder
    // similarity using the implicit Q theorem on (H-sigma_0 I)*(H-sigma_1 I),
    // the bulge has starting index winBeg.
    //
    // Initialize with the last bulge about to be introduced in the upper-left
    const Int sweepBeg = (winBeg-1) - 3*(numBulges-1);

    // The last bulge must be at least 3x3 in order to involve a 2x2
    // reflection, so it must start before winEnd-2
    const Int sweepEnd = winEnd-2;

    // Each chase of a packet involves shifting the first bulge one
    // column right of the starting position of the last bulge. This involves
    // shifting every bulge right 3*(numBulges-1) + 1 columns.
    const Int chaseStride = 3*(numBulges-1) + 1;

    // The total effected region of each movement is therefore (at most)
    // the span of numBulges 3x3 Householder reflections and the translation
    // distance of the packet (chaseStride)
    const Int maxSlabSize = 3*numBulges + chaseStride;

    EL_DEBUG_ONLY(
      if( H(winBeg+2,winBeg) != Field(0) )
          LogicError("H was not upper Hessenberg");
    )
    for( Int chaseBeg=sweepBeg; chaseBeg<sweepEnd; chaseBeg+=chaseStride )
    {
        const Int slabBeg = Max( chaseBeg, winBeg-1 );
        const Int slabEnd = Min( chaseBeg+maxSlabSize, winEnd );
        const Int slabSize = slabEnd - slabBeg;
        if( ctrl.accumulateReflections )
        {
            // The Householder transformations do not effect the first index
            Identity( U, slabSize-1, slabSize-1 );
        }

        Int innerTransformRowBeg;
        if( ctrl.accumulateReflections )
            innerTransformRowBeg = Max( winBeg, chaseBeg );
        else if( ctrl.fullTriangle )
            innerTransformRowBeg = 0;
        else
            innerTransformRowBeg = winBeg;

        Int innerTransformColEnd;
        if( ctrl.accumulateReflections )
            innerTransformColEnd = slabEnd;
        else if( ctrl.fullTriangle )
            innerTransformColEnd = n;
        else
            innerTransformColEnd = winEnd;

        const Int packetEnd = Min(chaseBeg+chaseStride,sweepEnd);
        for( Int packetBeg=chaseBeg; packetBeg<packetEnd; ++packetBeg )
        {
            const Int firstBulge = Max( 0, ((winBeg-1)-packetBeg+2)/3 );
            const Int numStepBulges =
              Min( numBulges, (winEnd-packetBeg)/3 ) - firstBulge;

            ComputeReflectors
            ( H, winBeg, winEnd, shifts, W, packetBeg, firstBulge,
              numStepBulges, ctrl.progress );

            ApplyReflectorsOpt
            ( H, winBeg, winEnd,
              chaseBeg, packetBeg,
              innerTransformRowBeg, innerTransformColEnd,
              Z, ctrl.wantSchurVecs, U, W,
              firstBulge, numStepBulges, ctrl.accumulateReflections,
              ctrl.progress );
        }

        if( ctrl.accumulateReflections )
        {
            const Int transformBeg = ( ctrl.fullTriangle ? 0 : winBeg );
            const Int transformEnd = ( ctrl.fullTriangle ? n : winEnd );

            // Horizontal far-from-diagonal application
            const Int rightIndBeg = slabEnd;
            const Int rightIndEnd = transformEnd;
            const auto rightInd = IR(rightIndBeg,rightIndEnd);
            const auto horzInd = IR( Max(chaseBeg+1,winBeg), slabEnd );
            auto HHorzFar = H( horzInd, rightInd );
            Gemm( ADJOINT, NORMAL, Field(1), U, HHorzFar, WAccum );
            HHorzFar = WAccum;

            // Vertical far-from-diagonal application
            auto vertInd = IR(transformBeg,Max(winBeg,chaseBeg));
            auto HVertFar = H( vertInd, horzInd );
            Gemm( NORMAL, NORMAL, Field(1), HVertFar, U, WAccum );
            HVertFar = WAccum;

            if( ctrl.wantSchurVecs )
            {
                auto ZSub = Z( ALL, horzInd );
                Gemm( NORMAL, NORMAL, Field(1), ZSub, U, WAccum );
                ZSub = WAccum;
            }
        }
    }
}

template<typename Field>
void SweepHelper
(       DistMatrix<Field,MC,MR,BLOCK>& H,
  const DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
        DistMatrix<Field,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    // TODO(poulson): Check that H is upper Hessenberg
    EL_DEBUG_ONLY(
      if( H.Height() < 2*H.BlockHeight() )
          LogicError("H spans less than two full blocks");
    )
    auto state = BuildDistChaseState( H, shifts, ctrl );

    while( state.bulgeEnd != 0 )
    {
        // Chase packets from the bottom-right corners of the even parity blocks
        // to the top-left corners of the odd parity blocks
        InterBlockChase( H, Z, shifts, true, state, ctrl );

        // Chase packets from the bottom-right corners of the odd parity blocks
        // to the top-left corners of the even parity blocks
        InterBlockChase( H, Z, shifts, false, state, ctrl );

        // Chase the packets from the top-left corners to the bottom-right
        // corners of each diagonal block.
        IntraBlockChase( H, Z, shifts, state, ctrl );

        AdvanceChaseState( H, state );
    }
}

template<typename Field>
void Sweep
(       Matrix<Field>& H,
        Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Field>& Z,
        Matrix<Field>& U,
        Matrix<Field>& W,
        Matrix<Field>& WAccum,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();

    const Int numShifts = shifts.Height();
    EL_DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts...");
    )
    if( numShifts % 2 != 0 )
        LogicError("Expected an even number of shifts");
    if( !IsComplex<Field>::value )
        PairShifts( shifts );

    const Int maxBulgesPerSweep = Max( n/6, 1 );
    const Int maxShiftsPerSweep = 2*maxBulgesPerSweep;
    for( Int shiftStart=0; shiftStart<numShifts; shiftStart+=maxShiftsPerSweep )
    {
        const Int numSweepShifts =
          Min( maxShiftsPerSweep, numShifts-shiftStart );
        auto sweepInd = IR(shiftStart,shiftStart+numSweepShifts);
        auto sweepShifts = shifts(sweepInd,ALL);
        SweepHelper( H, sweepShifts, Z, U, W, WAccum, ctrl );
    }
}

template<typename Field>
void Sweep
(       DistMatrix<Field,MC,MR,BLOCK>& H,
  const DistMatrix<Complex<Base<Field>>,STAR,STAR>& shifts,
        DistMatrix<Field,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( ctrl.wantSchurVecs && Z.DistData() != H.DistData() )
          LogicError("The distributions of H and Z should match");
    )
    // TODO(poulson): Check that H is upper Hessenberg
    const Int n = H.Height();
    EL_DEBUG_ONLY(
      if( n < 2*H.BlockHeight() )
          LogicError("H spans less than two full blocks");
    )

    const Int numShifts = shifts.Height();
    EL_DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts...");
    )
    if( numShifts % 2 != 0 )
        LogicError("Expected an even number of shifts");

    DistMatrix<Complex<Base<Field>>,STAR,STAR> validShifts(shifts);
    if( !IsComplex<Field>::value )
        PairShifts( validShifts.Matrix() );

    const Int maxBulgesPerSweep = Max( n/6, 1 );
    const Int maxShiftsPerSweep = 2*maxBulgesPerSweep;
    for( Int shiftStart=0; shiftStart<numShifts; shiftStart+=maxShiftsPerSweep )
    {
        const Int numSweepShifts =
          Min( maxShiftsPerSweep, numShifts-shiftStart );
        auto sweepInd = IR(shiftStart,shiftStart+numSweepShifts);
        auto sweepShifts = validShifts(sweepInd,ALL);
        SweepHelper( H, sweepShifts, Z, ctrl );
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_HPP
