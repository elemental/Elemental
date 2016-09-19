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

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename F>
void SweepHelper
(       Matrix<F>& H,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<F>& Z,
        Matrix<F>& U,
        Matrix<F>& W,
        Matrix<F>& WAccum,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    DEBUG_ONLY(
      if( winEnd-winBeg < 4 )
          LogicError
          ("multibulge::Sweep shouldn't be called for window sizes < 4"); 
    )

    const Int numShifts = shifts.Height();
    DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts..."); 
      if( numShifts % 2 != 0 )
          LogicError("Expected an even number of sweeps");
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
    const Int ghostBeg = (winBeg-1) - 3*(numBulges-1);

    // The last bulge must be at least 3x3 in order to involve a 2x2 
    // reflection, so it must start before winEnd-2
    const Int ghostEnd = winEnd-2;

    // Each movement of a packet involves shifting the first bulge one 
    // column right of the starting position of the last bulge. This involves
    // shifting every bulge right 3*(numBulges-1) + 1 columns.
    const Int ghostStride = 3*(numBulges-1) + 1;

    // The total effected region of each movement is therefore (at most)
    // the span of numBulges 3x3 Householder reflections and the translation
    // distance of the packet (ghostStride)
    const Int slabSize = 3*numBulges + ghostStride;

    DEBUG_ONLY(
      if( H(winBeg+2,winBeg) != F(0) )
          LogicError("H was not upper Hessenberg"); 
    )
    for( Int ghostCol=ghostBeg; ghostCol<ghostEnd; ghostCol+=ghostStride )
    {
        // Note that this slab endpoint may be past winEnd
        const Int slabEnd = ghostCol + slabSize;
        if( ctrl.accumulateReflections )
            Identity( U, slabSize-1, slabSize-1 );

        const Int packetEnd = Min(ghostCol+ghostStride,ghostEnd);
        for( Int packetBeg=ghostCol; packetBeg<packetEnd; ++packetBeg )
        {
            const Int firstBulge = Max( 0, ((winBeg-1)-packetBeg+2)/3 );
            const Int numFullBulges =
              Min( numBulges, ((winEnd-1)-packetBeg)/3 ) - firstBulge;

            // 2x2 reflectors can only occur if a bulge occupies the 3x3 in the
            // bottom-right corner (which begins at index winEnd-3)
            const bool haveSmallBulge =
              ( firstBulge+numFullBulges < numBulges &&
                packetBeg+3*(firstBulge+numFullBulges) == winEnd-3 );

            ComputeReflectors
            ( H, winBeg, shifts, W, packetBeg, firstBulge, numFullBulges,
              haveSmallBulge, ctrl.progress );

            Int transformBeg;
            if( ctrl.accumulateReflections )
                transformBeg = Max( winBeg, ghostCol );
            else if( ctrl.fullTriangle )
                transformBeg = 0;
            else
                transformBeg = winBeg;

            Int transformEnd;
            if( ctrl.accumulateReflections )
                transformEnd = Min( slabEnd, winEnd );
            else if( ctrl.fullTriangle )
                transformEnd = n;
            else
                transformEnd = winEnd;

            ApplyReflectorsOpt
            ( H, winBeg, winEnd,
              slabSize, ghostCol, packetBeg, transformBeg, transformEnd,
              Z, ctrl.wantSchurVecs, U, W,
              firstBulge, numFullBulges, haveSmallBulge,
              ctrl.accumulateReflections, ctrl.progress );
        }

        if( ctrl.accumulateReflections )
        {
            const Int transformBeg = ( ctrl.fullTriangle ? 0 : winBeg );
            const Int transformEnd = ( ctrl.fullTriangle ? n : winEnd );
 
            const Int slabRelBeg = Max(0,(winBeg-1)-ghostCol);
            const Int nU = (slabSize-1) - Max(0,slabEnd-winEnd) - slabRelBeg;

            auto contractInd = IR(0,nU) + slabRelBeg;
            auto UAccum = U( contractInd, contractInd );

            // Horizontal far-from-diagonal application
            const Int rightIndBeg = Min(ghostCol+slabSize,winEnd);
            const Int rightIndEnd = transformEnd;
            const auto rightInd = IR(rightIndBeg,rightIndEnd);
            auto horzInd = IR(0,nU) + (ghostCol+slabRelBeg+1);
            auto HHorzFar = H( horzInd, rightInd );
            Gemm( ADJOINT, NORMAL, F(1), UAccum, HHorzFar, WAccum );
            HHorzFar = WAccum;

            // Vertical far-from-diagonal application
            auto vertInd = IR(transformBeg,Max(winBeg,ghostCol));
            auto HVertFar = H( vertInd, horzInd );
            Gemm( NORMAL, NORMAL, F(1), HVertFar, UAccum, WAccum );
            HVertFar = WAccum;

            if( ctrl.wantSchurVecs )
            {
                auto ZSub = Z( ALL, horzInd );
                Gemm( NORMAL, NORMAL, F(1), ZSub, UAccum, WAccum );
                ZSub = WAccum;
            }
        }
    }
}

template<typename F>
void Sweep
(       Matrix<F>& H,
        Matrix<Complex<Base<F>>>& shifts,
        Matrix<F>& Z,
        Matrix<F>& U,
        Matrix<F>& W,
        Matrix<F>& WAccum,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();

    const Int numShifts = shifts.Height();
    DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts..."); 
    )
    if( numShifts % 2 != 0 )
        LogicError("Expected an even number of shifts");
    if( !IsComplex<F>::value )
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

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_HPP
