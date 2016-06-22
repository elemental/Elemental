/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESSQR_AED_NIBBLE_HPP
#define EL_SCHUR_HESSQR_AED_NIBBLE_HPP

#include "./SpikeDeflation.hpp"

namespace El {
namespace schur {
namespace hess_qr {
namespace aed {

template<typename Real>
AEDInfo Nibble
( Matrix<Real>& H,
  Int deflationSize,
  Matrix<Real>& wReal,
  Matrix<Real>& wImag,
  Matrix<Real>& Z,
  const HessenbergQRCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    AEDInfo info;

    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    if( winBeg > winEnd )
        return info;
    if( deflationSize < 1 )
        return info;

    Int blockSize = Min( deflationSize, winEnd-winBeg );
    const Int deflateBeg = winEnd-blockSize;

    // If the deflation window touches the beginning of the full window,
    // then there is no spike
    Real spikeValue =
      ( deflateBeg==winBeg ? zero : H(deflateBeg,deflateBeg-1) );

    if( blockSize == 1 )
    {
        wReal(deflateBeg) = H(deflateBeg,deflateBeg);
        wImag(deflateBeg) = zero;         
        if( Abs(spikeValue) <= Max( smallNum, ulp*Abs(wReal(deflateBeg)) ) )
        {
            // The offdiagonal entry was small enough to deflate
            info.numDeflated = 1;
            if( deflateBeg > winBeg )
            {
                // Explicitly deflate by zeroing the offdiagonal entry
                H(deflateBeg,deflateBeg-1) = zero;
            }
        }
        else
        {
            // The offdiagonal entry was too large to deflate
            info.numShiftCandidates = 1;
        }
        return info;
    }

    auto deflateInd = IR(deflateBeg,winEnd);
    auto H11 = H( deflateInd, deflateInd );

    // NOTE(poulson): We could only copy the upper-Hessenberg portion of H11
    Matrix<Real> T( H11 ); // TODO(poulson): Reuse this matrix?
    auto w1Real = wReal( deflateInd, ALL );
    auto w1Imag = wImag( deflateInd, ALL );
    Matrix<Real> V;
    Identity( V, blockSize, blockSize );
    auto ctrlSub( ctrl );
    ctrlSub.winBeg = 0;
    ctrlSub.winEnd = blockSize;
    ctrlSub.fullTriangle = true;
    ctrlSub.wantSchurVecs = true;
    ctrlSub.demandConverged = false;
    ctrlSub.useAED = ( ctrl.recursiveAED ? true : false );
    auto infoSub = HessenbergQR( T, w1Real, w1Imag, V, ctrlSub );
    DEBUG_ONLY(
      if( infoSub.numUnconverged != 0 )
          Output(infoSub.numUnconverged," eigenvalues did not converge");
    )

    // Clear the two diagonals below the upper-Hessenberg portion for
    // SchurExchange
    for( Int i=0; i<blockSize-3; ++i )
    {
        T(i+2,i) = zero;
        T(i+3,i) = zero;
    }
    if( blockSize >= 3 )
        T(blockSize-1,blockSize-3) = zero;

    vector<Real> work(2*blockSize);

    info = SpikeDeflation( T, V, spikeValue, infoSub.numUnconverged, work );
    if( ctrl.progress )
    {
        Output("AED info:");
        Output("  numUnconverged=",info.numUnconverged);
        Output("  numShiftCandidates=",info.numShiftCandidates);
        Output("  numDeflated=",info.numDeflated);
    }
    if( info.numUnconverged+info.numShiftCandidates == 0 )
    {
        // The entire spike has deflated
        spikeValue = zero;
    }
    if( info.numDeflated > 0 )
    {
        // TODO(poulson):
        //
        // Follow LAPACK's lead and sort the diagonal blocks of T from largest
        // to smallest magnitude to improve the accuracy of the algorithm on
        // graded matrices.
        //
        // Perhaps this should be its own subroutine that accepts a SortType,
        // which can take the values {ASCENDING, DESCENDING, UNSORTED}.
    }

    // Reform the eigenvalues and shift candidates by looping over the converged
    // eigenvalues from last to first
    for( Int i=blockSize-1; i>=info.numUnconverged; )
    {
        if( i == info.numUnconverged || T(i,i-1) == zero )
        {
            // 1x1 block
            wReal(deflateBeg+i) = T(i,i);
            wImag(deflateBeg+i) = zero;
            i -= 1;
        }
        else
        {
            // 2x2 block
            Real alpha00 = T(i-1,i-1);
            Real alpha10 = T(i,  i-1);
            Real alpha01 = T(i-1,i  );
            Real alpha11 = T(i,  i  );
            Real c, s;
            lapack::TwoByTwoSchur
            ( alpha00, alpha01,
              alpha10, alpha11, c, s,
              wReal(deflateBeg+i-1), wImag(deflateBeg+i-1),
              wReal(deflateBeg+i  ), wImag(deflateBeg+i  ) );
            i -= 2;
        }
    }

    const Int spikeSize = info.numUnconverged + info.numShiftCandidates;
    if( spikeSize < blockSize || spikeValue == zero )
    {
        // Either we deflated at least one eigenvalue or we can simply
        // rotate the deflation window into Schur form
        auto spikeInd = IR(0,spikeSize);
        auto TTL = T(spikeInd,spikeInd);
        auto TTR = T(spikeInd,IR(spikeSize,END));
        auto VL = V(ALL,spikeInd);
        Matrix<Real> phaseT;

        if( spikeSize > 1 && spikeValue != zero )
        {
            // The spike needs to be reduced to length one while maintaining
            // the Hessenberg form of the deflation window
            for( Int i=0; i<spikeSize; ++i )
                work[i] = V(0,i);

            // Compute a Householder reflector for condensing the spike
            Real beta = work[0];
            Real tau = lapack::Reflector( spikeSize, beta, &work[1], 1 );
            work[0] = Real(1);

            // Force T to be upper Hessenberg
            MakeTrapezoidal( UPPER, T, -1 );

            lapack::ApplyReflector
            ( true, spikeSize, blockSize,
              &work[0], 1, tau,
              T.Buffer(), T.LDim(),
              &work[blockSize] );

            lapack::ApplyReflector
            ( false, spikeSize, spikeSize,
              &work[0], 1, tau,
              TTL.Buffer(), TTL.LDim(),
              &work[blockSize] );

            lapack::ApplyReflector
            ( false, blockSize, spikeSize,
              &work[0], 1, tau,
              VL.Buffer(), VL.LDim(),
              &work[blockSize] );

            Hessenberg( UPPER, TTL, phaseT );
            hessenberg::ApplyQ( LEFT, UPPER, ADJOINT, TTL, phaseT, TTR );
        }

        if( deflateBeg >= 1 )
            H(deflateBeg,deflateBeg-1) = spikeValue*V(0,0);
        // NOTE(poulson): We could copy only the upper Hessenberg part of T
        H11 = T;
        MakeTrapezoidal( UPPER, H11, -1 );

        if( spikeSize > 1 && spikeValue != zero )
        {
            hessenberg::ApplyQ( RIGHT, UPPER, NORMAL, TTL, phaseT, VL );
        }

        Matrix<Real> WAccum; // TODO(poulson): Reuse this buffer?

        // TODO(poulson): Consider forming chunk-by-chunk to save memory
        Int applyBeg = ( ctrl.fullTriangle ? 0 : winBeg );
        auto H01 = H( IR(applyBeg,deflateBeg), deflateInd );
        Gemm( NORMAL, NORMAL, Real(1), H01, V, WAccum );
        H01 = WAccum;

        if( ctrl.fullTriangle ) 
        {
            // TODO(poulson): Consider forming chunk-by-chunk to save memory
            auto H12 = H( deflateInd, IR(winEnd,END) );
            Gemm( ADJOINT, NORMAL, Real(1), V, H12, WAccum );
            H12 = WAccum;
        }

        if( ctrl.wantSchurVecs )
        {
            // TODO(poulson): Consider forming chunk-by-chunk to save memory
            auto Z1 = Z(ALL,deflateInd);
            Gemm( NORMAL, NORMAL, Real(1), Z1, V, WAccum );
            Z1 = WAccum;
        }
    }

    return info;
}

} // namespace aed
} // namespace hess_qr
} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_HESSQR_AED_NIBBLE_HPP
