/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace hess_qr {

// The best references for the following Aggressive Early Deflation
// implementation are
//
//   Karen Braman, Ralph Byers, and Roy Mathias,
//   "The multishift QR algorithm. Part II: Aggressive Early Deflation",
//   SIAM J. Matrix Anal. Appl., Vol. 23, No. 4, pp. 948--973, 2002
//
// and the LAPACK implementation DLAQR2, which has several distinct differences
// from the suggestions of Braman et al., such as:
//
//   1) Solely using "nearby-diagonal deflation" instead of Braman et al.'s 
//      suggestion of also allowing for "window-Schur deflation".
//
//   2) Using the largest (in magnitude) eigenvalue of a 2x2 Schur block to 
//      determine whether it qualifies for "nearby-diagonal deflation" rather
//      that using the square-root of the absolute value of the determinant
//      (which would correspond to the geometric mean of the eigenvalue
//       magnitudes). 
//
// In both respects, the LAPACK implementation is significantly more
// conservative than the original suggestions of Braman et al.
//

struct AEDInfo
{
    Int numUnconverged=0;
    Int numShiftCandidates=0;
    Int numDeflated=0;
};

template<typename Real>
AEDInfo SpikeDeflation
(       Matrix<Real>& T,
        Matrix<Real>& V,
  const Real& eta,
        Int numUnconverged,
        vector<Real>& work )
{
    DEBUG_CSE

    const Int n = T.Height();
    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    work.resize( n );

    Int winBeg = numUnconverged;
    Int winEnd = n;
    while( winBeg < winEnd )
    {
        const bool twoByTwo =
          ( winEnd==1 ? false : T(winEnd-1,winEnd-2) != zero );

        if( twoByTwo )
        {
            // Follow LAPACK's suit (rather than Braman et al.) and use the 
            // eigenvalue of the 2x2 with largest magnitude in order to 
            // determine if this entry of the spike qualifies for 
            // "nearby-diagonal" deflation. Recall that the 2x2 block is assumed
            // to be in standard form,
            //
            //    | alpha, gamma |, where beta*gamma < 0,
            //    | beta,  alpha |
            //
            // so that the eigenvalues are alpha +- sqrt(beta*gamma) and the
            // spectral radius is |alpha| + sqrt(beta*gamma).
            //
            const Real& alpha = T(winEnd-2,winEnd-2);
            const Real& beta  = T(winEnd-1,winEnd-2); 
            const Real& gamma = T(winEnd-2,winEnd-1);
            const Real spectralRadius =
                Abs(alpha) + Sqrt(Abs(beta))*Sqrt(Abs(gamma));
            const Real scale =
              ( spectralRadius > 0 ? spectralRadius : Abs(eta) );

            // The relevant two entries of spike V^T [eta; zeros(n-1,1)]
            const Real sigma0 = eta*V(0,winEnd-2);
            const Real sigma1 = eta*V(0,winEnd-1);

            if( Max( Abs(sigma0), Abs(sigma1) ) <= Max( smallNum, ulp*scale ) )
            {
                // The two-by-two block satisfies the "nearby-diagonal" test
                winEnd -= 2;
            }
            else
            {
                // Move this undeflatable 2x2 block to the top of the window 
                lapack::SchurExchange
                ( n, &T(0,0), T.LDim(), &V(0,0), V.LDim(),
                  winEnd-2, winBeg, work.data() );
                winBeg += 2;
            }
        }
        else
        {
            const Real spectralRadius = Abs(T(winEnd-1,winEnd-1));
            const Real scale =
              ( spectralRadius > 0 ? spectralRadius : Abs(eta) );

            // The relevant entry of the spike V^T [eta; zeros(n-1,1)]
            const Real sigma = eta*V(0,winEnd-1);

            if( Abs(sigma) <= Max( smallNum, ulp*scale ) )
            {
                // The one-by-one block satisfies the "nearby-diagonal" test
                winEnd -= 1;
            }
            else
            {
                // Move the undeflatable 1x1 block to the top of the window
                lapack::SchurExchange
                ( n, &T(0,0), T.LDim(), &V(0,0), V.LDim(),
                  winEnd-1, winBeg, work.data() );
                winBeg += 1;
            }
        }
    }

    AEDInfo info;
    info.numUnconverged = numUnconverged;
    info.numShiftCandidates = winBeg-numUnconverged;
    info.numDeflated = n-winBeg;
    return info;
}

// TODO(poulson):
// Add the ability to recurse rather than only supporting a direct call to
// WindowedSingle
template<typename Real>
AEDInfo AggressiveEarlyDeflation
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  Int deflationSize,
  Matrix<Real>& wReal,
  Matrix<Real>& wImag,
  bool fullTriangle,
  Matrix<Real>& Z,
  bool wantSchurVecs )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    AEDInfo info;

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
    const bool fullTriangleSub = true;
    const bool wantSchurVecsSub = true;
    const bool demandConvergedSub = false;
    // TODO(poulson): Support recursively calling here AED as well
    Int numUnconverged =
      WindowedSingle
      ( T, 0, blockSize, w1Real, w1Imag,
        fullTriangleSub, V, wantSchurVecsSub, demandConvergedSub );
    DEBUG_ONLY(
      if( numUnconverged != 0 )
          Output
          (numUnconverged," eigenvalues did not converge in WindowedSingle");
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

    info = SpikeDeflation( T, V, spikeValue, numUnconverged, work );
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
            TwoByTwoSchur
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
            Real tau;
            lapack::Reflector( spikeSize, beta, &work[1], 1, tau );
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
        Int applyBeg = ( fullTriangle ? 0 : winBeg );
        auto H01 = H( IR(applyBeg,deflateBeg), deflateInd );
        Gemm( NORMAL, NORMAL, Real(1), H01, V, WAccum );
        H01 = WAccum;

        if( fullTriangle ) 
        {
            // TODO(poulson): Consider forming chunk-by-chunk to save memory
            auto H12 = H( deflateInd, IR(winEnd,END) );
            Gemm( ADJOINT, NORMAL, Real(1), V, H12, WAccum );
            H12 = WAccum;
        }

        if( wantSchurVecs )
        {
            // TODO(poulson): Consider forming chunk-by-chunk to save memory
            auto Z1 = Z(ALL,deflateInd);
            Gemm( NORMAL, NORMAL, Real(1), Z1, V, WAccum );
            Z1 = WAccum;
        }
    }

    return info;
}

} // namespace hess_qr
