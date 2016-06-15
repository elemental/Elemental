/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace hess_qr {

// The best references are
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
template<typename Real>
Int SpikeDeflation
(       Matrix<Real>& T,
        Matrix<Real>& V,
  const Real& eta,
        Int numUnconverged,
        vector<Real>& work )
{
    DEBUG_ONLY(CSE cse("hess_qr::SpikeDeflation"))

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
    // Return the number of deflated eigenvalues
    return n-winBeg;
}

template<typename Real>
void AggressiveEarlyDeflation
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
    const Real zero(0);
    if( winBeg > winEnd )
        return;
    if( deflationSize < 1 )
        return;

    Int blocksize = Min( deflationSize, winEnd-winBeg );
    const Int deflateBeg = winEnd-blocksize;

    auto deflateInd = IR(deflateBeg,winEnd);
    auto HBR = H( deflateInd, deflateInd );

    const Int n = H.Height();

    if( blocksize == 1 )
    {
        // TODO: Take a shortcut
    }

    // TODO: Only copy the upper-Hessenberg portion of HBR
    Matrix<Real> TBR( HBR );
    auto wBReal = wReal( IR(deflateBeg,winEnd), ALL );
    auto wBImag = wImag( IR(deflateBeg,winEnd), ALL );
    Matrix<Real> VBR;
    Identity( VBR, blocksize, blocksize );
    Int numUnconverged =
      WindowedSingle
      ( TBR, 0, blocksize, wBReal, wBImag, true, VBR, true, false );

    // Clear the two diagonals below the upper-Hessenberg portion for
    // SchurExchange
    for( Int i=0; i<blocksize-3; ++i )
    {
        TBR(i+2,i) = zero;
        TBR(i+3,i) = zero;
    }
    if( blocksize > 2 )
        TBR(blocksize-1,blocksize-3) = zero;

    vector<Real> work(blocksize);

    Real etaBR = ( deflateBeg == winBeg ? zero : H(deflateBeg,deflateBeg-1) );
    Int numDeflated = SpikeDeflation( TBR, VBR, etaBR, numUnconverged, work );
    if( numDeflated == blocksize )
    {
        // The entire spike has deflated
        etaBR = zero;
    }
    if( numDeflated > 0 )
    {
        // TODO: Sort diagonal blocks? 
    }

    // Reform eigenvalues
    for( Int i=blocksize; i>numUnconverged; )
    {
        if( i == numUnconverged+1 || TBR(i,i-1) == zero )
        {
            // 1x1 block
            wReal(deflateBeg+i) = TBR(i,i);
            wImag(deflateBeg+i) = zero;
            i -= 1;
        }
        else
        {
            // 2x2 block
            Real alpha00 = TBR(i-1,i-1);
            Real alpha10 = TBR(i,  i-1);
            Real alpha01 = TBR(i-1,i  );
            Real alpha11 = TBR(i,  i  );
            Real c, s;
            TwoByTwoSchur
            ( alpha00, alpha01,
              alpha10, alpha11, c, s,
              wReal(deflateBeg+i-1), wImag(deflateBeg+i-1),
              wReal(deflateBeg+i  ), wImag(deflateBeg+i  ) );
            i -= 2;
        }
    }

    const Int spikeSize = blocksize - numDeflated;
    if( spikeSize < blocksize || etaBR == zero )
    {
        // Either we deflated at least one eigenvalue or we can simply
        // rotate the deflation window into Schur form

        if( spikeSize > 1 && etaBR != zero )
        {
            // The spike needs to be reduced to length one while maintaining
            // the Hessenberg form of the deflation window
            for( Int i=0; i<spikeSize; ++i )
                work[i] = VBR(0,i);
            // TODO
        }
    }
}

} // namespace hess_qr
