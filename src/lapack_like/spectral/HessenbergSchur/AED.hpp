/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_AED_HPP
#define EL_SCHUR_HESS_AED_HPP

#include "./DoubleShift.hpp"

#include "./AED/SpikeDeflation.hpp"
#include "./AED/Nibble.hpp"
#include "./AED/Sweep.hpp"

namespace El {
namespace hess_schur {

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

template<typename Real>
HessenbergSchurInfo
AED
( Matrix<Real>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Real>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE 
    const Int n = H.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const Real zero(0);
    const Real exceptShift0(Real(4)/Real(3)),
               exceptShift1(-Real(7)/Real(16));
    HessenbergSchurInfo info;

    if( n < ctrl.minAEDSize )
    {
        // Run the double-shift QR algorithm
        return DoubleShift( H, w, Z, ctrl );
    }

    w.Resize( n, 1 );

    const Int numShiftsRec = ctrl.numShifts( n, winSize );
    const Int deflationSizeRec = ctrl.deflationSize( n, winSize, numShiftsRec );
    if( ctrl.progress )
    {
        Output
        ("Recommending ",numShiftsRec," shifts and a deflation window of size ",
         deflationSizeRec);
    }
    Int deflationSize = deflationSizeRec;

    // For aed::Sweep
    Matrix<Real> U, W, WAccum;
    auto ctrlSub( ctrl );

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);
    Int decreaseLevel = -1;
    while( winBeg < winEnd )
    {
        if( info.numIterations >= maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("AED QR iteration did not converge");
            else
                break;
        }

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        Int iterBeg=winEnd-1;
        for( ; iterBeg>winBeg; --iterBeg )
            if( H(iterBeg,iterBeg-1) == zero ) 
                break;
        if( ctrl.progress )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
        }
         
        // Intelligently choose a deflation window size
        // --------------------------------------------
        // Cf. LAPACK's DLAQR0 for the high-level approach
        const Int iterWinSize = winEnd-iterBeg;
        if( numIterSinceDeflation < numStaleIterBeforeExceptional )
        {
            // Use the recommendation if possible
            deflationSize = Min( iterWinSize, deflationSizeRec );
        }
        else
        {
            // Double the size if possible
            deflationSize = Min( iterWinSize, 2*deflationSize );
        }
        if( deflationSize >= iterWinSize-1 )
        {
            // Go ahead and increase by at most one to use the full window
            deflationSize = iterWinSize;
        }
        else
        {
            const Int deflationBeg = winEnd - deflationSize;
            if( Abs(H(deflationBeg,  deflationBeg-1)) >
                Abs(H(deflationBeg-1,deflationBeg-2)) )
            {
                ++deflationSize;
            }
        }
        if( numIterSinceDeflation < numStaleIterBeforeExceptional )
        {
            decreaseLevel = -1; 
        }
        else if( decreaseLevel >= 0 || deflationSize == iterWinSize )
        {
            ++decreaseLevel;
            if( deflationSize-decreaseLevel < 2 )
                decreaseLevel = 0;
            deflationSize -= decreaseLevel;
        }

        // Run AED on the bottom-right window of size deflationSize
        ctrlSub.winBeg = iterBeg;
        ctrlSub.winEnd = winEnd;
        auto deflateInfo = aed::Nibble( H, deflationSize, w, Z, ctrlSub );
        const Int numDeflated = deflateInfo.numDeflated;
        winEnd -= numDeflated;
        Int shiftBeg = winEnd - deflateInfo.numShiftCandidates;

        const Int newIterWinSize = winEnd-iterBeg;
        if( numDeflated == 0 ||
          (numDeflated <= ctrl.sufficientDeflation(deflationSize) && 
           newIterWinSize >= ctrl.minAEDSize) )
        {
            Int numShifts = Min( numShiftsRec, Max(2,newIterWinSize-1) );
            numShifts = numShifts - Mod(numShifts,2); 

            if( numIterSinceDeflation > 0 &&
                Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
            {
                // Use exceptional shifts
                shiftBeg = winEnd - numShifts;
                for( Int i=winEnd-1; i>=Max(shiftBeg+1,winBeg+2); i-=2 ) 
                {
                    const Real scale = Abs(H(i,i-1)) + Abs(H(i-1,i-2));
                    Real eta00 = exceptShift0*scale + H(i,i);
                    Real eta01 = scale;
                    Real eta10 = exceptShift1*scale;
                    Real eta11 = eta00;
                    schur::TwoByTwo
                    ( eta00, eta01,
                      eta10, eta11,
                      w(i-1), w(i) );
                }
                if( shiftBeg == winBeg )
                {
                    w(shiftBeg) = w(shiftBeg+1) = H(shiftBeg+1,shiftBeg+1);
                }
            }
            else
            {
                if( winEnd-shiftBeg <= numShifts/2 )
                {
                    // Grab more shifts from another trailing submatrix
                    shiftBeg = winEnd - numShifts;
                    auto shiftsInd = IR(shiftBeg,shiftBeg+numShifts);
                    auto HShifts = H(shiftsInd,shiftsInd);
                    auto wShifts = w(shiftsInd,ALL);
                    auto HShiftsCopy( HShifts );

                    auto ctrlShifts( ctrl );
                    ctrlShifts.winBeg = 0;
                    ctrlShifts.winEnd = numShifts;
                    ctrlShifts.fullTriangle = false;
                    ctrlShifts.demandConverged = false;
                    auto infoShifts =
                      HessenbergSchur( HShiftsCopy, wShifts, ctrlShifts );

                    shiftBeg += infoShifts.numUnconverged;
                    if( shiftBeg >= winEnd-1 )
                    {
                        // This should be very rare; use eigenvalues of 2x2
                        Real eta00 = H(winEnd-2,winEnd-2);
                        Real eta01 = H(winEnd-2,winEnd-1);
                        Real eta10 = H(winEnd-1,winEnd-2);
                        Real eta11 = H(winEnd-1,winEnd-1);
                        schur::TwoByTwo
                        ( eta00, eta01,
                          eta10, eta11,
                          w(winEnd-2), w(winEnd-1) );
                        shiftBeg = winEnd-2;
                    }
                }
                if( winEnd-shiftBeg > numShifts )
                {
                    bool sorted = false;
                    for( Int k=winEnd-1; k>shiftBeg; --k )
                    {
                        if( sorted )
                            break;
                        sorted = true;
                        for( Int i=shiftBeg; i<k; ++i )
                        {
                            if( OneAbs(w(i)) < OneAbs(w(i+1)) )
                            {
                                sorted = false;
                                RowSwap( w, i, i+1 );
                            }
                        }
                    }
                }
                // Pair together the real shifts
                auto wSub = w(IR(shiftBeg,winEnd),ALL); 
                aed::PairShifts( wSub );
            }

            if( winBeg-shiftBeg == 2 )
            {
                // Use a single real shift twice instead of using two separate
                // real shifts; we choose the one closest to the bottom-right
                // entry, as it is our best guess as to the smallest eigenvalue
                if( w(winEnd-1).imag() == zero ) 
                {
                    if( Abs(w(winEnd-1).real()-H(winEnd-1,winEnd-1)) <
                        Abs(w(winEnd-2).real()-H(winEnd-1,winEnd-1)) )
                    {
                        w(winEnd-2) = w(winEnd-1);
                    }
                    else
                    {
                        w(winEnd-1) = w(winEnd-2);
                    }
                }
            }

            // Use the smallest magnitude shifts
            numShifts = Min( numShifts, winEnd-shiftBeg );
            numShifts = numShifts - Mod(numShifts,2);
            shiftBeg = winEnd - numShifts;

            // Perform a small-bulge sweep
            auto wSub = w(IR(shiftBeg,winEnd),ALL); 
            ctrlSub.winBeg = iterBeg;
            ctrlSub.winEnd = winEnd;
            aed::Sweep( H, wSub, Z, U, W, WAccum, ctrlSub );
        }
        else if( ctrl.progress )
            Output("  Skipping QR sweep");

        ++info.numIterations;
        if( numDeflated > 0 )
            numIterSinceDeflation = 0;
        else
            ++numIterSinceDeflation;
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

template<typename Real>
HessenbergSchurInfo
AED
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE 
    typedef Complex<Real> F;

    const Int n = H.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const Real zero(0);
    // For some reason, LAPACK suggests only using a single exceptional shift
    // for complex matrices.
    const Real exceptShift0(Real(4)/Real(3));
    HessenbergSchurInfo info;

    if( n < ctrl.minAEDSize )
    {
        // Run the single-shift QR algorithm
        return SingleShift( H, w, Z, ctrl );
    }

    w.Resize( n, 1 );

    const Int numShiftsRec = ctrl.numShifts( n, winSize );
    const Int deflationSizeRec = ctrl.deflationSize( n, winSize, numShiftsRec );
    if( ctrl.progress )
    {
        Output
        ("Recommending ",numShiftsRec," shifts and a deflation window of size ",
         deflationSizeRec);
    }
    Int deflationSize = deflationSizeRec;

    // For aed::Sweep
    Matrix<F> U, W, WAccum;
    auto ctrlSub( ctrl );

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);
    Int decreaseLevel = -1;
    while( winBeg < winEnd )
    {
        if( info.numIterations >= maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("AED QR iteration did not converge");
            else
                break;
        }

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        Int iterBeg=winEnd-1;
        for( ; iterBeg>winBeg; --iterBeg )
            if( H(iterBeg,iterBeg-1) == zero ) 
                break;
        if( ctrl.progress )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
        }
         
        // Intelligently choose a deflation window size
        // --------------------------------------------
        // Cf. LAPACK's DLAQR0 for the high-level approach
        const Int iterWinSize = winEnd-iterBeg;
        if( numIterSinceDeflation < numStaleIterBeforeExceptional )
        {
            // Use the recommendation if possible
            deflationSize = Min( iterWinSize, deflationSizeRec );
        }
        else
        {
            // Double the size if possible
            deflationSize = Min( iterWinSize, 2*deflationSize );
        }
        if( deflationSize >= iterWinSize-1 )
        {
            // Go ahead and increase by at most one to use the full window
            deflationSize = iterWinSize;
        }
        else
        {
            const Int deflationBeg = winEnd - deflationSize;
            if( OneAbs(H(deflationBeg,  deflationBeg-1)) >
                OneAbs(H(deflationBeg-1,deflationBeg-2)) )
            {
                ++deflationSize;
            }
        }
        if( numIterSinceDeflation < numStaleIterBeforeExceptional )
        {
            decreaseLevel = -1; 
        }
        else if( decreaseLevel >= 0 || deflationSize == iterWinSize )
        {
            ++decreaseLevel;
            if( deflationSize-decreaseLevel < 2 )
                decreaseLevel = 0;
            deflationSize -= decreaseLevel;
        }

        // Run AED on the bottom-right window of size deflationSize
        ctrlSub.winBeg = iterBeg;
        ctrlSub.winEnd = winEnd;
        auto deflateInfo = aed::Nibble( H, deflationSize, w, Z, ctrlSub );
        const Int numDeflated = deflateInfo.numDeflated;
        winEnd -= numDeflated;
        Int shiftBeg = winEnd - deflateInfo.numShiftCandidates;

        const Int newIterWinSize = winEnd-iterBeg;
        if( numDeflated == 0 ||
          (numDeflated <= ctrl.sufficientDeflation(deflationSize) && 
           newIterWinSize >= ctrl.minAEDSize) )
        {
            Int numShifts = Min( numShiftsRec, Max(2,newIterWinSize-1) );
            numShifts = numShifts - Mod(numShifts,2); 

            if( numIterSinceDeflation > 0 &&
                Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
            {
                // Use exceptional shifts
                shiftBeg = winEnd - numShifts;
                for( Int i=winEnd-1; i>=shiftBeg+1; i-=2 )
                {
                    w(i-1) = w(i) = H(i,i) + exceptShift0*OneAbs(H(i,i-1));
                }
            }
            else
            {
                if( winEnd-shiftBeg <= numShifts/2 )
                {
                    // Grab more shifts from another trailing submatrix
                    shiftBeg = winEnd - numShifts;
                    auto shiftsInd = IR(shiftBeg,shiftBeg+numShifts);
                    auto HShifts = H(shiftsInd,shiftsInd);
                    auto wShifts = w(shiftsInd,ALL);
                    auto HShiftsCopy( HShifts );

                    auto ctrlShifts( ctrl );
                    ctrlShifts.winBeg = 0;
                    ctrlShifts.winEnd = numShifts;
                    ctrlShifts.fullTriangle = false;
                    ctrlShifts.demandConverged = false;
                    auto infoShifts =
                      HessenbergSchur( HShiftsCopy, wShifts, ctrlShifts );

                    shiftBeg += infoShifts.numUnconverged;
                    if( shiftBeg >= winEnd-1 )
                    {
                        // This should be very rare; use eigenvalues of 2x2
                        F eta00 = H(winEnd-2,winEnd-2);
                        F eta01 = H(winEnd-2,winEnd-1);
                        F eta10 = H(winEnd-1,winEnd-2);
                        F eta11 = H(winEnd-1,winEnd-1);
                        schur::TwoByTwo
                        ( eta00, eta01,
                          eta10, eta11,
                          w(winEnd-2), w(winEnd-1) );
                        shiftBeg = winEnd-2;
                    }
                }
                if( winEnd-shiftBeg > numShifts )
                {
                    bool sorted = false;
                    for( Int k=winEnd-1; k>shiftBeg; --k )
                    {
                        if( sorted )
                            break;
                        sorted = true;
                        for( Int i=shiftBeg; i<k; ++i )
                        {
                            if( OneAbs(w(i)) < OneAbs(w(i+1)) )
                            {
                                sorted = false;
                                RowSwap( w, i, i+1 );
                            }
                        }
                    }
                }
            }

            if( winBeg-shiftBeg == 2 )
            {
                // Use a single real shift twice instead of using two separate
                // real shifts; we choose the one closest to the bottom-right
                // entry, as it is our best guess as to the smallest eigenvalue
                if( w(winEnd-1).imag() == zero ) 
                {
                    if( Abs(w(winEnd-1).real()-H(winEnd-1,winEnd-1)) <
                        Abs(w(winEnd-2).real()-H(winEnd-1,winEnd-1)) )
                    {
                        w(winEnd-2) = w(winEnd-1);
                    }
                    else
                    {
                        w(winEnd-1) = w(winEnd-2);
                    }
                }
            }

            // Use the smallest magnitude shifts
            numShifts = Min( numShifts, winEnd-shiftBeg );
            numShifts = numShifts - Mod(numShifts,2);
            shiftBeg = winEnd - numShifts;

            // Perform a small-bulge sweep
            auto wSub = w(IR(shiftBeg,winEnd),ALL); 
            ctrlSub.winBeg = iterBeg;
            ctrlSub.winEnd = winEnd;
            aed::Sweep( H, wSub, Z, U, W, WAccum, ctrlSub );
        }
        else if( ctrl.progress )
            Output("  Skipping QR sweep");

        ++info.numIterations;
        if( numDeflated > 0 )
            numIterSinceDeflation = 0;
        else
            ++numIterSinceDeflation;
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_AED_HPP
