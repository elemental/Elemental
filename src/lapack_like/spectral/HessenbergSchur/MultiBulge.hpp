/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_HPP
#define EL_HESS_SCHUR_MULTIBULGE_HPP

#include "./MultiBulge/Sweep.hpp"

namespace El {
namespace hess_schur {

template<typename Real>
HessenbergSchurInfo
MultiBulge
( Matrix<Real>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Real>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE 
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const Real zero(0);
    const Real exceptShift0(Real(4)/Real(3)),
               exceptShift1(-Real(7)/Real(16));
    HessenbergSchurInfo info;

    if( n < ctrl.minMultiBulgeSize )
    {
        // Run the double-shift QR algorithm
        return DoubleShift( H, w, Z, ctrl );
    }

    w.Resize( n, 1 );

    const Int numShiftsRec = ctrl.numShifts( n, winSize );
    if( ctrl.progress )
        Output("Recommending ",numShiftsRec," shifts");

    Matrix<Real> U, W, WAccum, shifts;
    auto ctrlSweep( ctrl ), ctrlShifts( ctrl );
    ctrlShifts.fullTriangle = false;

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);
    Int iterBegLast=-1, winEndLast=-1;
    while( winBeg < winEnd )
    {
        if( info.numIterations >= maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("MultiBulge QR iteration did not converge");
            else
                break;
        }

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        Int iterBeg = winBeg;
        {
            auto winInd = IR(iterBeg,winEnd);
            iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
        }
        if( iterBeg > winBeg )
        {
            H(iterBeg,iterBeg-1) = zero;
        }
        if( iterBeg == winEnd-1 )
        {
            w(iterBeg) = H(iterBeg,iterBeg);
            --winEnd;
            numIterSinceDeflation = 0;
            continue;
        }
        else if( iterBeg == winEnd-2 )
        {
            Real c, s;
            schur::TwoByTwo
            ( H(winEnd-2,winEnd-2), H(winEnd-2,winEnd-1),
              H(winEnd-1,winEnd-2), H(winEnd-1,winEnd-1),
              w(iterBeg), w(iterBeg+1),
              c, s );
            if( ctrl.fullTriangle )
            {
                if( n > winEnd )
                    blas::Rot
                    ( n-winEnd,
                      &H(winEnd-2,winEnd), H.LDim(),
                      &H(winEnd-1,winEnd), H.LDim(),
                      c, s );
                blas::Rot
                ( winEnd-2, &H(0,winEnd-2), 1, &H(0,winEnd-1), 1, c, s );
            }
            if( ctrl.wantSchurVecs )
                blas::Rot( nZ, &Z(0,winEnd-2), 1, &Z(0,winEnd-1), 1, c, s );
            winEnd -= 2;
            numIterSinceDeflation = 0;
            continue;
        }

        const Int iterWinSize = winEnd-iterBeg;
        if( ctrl.progress )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
        }
        const Int shiftBeg = Max(iterBeg,winEnd-numShiftsRec);
        const Int numShifts = winEnd - shiftBeg;
        auto shiftInd = IR(shiftBeg,winEnd);
        auto wShifts = w(shiftInd,ALL); 
         
        if( numIterSinceDeflation > 0 &&
            Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
        {
            // Use exceptional shifts
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
            // Compute the eigenvalues of the bottom-right window
            auto HShifts = H(shiftInd,shiftInd);
            auto shiftInfo = HessenbergSchur( HShifts, wShifts, ctrlShifts );
            multibulge::PairShifts( wShifts );
        }

        if( winBeg-shiftBeg == 2 )
        {
            // Use a single real shift twice instead of using two separate
            // real shifts; we choose the one closest to the bottom-right
            // entry, as it is our best guess as to the smallest eigenvalue
            if( wShifts(numShifts-1).imag() == zero ) 
            {
                if( Abs(wShifts(numShifts-1).real()-H(winEnd-1,winEnd-1)) <
                    Abs(wShifts(numShifts-2).real()-H(winEnd-1,winEnd-1)) )
                {
                    wShifts(numShifts-2) = wShifts(numShifts-1);
                }
                else
                {
                    wShifts(numShifts-1) = wShifts(numShifts-2);
                }
            }
        }

        // Perform a small-bulge sweep
        ctrlSweep.winBeg = iterBeg;
        ctrlSweep.winEnd = winEnd;
        multibulge::Sweep( H, wShifts, Z, U, W, WAccum, ctrlSweep );

        ++info.numIterations;
        if( iterBeg == iterBegLast && winEnd == winEndLast )
            ++numIterSinceDeflation;
        iterBegLast = iterBeg;
        winEndLast = winEnd;
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

template<typename Real>
HessenbergSchurInfo
MultiBulge
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE 
    typedef Complex<Real> F;

    const Int n = H.Height();
    const Int nZ = Z.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const Real zero(0);
    // For some reason, LAPACK suggests only using a single exceptional shift
    // for complex matrices.
    const Real exceptShift0(Real(4)/Real(3));
    HessenbergSchurInfo info;

    if( n < ctrl.minMultiBulgeSize )
    {
        // Run the single-shift QR algorithm
        return SingleShift( H, w, Z, ctrl );
    }

    w.Resize( n, 1 );

    const Int numShiftsRec = ctrl.numShifts( n, winSize );
    if( ctrl.progress )
        Output("Recommending ",numShiftsRec," shifts");

    Matrix<F> U, W, WAccum;
    auto ctrlSweep( ctrl ), ctrlShifts( ctrl );
    ctrlShifts.fullTriangle = false;

    Int numIterSinceDeflation = 0;
    const Int numStaleIterBeforeExceptional = 5;
    // Cf. LAPACK's DLAQR0 for this choice
    const Int maxIter =
      Max(30,2*numStaleIterBeforeExceptional) * Max(10,winSize);
    Int iterBegLast=-1, winEndLast=-1;
    while( winBeg < winEnd )
    {
        if( info.numIterations >= maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("MultiBulge QR iteration did not converge");
            else
                break;
        }

        // Detect an irreducible Hessenberg window, [iterBeg,winEnd)
        // ---------------------------------------------------------
        Int iterBeg = winBeg;
        {
            auto winInd = IR(iterBeg,winEnd);
            iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
        }
        if( iterBeg > winBeg )
        {
            H(iterBeg,iterBeg-1) = zero;
        }
        if( iterBeg == winEnd-1 )
        {
            w(iterBeg) = H(iterBeg,iterBeg);
            --winEnd;
            numIterSinceDeflation = 0;
            continue;
        }
        if( ctrl.progress )
        {
            Output("Iter. ",info.numIterations,": ");
            Output("  window is [",iterBeg,",",winEnd,")");
        }
        const Int shiftBeg = Max(iterBeg,winEnd-numShiftsRec);
        const Int numShifts = winEnd - shiftBeg;
        auto shiftInd = IR(shiftBeg,winEnd);
        auto wShifts = w(shiftInd,ALL);

        if( numIterSinceDeflation > 0 &&
            Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
        {
            for( Int i=winEnd-1; i>=shiftBeg+1; i-=2 )
                w(i-1) = w(i) = H(i,i) + exceptShift0*OneAbs(H(i,i-1));
        }
        else
        {
            // Compute the eigenvalues of the bottom-right window
            auto HShifts = H(shiftInd,shiftInd);
            auto shiftInfo = HessenbergSchur( HShifts, wShifts, ctrlShifts );
            multibulge::PairShifts( wShifts );
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

        // Perform a small-bulge sweep
        ctrlSweep.winBeg = iterBeg;
        ctrlSweep.winEnd = winEnd;
        multibulge::Sweep( H, wShifts, Z, U, W, WAccum, ctrlSweep );

        ++info.numIterations;
        if( iterBeg == iterBegLast && winEnd == winEndLast )
            ++numIterSinceDeflation;
        iterBegLast = iterBeg;
        winEndLast = winEnd;
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_HPP
