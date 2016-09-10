/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_AED_MODIFY_SHIFTS_HPP
#define EL_HESS_SCHUR_AED_MODIFY_SHIFTS_HPP

namespace El {
namespace hess_schur {
namespace aed {

template<typename Real>
Int ModifyShifts
( Int numShiftsRec,
  Int newIterWinSize,
  Int numIterSinceDeflation,
  Int numStaleIterBeforeExceptional, 
  Int winBeg,
  Int winEnd,
  Int shiftBeg,
  const Matrix<Real>& H,
        Matrix<Complex<Real>>& w,
  const HessenbergSchurCtrl& ctrl )
{
    const Real zero(0);
    const Real exceptShift0(Real(4)/Real(3)),
               exceptShift1(-Real(7)/Real(16));

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
    return shiftBeg;
}

template<typename Real>
Int ModifyShifts
( Int numShiftsRec,
  Int newIterWinSize,
  Int numIterSinceDeflation,
  Int numStaleIterBeforeExceptional, 
  Int winBeg,
  Int winEnd,
  Int shiftBeg,
  const Matrix<Complex<Real>>& H,
        Matrix<Complex<Real>>& w,
  const HessenbergSchurCtrl& ctrl )
{
    typedef Complex<Real> F;
    const Real zero(0);
    // For some reason, LAPACK suggests only using a single exceptional shift
    // for complex matrices.
    const Real exceptShift0(Real(4)/Real(3));

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
    return shiftBeg;
}

} // namespace aed
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_AED_MODIFY_SHIFTS_HPP
