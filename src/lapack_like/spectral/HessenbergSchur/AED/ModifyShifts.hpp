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
    EL_DEBUG_CSE
    // Compute the ideal number of (even) shifts to apply this sweep
    Int numShiftsIdeal = Min( numShiftsRec, Max(2,newIterWinSize-1) );
    numShiftsIdeal = numShiftsIdeal - Mod(numShiftsIdeal,2);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts to attempt to break out of a convergence rut
        if( ctrl.progress )
            Output("Using exceptional shifts");
        shiftBeg = winEnd - numShiftsIdeal;
        const Real exceptShift0(Real(4)/Real(3)),
                   exceptShift1(-Real(7)/Real(16));
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
            w(shiftBeg) = w(shiftBeg+1) = H(shiftBeg+1,shiftBeg+1);
    }
    else
    {
        // TODO(poulson): Come up with a principled switching point
        if( winEnd-shiftBeg <= numShiftsIdeal/2 )
        {
            // Grab more shifts from another trailing submatrix
            if( ctrl.progress )
                Output("Grabbing more shifts");
            shiftBeg = winEnd - numShiftsIdeal;
            auto shiftsInd = IR(shiftBeg,shiftBeg+numShiftsIdeal);
            auto HShifts = H(shiftsInd,shiftsInd);
            auto wShifts = w(shiftsInd,ALL);
            auto HShiftsCopy( HShifts );

            auto ctrlShifts( ctrl );
            ctrlShifts.winBeg = 0;
            ctrlShifts.winEnd = numShiftsIdeal;
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
    }
    if( winEnd-shiftBeg == 2 )
    {
        // Use a single real shift twice instead of using two separate
        // real shifts; we choose the one closest to the bottom-right
        // entry, as it is our best guess as to the smallest eigenvalue
        if( w(winEnd-1).imag() == Real(0) )
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

    auto wShifts = w(IR(shiftBeg,winEnd),ALL);
    multibulge::PairShifts( wShifts );
    if( winEnd-shiftBeg > numShiftsIdeal )
    {
        // We have too many shifts and need to eliminate some
        if( ctrl.sortShifts )
        {
            // Sort the shifts into descending order in a manner which preserves
            // consecutive pairs
            if( ctrl.progress )
                Output("Sorting shifts");
            typedef std::pair<Complex<Real>,Complex<Real>> ShiftPair;
            const Int numShifts = winEnd - shiftBeg;
            const Int remainder = numShifts % 2;
            const Int numShiftPairs = (numShifts/2) + remainder;
            std::vector<ShiftPair> shiftPairs(numShiftPairs);
            if( remainder == 1 )
            {
                shiftPairs[0].first = w(shiftBeg);
                shiftPairs[0].second = Conj(w(shiftBeg));
            }
            for( Int i=0; i<numShiftPairs-remainder; ++i )
            {
                shiftPairs[i+remainder].first = w(shiftBeg+2*i+remainder);
                shiftPairs[i+remainder].second = w(shiftBeg+2*i+1+remainder);
            }
            std::sort
            ( shiftPairs.begin(), shiftPairs.end(),
              [](const ShiftPair& a, const ShiftPair& b)
              { return OneAbs(b.first) < OneAbs(a.first); } );
            const Int numKeptPairs = numShiftsIdeal / 2;
            const Int numDroppedPairs = numShiftPairs - numKeptPairs;
            const Int shiftBegNew = winEnd - numShiftsIdeal;
            for( Int i=0; i<numKeptPairs; ++i )
            {
                w(shiftBegNew+2*i) = shiftPairs[i+numDroppedPairs].first;
                w(shiftBegNew+2*i+1) = shiftPairs[i+numDroppedPairs].second;
            }
        }
        shiftBeg = winEnd - numShiftsIdeal;
    }
    else
    {
        Int numShifts = winEnd-shiftBeg;
        numShifts = numShifts - Mod(numShifts,2);
        shiftBeg = winEnd - numShifts;
    }
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
  const DistMatrix<Real,MC,MR,BLOCK>& H,
  const Matrix<Real>& hMainWin,
  const Matrix<Real>& hSubWin,
        DistMatrix<Complex<Real>,STAR,STAR>& w,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = H.Grid();
    auto& wLoc = w.Matrix();

    // Compute the ideal number of (even) shifts to apply this sweep
    Int numShiftsIdeal = Min( numShiftsRec, Max(2,newIterWinSize-1) );
    numShiftsIdeal = numShiftsIdeal - Mod(numShiftsIdeal,2);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts to attempt to break out of a convergence rut
        if( ctrl.progress && grid.Rank() == 0 )
            Output("Using exceptional shifts");
        shiftBeg = winEnd - numShiftsIdeal;
        const Real exceptShift0(Real(4)/Real(3)),
                   exceptShift1(-Real(7)/Real(16));
        for( Int i=winEnd-1; i>=Max(shiftBeg+1,winBeg+2); i-=2 )
        {
            const Real scale = Abs(hSubWin(i-1-winBeg)) +
              Abs(hSubWin(i-2-winBeg));
            Real eta00 = exceptShift0*scale + hMainWin(i-winBeg);
            Real eta01 = scale;
            Real eta10 = exceptShift1*scale;
            Real eta11 = eta00;
            if( grid.VCRank() == 0 )
            {
                schur::TwoByTwo
                ( eta00, eta01,
                  eta10, eta11,
                  wLoc(i-1), wLoc(i) );
            }
        }
        if( shiftBeg == winBeg )
            wLoc(shiftBeg) = wLoc(shiftBeg+1) = hMainWin(shiftBeg+1-winBeg);
        auto wLocSub = wLoc( IR(shiftBeg,winEnd), ALL );
        El::Broadcast( wLocSub, grid.VCComm(), 0 );
    }
    else
    {
        // TODO(poulson): Come up with a principled switching point
        if( winEnd-shiftBeg <= numShiftsIdeal/2 )
        {
            // Grab more shifts from another trailing submatrix
            if( ctrl.progress && grid.Rank() == 0 )
                Output("Grabbing more shifts");

            auto ctrlShifts( ctrl );
            ctrlShifts.winBeg = 0;
            ctrlShifts.winEnd = numShiftsIdeal;
            ctrlShifts.fullTriangle = false;
            ctrlShifts.demandConverged = false;

            shiftBeg = winEnd - numShiftsIdeal;
            auto shiftsInd = IR(shiftBeg,shiftBeg+numShiftsIdeal);
            auto HShifts = H(shiftsInd,shiftsInd);
            auto wShifts = w(shiftsInd,ALL);
            const Int numUnconverged =
              multibulge::ConsistentlyComputeEigenvalues
              ( HShifts, wShifts, ctrlShifts );
            shiftBeg += numUnconverged;
            if( grid.Rank() == 0 )
                Output("shiftBeg=",shiftBeg,", winEnd=",winEnd);
            if( shiftBeg >= winEnd-1 )
            {
                // This should be very rare; use eigenvalues of 2x2
                Real eta00 = hMainWin(winEnd-2-winBeg);
                Real eta10 = hSubWin(winEnd-2-winBeg);
                Real eta01 = H.Get(winEnd-2,winEnd-1); // This involves a bcast
                Real eta11 = hMainWin(winEnd-1-winBeg);
                if( grid.VCRank() == 0 )
                {
                    schur::TwoByTwo
                    ( eta00, eta01,
                      eta10, eta11,
                      wLoc(winEnd-2), wLoc(winEnd-1) );
                }
                auto wSub = wLoc( IR(winEnd-2,winEnd), ALL );
                El::Broadcast( wSub, grid.VCComm(), 0 );
                shiftBeg = winEnd-2;
            }
        }
    }
    if( winEnd-shiftBeg == 2 )
    {
        if( ctrl.progress && grid.Rank() == 0 )
            Output("Forcing the same shift (if real)");
        // Use a single real shift twice instead of using two separate
        // real shifts; we choose the one closest to the bottom-right
        // entry, as it is our best guess as to the smallest eigenvalue
        if( wLoc(winEnd-1).imag() == Real(0) )
        {
            if( Abs(wLoc(winEnd-1).real()-hMainWin(winEnd-1-winBeg)) <
                Abs(wLoc(winEnd-2).real()-hMainWin(winEnd-1-winBeg)) )
            {
                wLoc(winEnd-2) = wLoc(winEnd-1);
            }
            else
            {
                wLoc(winEnd-1) = wLoc(winEnd-2);
            }
        }
    }

    auto wLocShifts = wLoc(IR(shiftBeg,winEnd),ALL);
    multibulge::PairShifts( wLocShifts );
    if( winEnd-shiftBeg > numShiftsIdeal )
    {
        // We have too many shifts and need to eliminate some
        if( ctrl.sortShifts )
        {
            // Sort the shifts into descending order in a manner which preserves
            // consecutive pairs
            if( ctrl.progress && grid.Rank() == 0 )
                Output("Sorting shifts");
            typedef std::pair<Complex<Real>,Complex<Real>> ShiftPair;
            const Int numShifts = winEnd - shiftBeg;
            const Int remainder = numShifts % 2;
            const Int numShiftPairs = (numShifts/2) + remainder;
            std::vector<ShiftPair> shiftPairs(numShiftPairs);
            if( remainder == 1 )
            {
                shiftPairs[0].first = wLoc(shiftBeg);
                shiftPairs[0].second = Conj(wLoc(shiftBeg));
            }
            for( Int i=0; i<numShiftPairs-remainder; ++i )
            {
                shiftPairs[i+remainder].first = wLoc(shiftBeg+2*i+remainder);
                shiftPairs[i+remainder].second = wLoc(shiftBeg+2*i+1+remainder);
            }
            std::sort
            ( shiftPairs.begin(), shiftPairs.end(),
              [](const ShiftPair& a, const ShiftPair& b)
              { return OneAbs(b.first) < OneAbs(a.first); } );
            const Int numKeptPairs = numShiftsIdeal / 2;
            const Int numDroppedPairs = numShiftPairs - numKeptPairs;
            const Int shiftBegNew = winEnd - numShiftsIdeal;
            for( Int i=0; i<numKeptPairs; ++i )
            {
                wLoc(shiftBegNew+2*i) = shiftPairs[i+numDroppedPairs].first;
                wLoc(shiftBegNew+2*i+1) = shiftPairs[i+numDroppedPairs].second;
            }
        }
        shiftBeg = winEnd - numShiftsIdeal;
    }
    else
    {
        Int numShifts = winEnd-shiftBeg;
        numShifts = numShifts - Mod(numShifts,2);
        shiftBeg = winEnd - numShifts;
    }
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
    EL_DEBUG_CSE
    typedef Complex<Real> F;

    // Compute the ideal number of (even) shifts to apply this sweep
    Int numShiftsIdeal = Min( numShiftsRec, Max(2,newIterWinSize-1) );
    numShiftsIdeal = numShiftsIdeal - Mod(numShiftsIdeal,2);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts.
        // For some reason, LAPACK suggests only using a single exceptional
        // shift for complex matrices.
        shiftBeg = winEnd-numShiftsIdeal;
        const Real exceptShift0(Real(4)/Real(3));
        for( Int i=winEnd-1; i>=shiftBeg+1; i-=2 )
            w(i-1) = w(i) = H(i,i) + exceptShift0*OneAbs(H(i,i-1));
    }
    else
    {
        // TODO(poulson): Come up with a principled switching point
        if( winEnd-shiftBeg <= numShiftsIdeal/2 )
        {
            // Grab more shifts from another trailing submatrix
            shiftBeg = winEnd-numShiftsIdeal;
            auto shiftsInd = IR(shiftBeg,shiftBeg+numShiftsIdeal);
            auto HShifts = H(shiftsInd,shiftsInd);
            auto wShifts = w(shiftsInd,ALL);
            auto HShiftsCopy( HShifts );

            auto ctrlShifts( ctrl );
            ctrlShifts.winBeg = 0;
            ctrlShifts.winEnd = numShiftsIdeal;
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
    }
    if( winEnd-shiftBeg == 2 )
    {
        // Use the same shift twice; we choose the one closest to the
        // bottom-right entry, as it is our best guess as to the smallest
        // eigenvalue
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

    if( winEnd-shiftBeg > numShiftsIdeal )
    {
        // We have too many shifts and need to eliminate some
        if( ctrl.sortShifts )
        {
            // Sort the shifts into descending order
            std::sort
            ( w.Buffer()+shiftBeg, w.Buffer()+winEnd,
              [](const Complex<Real>& a, const Complex<Real>& b)
              { return OneAbs(b) < OneAbs(a); } );
        }
        shiftBeg = winEnd-numShiftsIdeal;
    }
    else
    {
        Int numShifts = winEnd-shiftBeg;
        numShifts = numShifts - Mod(numShifts,2);
        shiftBeg = winEnd - numShifts;
    }
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
  const DistMatrix<Complex<Real>,MC,MR,BLOCK>& H,
  const Matrix<Complex<Real>>& hMainWin,
  const Matrix<Complex<Real>>& hSubWin,
        DistMatrix<Complex<Real>,STAR,STAR>& w,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> F;
    auto& wLoc = w.Matrix();
    const Grid& grid = H.Grid();

    // Compute the ideal number of (even) shifts to apply this sweep
    Int numShiftsIdeal = Min( numShiftsRec, Max(2,newIterWinSize-1) );
    numShiftsIdeal = numShiftsIdeal - Mod(numShiftsIdeal,2);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts.
        // For some reason, LAPACK suggests only using a single exceptional
        // shift for complex matrices.
        shiftBeg = winEnd-numShiftsIdeal;
        const Real exceptShift0(Real(4)/Real(3));
        for( Int i=winEnd-1; i>=shiftBeg+1; i-=2 )
            wLoc(i-1) = wLoc(i) = hMainWin(i-winBeg) +
              exceptShift0*OneAbs(hSubWin(i-1-winBeg));
    }
    else
    {
        // TODO(poulson): Come up with a principled switching point
        if( winEnd-shiftBeg <= numShiftsIdeal/2 )
        {
            // Grab more shifts from another trailing submatrix
            shiftBeg = winEnd-numShiftsIdeal;
            auto shiftsInd = IR(shiftBeg,shiftBeg+numShiftsIdeal);
            auto HShifts = H(shiftsInd,shiftsInd);
            auto wShifts = w(shiftsInd,ALL);

            auto ctrlShifts( ctrl );
            ctrlShifts.winBeg = 0;
            ctrlShifts.winEnd = numShiftsIdeal;
            ctrlShifts.fullTriangle = false;
            ctrlShifts.demandConverged = false;

            const Int numUnconverged =
              multibulge::ConsistentlyComputeEigenvalues
              ( HShifts, wShifts, ctrlShifts );
            shiftBeg += numUnconverged;
            if( shiftBeg >= winEnd-1 )
            {
                // This should be very rare; use eigenvalues of 2x2
                F eta00 = hMainWin(winEnd-2-winBeg);
                F eta01 = H.Get(winEnd-2,winEnd-1); // This involves a bcast
                F eta10 = hSubWin(winEnd-2-winBeg);
                F eta11 = hMainWin(winEnd-1-winBeg);
                if( grid.VCRank() == 0 )
                {
                    schur::TwoByTwo
                    ( eta00, eta01,
                      eta10, eta11,
                      wLoc(winEnd-2), wLoc(winEnd-1) );
                }
                auto wSub = wLoc( IR(winEnd-2,winEnd), ALL );
                El::Broadcast( wSub, grid.VCComm(), 0 );
                shiftBeg = winEnd-2;
            }
        }
    }
    if( winEnd-shiftBeg == 2 )
    {
        // Use the same shift twice; we choose the one closest to the
        // bottom-right entry, as it is our best guess as to the smallest
        // eigenvalue
        if( Abs(wLoc(winEnd-1).real()-hMainWin(winEnd-1-winBeg)) <
            Abs(wLoc(winEnd-2).real()-hMainWin(winEnd-1-winBeg)) )
        {
            wLoc(winEnd-2) = wLoc(winEnd-1);
        }
        else
        {
            wLoc(winEnd-1) = wLoc(winEnd-2);
        }
    }

    if( winEnd-shiftBeg > numShiftsIdeal )
    {
        // We have too many shifts and need to eliminate some
        if( ctrl.sortShifts )
        {
            // Sort the shifts into descending order
            std::sort
            ( w.Buffer()+shiftBeg, w.Buffer()+winEnd,
              [](const Complex<Real>& a, const Complex<Real>& b)
              { return OneAbs(b) < OneAbs(a); } );
        }
        shiftBeg = winEnd-numShiftsIdeal;
    }
    else
    {
        Int numShifts = winEnd-shiftBeg;
        numShifts = numShifts - Mod(numShifts,2);
        shiftBeg = winEnd - numShifts;
    }
    return shiftBeg;
}

} // namespace aed
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_AED_MODIFY_SHIFTS_HPP
