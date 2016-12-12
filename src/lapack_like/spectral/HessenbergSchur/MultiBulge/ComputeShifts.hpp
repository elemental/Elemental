/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_COMPUTE_SHIFTS_HPP
#define EL_HESS_SCHUR_MULTIBULGE_COMPUTE_SHIFTS_HPP

#include "../Util/Gather.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

// Return the number of unconverged eigenvalues
template<typename Field>
Int ConsistentlyComputeEigenvalues
( const DistMatrix<Field,MC,MR,BLOCK>& H,
        DistMatrix<Complex<Base<Field>>,STAR,STAR>& w,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    // Because double-precision floating-point computation is often
    // non-deterministic due to extra-precision computation being frequent but
    // not guaranteed, we must be careful to not allow this non-determinism to
    // be amplified by the forward instability of Francis sweeps.
    const Grid& grid = H.Grid();
    const int owner = H.Owner(0,0);
    DistMatrix<Field,CIRC,CIRC> H_CIRC_CIRC( grid, owner );
    H_CIRC_CIRC = H;
    w.Resize( H.Height(), 1 );
    Int numUnconverged = 0;
    if( H_CIRC_CIRC.CrossRank() == H_CIRC_CIRC.Root() )
    {
        auto info = HessenbergSchur( H_CIRC_CIRC.Matrix(), w.Matrix(), ctrl );
        numUnconverged = info.numUnconverged;
    }
    // TODO(poulson): Combine the two broadcasts to reduce the latency cost?
    El::Broadcast( w.Matrix(), H_CIRC_CIRC.CrossComm(), H_CIRC_CIRC.Root() );
    if( !ctrl.demandConverged )
        mpi::Broadcast
        ( numUnconverged, H_CIRC_CIRC.Root(), H_CIRC_CIRC.CrossComm() );
    return numUnconverged;
}

template<typename Real>
Int ComputeShifts
( const Matrix<Real>& H,
        Matrix<Complex<Real>>& w,
        Int iterBeg,
        Int winBeg,
        Int winEnd,
        Int numShiftsRec,
        Int numIterSinceDeflation,
        Int numStaleIterBeforeExceptional,
  const HessenbergSchurCtrl& ctrlShifts )
{
    EL_DEBUG_CSE

    const Int shiftBeg = Max(iterBeg,winEnd-numShiftsRec);
    const Int numShifts = winEnd - shiftBeg;
    auto shiftInd = IR(shiftBeg,winEnd);
    auto wShifts = w(shiftInd,ALL);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts
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
        // Compute the eigenvalues of the bottom-right window
        auto HShifts = H(shiftInd,shiftInd);
        auto HShiftsCopy( HShifts );
        HessenbergSchur( HShiftsCopy, wShifts, ctrlShifts );
    }

    if( winBeg-shiftBeg == 2 )
    {
        // Use a single real shift twice instead of using two separate
        // real shifts; we choose the one closest to the bottom-right
        // entry, as it is our best guess as to the smallest eigenvalue
        if( wShifts(numShifts-1).imag() == Real(0) )
        {
            const Real eta11 = H(winEnd-1,winEnd-1);
            if( Abs(wShifts(1).real()-eta11) < Abs(wShifts(0).real()-eta11) )
                wShifts(0) = wShifts(1);
            else
                wShifts(1) = wShifts(0);
        }
    }

    return shiftBeg;
}

template<typename Real>
Int ComputeShifts
( const DistMatrix<Real,MC,MR,BLOCK>& H,
        DistMatrix<Complex<Real>,STAR,STAR>& w,
        Int iterBeg,
        Int winBeg,
        Int winEnd,
        Int numShiftsRec,
        Int numIterSinceDeflation,
        Int numStaleIterBeforeExceptional,
  const HessenbergSchurCtrl& ctrlShifts )
{
    EL_DEBUG_CSE

    const Int shiftBeg = Max(iterBeg,winEnd-numShiftsRec);
    const Int numShifts = winEnd - shiftBeg;
    auto shiftInd = IR(shiftBeg,winEnd);
    auto wShifts = w(shiftInd,ALL);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts
        const Real exceptShift0(Real(4)/Real(3)),
                   exceptShift1(-Real(7)/Real(16));

        // Get a full copy of the bidiagonal of the bottom-right section of H
        const Int subStart = Max(shiftBeg-1,winBeg);
        auto subInd = IR(subStart,winEnd);
        DistMatrix<Real,STAR,STAR> hMain(H.Grid()), hSub(H.Grid());
        util::GatherBidiagonal( H, subInd, hMain, hSub );
        const auto& hMainLoc = hMain.LockedMatrix();
        const auto& hSubLoc = hSub.LockedMatrix();

        for( Int i=winEnd-1; i>=subStart+2; i-=2 )
        {
            const Int iRel = i - subStart;
            const Real scale = Abs(hSubLoc(iRel-1)) + Abs(hSubLoc(iRel-2));
            Real eta00 = exceptShift0*scale + hMainLoc(iRel);
            Real eta01 = scale;
            Real eta10 = exceptShift1*scale;
            Real eta11 = eta00;
            Complex<Real> omega0, omega1;
            schur::TwoByTwo
            ( eta00, eta01,
              eta10, eta11,
              omega0, omega1 );
            w.Set( i-1, 0, omega0 );
            w.Set( i,   0, omega1 );
        }
        if( shiftBeg == winBeg )
        {
            w.Set( shiftBeg,   0, hMainLoc(1) );
            w.Set( shiftBeg+1, 0, hMainLoc(1) );
        }
    }
    else
    {
        ConsistentlyComputeEigenvalues
        ( H(shiftInd,shiftInd), wShifts, ctrlShifts );
    }

    if( numShifts == 2 )
    {
        // Use a single real shift twice instead of using two separate
        // real shifts; we choose the one closest to the bottom-right
        // entry, as it is our best guess as to the smallest eigenvalue
        const Complex<Real> omega0 = wShifts.GetLocal(0,0);
        const Complex<Real> omega1 = wShifts.GetLocal(1,0);
        if( omega1.imag() == Real(0) )
        {
            const Real eta11 = H.Get(winEnd-1,winEnd-1);
            if( Abs(omega1.real()-eta11) < Abs(omega0.real()-eta11) )
                wShifts.Set( 0, 0, omega1 );
            else
                wShifts.Set( 1, 0, omega0 );
        }
    }

    return shiftBeg;
}

template<typename Real>
Int ComputeShifts
( const Matrix<Complex<Real>>& H,
        Matrix<Complex<Real>>& w,
        Int iterBeg,
        Int winBeg,
        Int winEnd,
        Int numShiftsRec,
        Int numIterSinceDeflation,
        Int numStaleIterBeforeExceptional,
  const HessenbergSchurCtrl& ctrlShifts )
{
    EL_DEBUG_CSE

    const Int shiftBeg = Max(iterBeg,winEnd-numShiftsRec);
    const Int numShifts = winEnd - shiftBeg;
    auto shiftInd = IR(shiftBeg,winEnd);
    auto wShifts = w(shiftInd,ALL);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts
        //
        // For some reason, LAPACK suggests only using a single exceptional
        // shift for complex matrices.
        const Real exceptShift0(Real(4)/Real(3));
        for( Int i=winEnd-1; i>=shiftBeg+1; i-=2 )
            w(i-1) = w(i) = H(i,i) + exceptShift0*OneAbs(H(i,i-1));
    }
    else
    {
        // Compute the eigenvalues of the bottom-right window
        auto HShifts = H(shiftInd,shiftInd);
        auto HShiftsCopy( HShifts );
        HessenbergSchur( HShiftsCopy, wShifts, ctrlShifts );
    }

    if( numShifts == 2 )
    {
        // Use the same shift twice; we choose the one closest to the
        // bottom-right entry, as it is our best guess as to the smallest
        // eigenvalue
        const Complex<Real> eta11 = H(winEnd-1,winEnd-1);
        if( Abs(wShifts(1)-eta11) < Abs(wShifts(0)-eta11) )
            wShifts(0) = wShifts(1);
        else
            wShifts(1) = wShifts(0);
    }

    return shiftBeg;
}

template<typename Real>
Int ComputeShifts
( const DistMatrix<Complex<Real>,MC,MR,BLOCK>& H,
        DistMatrix<Complex<Real>,STAR,STAR>& w,
        Int iterBeg,
        Int winBeg,
        Int winEnd,
        Int numShiftsRec,
        Int numIterSinceDeflation,
        Int numStaleIterBeforeExceptional,
  const HessenbergSchurCtrl& ctrlShifts )
{
    EL_DEBUG_CSE

    const Int shiftBeg = Max(iterBeg,winEnd-numShiftsRec);
    const Int numShifts = winEnd - shiftBeg;
    auto shiftInd = IR(shiftBeg,winEnd);
    auto wShifts = w(shiftInd,ALL);

    if( numIterSinceDeflation > 0 &&
        Mod(numIterSinceDeflation,numStaleIterBeforeExceptional) == 0 )
    {
        // Use exceptional shifts
        //
        // For some reason, LAPACK suggests only using a single exceptional
        // shift for complex matrices.
        const Real exceptShift0(Real(4)/Real(3));

        // Gather the relevant bidiagonal of H
        DistMatrix<Complex<Real>,STAR,STAR> hMain(H.Grid()), hSub(H.Grid());
        util::GatherBidiagonal( H, shiftInd, hMain, hSub );
        const auto& hMainLoc = hMain.LockedMatrix();
        const auto& hSubLoc = hSub.LockedMatrix();

        for( Int i=winEnd-1; i>=shiftBeg+1; i-=2 )
        {
            const Int iRel = i - shiftBeg;
            const Complex<Real> shift = hMainLoc(iRel) +
              exceptShift0*OneAbs(hSubLoc(iRel-1));
            wShifts.Set( iRel-1, 0, shift );
            wShifts.Set( iRel,   0, shift );
        }
    }
    else
    {
        ConsistentlyComputeEigenvalues
        ( H(shiftInd,shiftInd), wShifts, ctrlShifts );
    }

    if( numShifts == 2 )
    {
        // Use the same shift twice; we choose the one closest to the
        // bottom-right entry, as it is our best guess as to the smallest
        // eigenvalue
        const Complex<Real> omega0 = wShifts.GetLocal(0,0);
        const Complex<Real> omega1 = wShifts.GetLocal(1,0);
        const Complex<Real> eta11 = H.Get(winEnd-1,winEnd-1);
        if( Abs(omega1-eta11) < Abs(omega0-eta11) )
            wShifts.Set( 0, 0, omega1 );
        else
            wShifts.Set( 1, 0, omega0 );
    }

    return shiftBeg;
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_COMPUTE_SHIFTS_HPP
