/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_COMPUTE_SHIFTS_HPP
#define EL_HESS_SCHUR_MULTIBULGE_COMPUTE_SHIFTS_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename F>
void GatherBidiagonal
( const DistMatrix<F,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<F,STAR,STAR>& hMainWin,
        DistMatrix<Base<F>,STAR,STAR>& hSubWin )
{
    DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    DEBUG_ONLY(
      if( H.BlockHeight() != H.BlockWidth() )
          LogicError("Assumed square distribution blocks");
      if( H.ColCut() != H.RowCut() )
          LogicError("Assumed symmetric cuts");
      if( blockSize < 2 )
          LogicError("Assumed blocks of size at least two");
    )
    Zeros( hMainWin, winSize, 1 );
    Zeros( hSubWin, winSize-1, 1 );

    // TODO(poulson): Optimize these redistributions
    hMainWin.Reserve( winSize );
    hSubWin.Reserve( winSize-1 );

    Int localRowOffset = H.LocalRowOffset( winInd.beg );
    Int localColOffset = H.LocalColOffset( winInd.beg );
    int rowOwner = H.RowOwner( winInd.beg );
    int colOwner = H.ColOwner( winInd.beg );

    Int thisCut = Mod( H.ColCut()+winInd.beg, H.BlockHeight() );
    Int winIndex = 0;
    while( winIndex < winSize )
    {
        const Int thisBlockSize = Min( blockSize-thisCut, winSize-winIndex );
        const int nextRowOwner = Mod( rowOwner+1, grid.Height() );
        const int nextColOwner = Mod( colOwner+1, grid.Width() );

        if( colOwner == grid.Col() )
        {
            const bool haveLastSub = (winSize > winIndex+thisBlockSize);
            if( rowOwner == grid.Row() )
            {
                for( Int offset=0; offset<thisBlockSize; ++offset )
                {
                    hMainWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset,localColOffset+offset) );
                    if( offset < thisBlockSize-1 ||
                        (haveLastSub && grid.Height()==1) )
                        hSubWin.QueueUpdate
                        ( winIndex+offset, 0,
                          RealPart(HLoc(localRowOffset+offset+1,
                                        localColOffset+offset)) );
                }
            }
            else if( haveLastSub && nextRowOwner == grid.Row() )
            {
                // Grab the entry in the top-right of this block
                hSubWin.QueueUpdate
                ( winIndex+thisBlockSize-1, 0,
                  RealPart(HLoc(localRowOffset,
                                localColOffset+thisBlockSize-1)) );
            }
        }

        winIndex += thisBlockSize;
        thisCut = 0;
        if( colOwner == grid.Col() )
            localColOffset += thisBlockSize;
        if( rowOwner == grid.Row() )
            localRowOffset += thisBlockSize;
        rowOwner = nextRowOwner;
        colOwner = nextColOwner;
    }
    hMainWin.ProcessQueues();
    hSubWin.ProcessQueues();
}

template<typename F>
void GatherTridiagonal
( const DistMatrix<F,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<F,STAR,STAR>& hMainWin,
        DistMatrix<Base<F>,STAR,STAR>& hSubWin,
        DistMatrix<F,STAR,STAR>& hSuperWin )
{
    DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    DEBUG_ONLY(
      if( H.BlockHeight() != H.BlockWidth() )
          LogicError("Assumed square distribution blocks");
      if( H.ColCut() != H.RowCut() )
          LogicError("Assumed symmetric cuts");
      if( blockSize < 2 )
          LogicError("Assumed blocks of size at least two");
    )
    Zeros( hMainWin, winSize, 1 );
    Zeros( hSubWin, winSize-1, 1 );

    // TODO(poulson): Optimize these redistributions
    hMainWin.Reserve( winSize );
    hSubWin.Reserve( winSize-1 );
    hSuperWin.Reserve( winSize-1 );

    Int localRowOffset = H.LocalRowOffset( winInd.beg );
    Int localColOffset = H.LocalColOffset( winInd.beg );
    int rowOwner = H.RowOwner( winInd.beg );
    int colOwner = H.ColOwner( winInd.beg );

    Int thisCut = Mod( H.ColCut()+winInd.beg, H.BlockHeight() );
    Int winIndex = 0;
    while( winIndex < winSize )
    {
        const Int thisBlockSize = Min( blockSize-thisCut, winSize-winIndex );
        const int nextRowOwner = Mod( rowOwner+1, grid.Height() );
        const int nextColOwner = Mod( colOwner+1, grid.Width() );

        if( colOwner == grid.Col() )
        {
            const bool haveLastOff = (winSize > winIndex+thisBlockSize);
            if( rowOwner == grid.Row() )
            {
                for( Int offset=0; offset<thisBlockSize; ++offset )
                {
                    hMainWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset,localColOffset+offset) );
                    if( offset < thisBlockSize-1 ||
                        (haveLastOff && grid.Height()==1) )
                        hSubWin.QueueUpdate
                        ( winIndex+offset, 0,
                          RealPart(HLoc(localRowOffset+offset+1,
                                        localColOffset+offset)) );
                    if( offset < thisBlockSize-1 ||
                        (haveLastOff && grid.Width()==1) )
                        hSuperWin.QueueUpdate
                        ( winIndex+offset, 0,
                          HLoc(localRowOffset+offset,localColOffset+offset+1) );
                }
            }
            else if( haveLastOff && nextRowOwner == grid.Row() )
            {
                // Grab the entry in the top-right of this block
                hSubWin.QueueUpdate
                ( winIndex+thisBlockSize-1, 0,
                  RealPart(HLoc(localRowOffset,
                                localColOffset+thisBlockSize-1)) );
            }
        }
        else if( rowOwner == grid.Row() && nextColOwner == grid.Col() )
        {
            const bool haveLastOff = (winSize > winIndex+thisBlockSize);
            if( haveLastOff )
            {
                // Grab the entry in the bottom-left of this block
                hSuperWin.QueueUpdate
                ( winIndex+thisBlockSize-1, 0,
                  HLoc(localRowOffset+thisBlockSize-1,localColOffset) );
            }
        }

        winIndex += thisBlockSize;
        thisCut = 0;
        if( colOwner == grid.Col() )
            localColOffset += thisBlockSize;
        if( rowOwner == grid.Row() )
            localRowOffset += thisBlockSize;
        rowOwner = nextRowOwner;
        colOwner = nextColOwner;
    }
    hMainWin.ProcessQueues();
    hSubWin.ProcessQueues();
    hSuperWin.ProcessQueues();
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
    DEBUG_CSE

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
    DEBUG_CSE

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
        DistMatrix<Real,STAR,STAR> hMain(H.Grid());
        DistMatrix<Real,STAR,STAR> hSub(H.Grid());
        GatherBidiagonal( H, shiftInd, hMain, hSub );
        const auto& hMainLoc = hMain.LockedMatrix();
        const auto& hSubLoc = hSub.LockedMatrix();

        for( Int i=winEnd-1; i>=Max(shiftBeg+1,winBeg+2); i-=2 )
        {
            const Int iRel = i - shiftBeg;
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
        // Compute the eigenvalues of the bottom-right window
        DistMatrix<Real,STAR,STAR> HShifts( H(shiftInd,shiftInd) );
        HessenbergSchur( HShifts.Matrix(), wShifts.Matrix(), ctrlShifts );
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
    DEBUG_CSE

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
    DEBUG_CSE

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
        DistMatrix<Complex<Real>,STAR,STAR> hMain(H.Grid());
        DistMatrix<Real,STAR,STAR> hSub(H.Grid());
        GatherBidiagonal( H, shiftInd, hMain, hSub );
        const auto& hMainLoc = hMain.LockedMatrix();
        const auto& hSubLoc = hSub.LockedMatrix();

        for( Int i=winEnd-1; i>=shiftBeg+1; i-=2 )
        {
            const Int iRel = i - shiftBeg;
            const Complex<Real> shift = hMainLoc(iRel) +
              exceptShift0*Abs(hSubLoc(iRel-1));
            wShifts.Set( iRel-1, 0, shift );
            wShifts.Set( iRel,   0, shift );
        }
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
