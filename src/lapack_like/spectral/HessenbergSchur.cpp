/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./HessenbergSchur/Simple.hpp"
#include "./HessenbergSchur/MultiBulge.hpp"
#include "./HessenbergSchur/AED.hpp"

namespace El {

template<typename F>
HessenbergSchurInfo
HessenbergSchur
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    const Int n = H.Height();
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    ctrlMod.winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    ctrlMod.wantSchurVecs = true;

    if( ctrl.alg == HESSENBERG_SCHUR_AED )
    {
        return hess_schur::AED( H, w, Z, ctrlMod );
    }
    else if( ctrl.alg == HESSENBERG_SCHUR_MULTIBULGE )
    {
        return hess_schur::MultiBulge( H, w, Z, ctrlMod );
    }
    else
    {
        return hess_schur::Simple( H, w, Z, ctrlMod );
    }
}

template<typename F>
HessenbergSchurInfo
HessenbergSchur
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& w,
  const HessenbergSchurCtrl& ctrl )
{
    const Int n = H.Height();
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    ctrlMod.winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    ctrlMod.wantSchurVecs = false;

    Matrix<F> Z;
    if( ctrl.alg == HESSENBERG_SCHUR_AED )
    {
        return hess_schur::AED( H, w, Z, ctrlMod );
    }
    else if( ctrl.alg == HESSENBERG_SCHUR_MULTIBULGE )
    {
        return hess_schur::MultiBulge( H, w, Z, ctrlMod );
    }
    else
    {
        return hess_schur::Simple( H, w, Z, ctrlMod );
    }
}

namespace hess_schur {

template<typename Real>
void SweepHelper
( Matrix<Real>& H,
  Matrix<Complex<Real>>& shifts,
  Matrix<Real>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd ); 
    const Int winSize = winEnd-winBeg;
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = winBeg;
    ctrlMod.winEnd = winEnd;

    const Int numShifts = shifts.Height();
    multibulge::PairShifts( shifts );

    Matrix<Real> U, W, WAccum;
    const Int remainder = (numShifts % 2);
    if( remainder == 1 )
    {
        // We will separately apply the odd shift and its conjugate
        double_shift::SweepOpt( H, shifts(0), Conj(shifts(0)), Z, ctrlMod );
    }
    auto shiftsEven = shifts(IR(remainder,END),ALL);

    if( winSize >= 4 )
    {
        multibulge::Sweep( H, shiftsEven, Z, U, W, WAccum, ctrlMod );
    }
    else
    {
        // Sweep in pairs
        for( Int i=0; i<numShifts/2; ++i )
            double_shift::SweepOpt
            ( H, shiftsEven(2*i), shiftsEven(2*i+1), Z, ctrlMod );
    }
}

template<typename Real>
void SweepHelper
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& shifts,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = H.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd ); 
    const Int winSize = winEnd-winBeg;
    const Int numShifts = shifts.Height();
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = winBeg;
    ctrlMod.winEnd = winEnd;

    const Int remainder = (numShifts % 2);
    // We separately apply an odd shift if one exists
    if( remainder == 1 )
    {
        // The optimized sweeping procedure assumes/maintains a real subdiagonal
        MakeSubdiagonalReal( H, Z, ctrl );
        single_shift::SweepOpt( H, shifts(0), Z, ctrlMod );
    }

    auto shiftsEven = shifts(IR(remainder,END),ALL);
    if( winSize >= 4 )
    {
        Matrix<Complex<Real>> U, W, WAccum;
        multibulge::Sweep( H, shiftsEven, Z, U, W, WAccum, ctrlMod );
    }
    else
    {
        if( remainder == 0 )
        {
            // We have not already forced the subdiagonal to be real, which is
            // required for single_shift::SweepOpt to properly function
            MakeSubdiagonalReal( H, Z, ctrl );
        }
        for( Int i=0; i<numShifts-remainder; ++i )
            single_shift::SweepOpt( H, shiftsEven(i), Z, ctrlMod );
    }
}

template<typename F>
void Sweep
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& shifts,
  Matrix<F>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    SweepHelper( H, shifts, Z, ctrl );
}

} // namespace hess_schur

#define PROTO(F) \
  template HessenbergSchurInfo HessenbergSchur \
  ( Matrix<F>& H, \
    Matrix<Complex<Base<F>>>& w, \
    const HessenbergSchurCtrl& ctrl ); \
  template HessenbergSchurInfo HessenbergSchur \
  ( Matrix<F>& H, \
    Matrix<Complex<Base<F>>>& w, \
    Matrix<F>& Z, \
    const HessenbergSchurCtrl& ctrl ); \
  template void hess_schur::Sweep \
  ( Matrix<F>& H, \
    Matrix<Complex<Base<F>>>& shifts, \
    Matrix<F>& Z, \
    const HessenbergSchurCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
