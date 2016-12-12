/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./HessenbergSchur/Util/MakeSubdiagonalReal.hpp"
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
    EL_DEBUG_CSE
    const Int n = H.Height();
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    ctrlMod.winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    ctrlMod.wantSchurVecs = true;
    if( !ctrl.accumulateSchurVecs )
    {
        Identity( Z, n, n );
        ctrlMod.accumulateSchurVecs = true;
    }

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

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
HessenbergSchurInfo
ScaLAPACKHelper
( AbstractDistMatrix<F>& HPre,
  AbstractDistMatrix<Complex<Base<F>>>& wPre,
  AbstractDistMatrix<F>& ZPre,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    AssertScaLAPACKSupport();
    HessenbergSchurInfo info;

    ProxyCtrl proxyCtrl;
    proxyCtrl.colConstrain = true;
    proxyCtrl.rowConstrain = true;
    proxyCtrl.blockHeight = ctrl.blockHeight;
    proxyCtrl.blockWidth = ctrl.blockHeight;
    proxyCtrl.colAlign = 0;
    proxyCtrl.rowAlign = 0;
    proxyCtrl.colCut = 0;
    proxyCtrl.rowCut = 0;

    DistMatrixReadWriteProxy<F,F,MC,MR,BLOCK> HProx( HPre, proxyCtrl );
    auto& H = HProx.Get();
    const Int n = H.Height();

    DistMatrixReadWriteProxy<F,F,MC,MR,BLOCK> ZProx( ZPre, proxyCtrl );
    auto& Z = ZProx.Get();
    if( !ctrl.accumulateSchurVecs )
        Identity( Z, n, n );

    DistMatrixWriteProxy<Complex<Real>,Complex<Real>,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

#ifdef EL_HAVE_SCALAPACK
    bool accumulateSchurVecs=true;
    bool useAED = ( ctrl.alg == HESSENBERG_SCHUR_AED );

    auto descH = FillDesc( H );
    scalapack::HessenbergSchur
    ( n,
      H.Buffer(), descH.data(),
      w.Buffer(),
      Z.Buffer(), descH.data(),
      ctrl.fullTriangle,
      accumulateSchurVecs,
      useAED );
    // TODO(poulson): Find a way to fill in info?
#endif

    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
HessenbergSchurInfo
ScaLAPACKHelper
( AbstractDistMatrix<F>& HPre,
  AbstractDistMatrix<Complex<Base<F>>>& wPre,
  AbstractDistMatrix<F>& ZPre,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    LogicError("ScalapackHelper should never be called for this datatype");
    HessenbergSchurInfo info;
    return info;
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
HessenbergSchurInfo
ScaLAPACKHelper
( AbstractDistMatrix<F>& HPre,
  AbstractDistMatrix<Complex<Base<F>>>& wPre,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    AssertScaLAPACKSupport();
    HessenbergSchurInfo info;

    ProxyCtrl proxyCtrl;
    proxyCtrl.colConstrain = true;
    proxyCtrl.rowConstrain = true;
    proxyCtrl.blockHeight = ctrl.blockHeight;
    proxyCtrl.blockWidth = ctrl.blockHeight;
    proxyCtrl.colAlign = 0;
    proxyCtrl.rowAlign = 0;
    proxyCtrl.colCut = 0;
    proxyCtrl.rowCut = 0;

    DistMatrixReadWriteProxy<F,F,MC,MR,BLOCK> HProx( HPre, proxyCtrl );
    auto& H = HProx.Get();
    const Int n = H.Height();

    DistMatrixWriteProxy<Complex<Real>,Complex<Real>,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

#ifdef EL_HAVE_SCALAPACK
    bool useAED = ( ctrl.alg == HESSENBERG_SCHUR_AED );

    auto descH = FillDesc( H );
    scalapack::HessenbergSchur
    ( n,
      H.Buffer(), descH.data(),
      w.Buffer(),
      ctrl.fullTriangle, useAED );
    // TODO(poulson): Find a way to fill in info?
#endif

    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
HessenbergSchurInfo
ScaLAPACKHelper
( AbstractDistMatrix<F>& HPre,
  AbstractDistMatrix<Complex<Base<F>>>& wPre,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    LogicError("ScalapackHelper should never be called for this datatype");
    HessenbergSchurInfo info;
    return info;
}

} // namespace hess_schur

template<typename F>
HessenbergSchurInfo
HessenbergSchur
( AbstractDistMatrix<F>& HPre,
  AbstractDistMatrix<Complex<Base<F>>>& wPre,
  AbstractDistMatrix<F>& ZPre,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadWriteProxy<F,F,MC,MR,BLOCK> HProx( HPre );
    auto& H = HProx.Get();

    DistMatrixWriteProxy<Complex<Real>,Complex<Real>,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = H.ColAlign();
    proxCtrl.rowAlign = H.RowAlign();
    proxCtrl.colCut = H.ColCut();
    proxCtrl.rowCut = H.RowCut();
    // Technically, the 'Read' portion of the proxy is only needed if
    // ctrl.accumulateSchurVecs is true.
    DistMatrixReadWriteProxy<F,F,MC,MR,BLOCK> ZProx( ZPre, proxCtrl );
    auto& Z = ZProx.Get();

    const Int n = H.Height();
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    ctrlMod.winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    ctrlMod.wantSchurVecs = true;
    if( !ctrl.accumulateSchurVecs )
    {
        Identity( Z, n, n );
        ctrlMod.accumulateSchurVecs = true;
    }

    if( ctrlMod.scalapack )
    {
#ifdef EL_HAVE_SCALAPACK
        if( IsBlasScalar<F>::value )
            return hess_schur::ScaLAPACKHelper( H, w, Z, ctrlMod );
        else
            Output("Warning: ScaLAPACK is not supported for this datatype");
#else
        Output("Warning: Elemental was not configured with ScaLAPACK support");
#endif
    }

    if( ctrl.alg == HESSENBERG_SCHUR_AED )
    {
        return hess_schur::AED( H, w, Z, ctrlMod );
    }
    else
    {
        return hess_schur::MultiBulge( H, w, Z, ctrlMod );
    }
}

template<typename F>
HessenbergSchurInfo
HessenbergSchur
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& w,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
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

template<typename F>
HessenbergSchurInfo
HessenbergSchur
( AbstractDistMatrix<F>& HPre,
  AbstractDistMatrix<Complex<Base<F>>>& wPre,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadWriteProxy<F,F,MC,MR,BLOCK> HProx( HPre );
    auto& H = HProx.Get();

    DistMatrixWriteProxy<Complex<Real>,Complex<Real>,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    DistMatrix<F,MC,MR,BLOCK> Z(H.Grid());
    Z.AlignWith( H );

    const Int n = H.Height();
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    ctrlMod.winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    ctrlMod.wantSchurVecs = false;

    if( ctrlMod.scalapack )
    {
#ifdef EL_HAVE_SCALAPACK
        if( IsBlasScalar<F>::value )
            return hess_schur::ScaLAPACKHelper( H, w, ctrlMod );
        else
            Output("Warning: ScaLAPACK is not supported for this datatype");
#else
        Output("Warning: Elemental was not configured with ScaLAPACK support");
#endif
    }

    if( ctrl.alg == HESSENBERG_SCHUR_AED )
    {
        return hess_schur::AED( H, w, Z, ctrlMod );
    }
    else
    {
        return hess_schur::MultiBulge( H, w, Z, ctrlMod );
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
    EL_DEBUG_CSE
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
( DistMatrix<Real,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Real>,STAR,STAR>& shifts,
  DistMatrix<Real,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd-winBeg;
    auto ctrlMod( ctrl );
    ctrlMod.winBeg = winBeg;
    ctrlMod.winEnd = winEnd;

    const Int numShifts = shifts.Height();
    multibulge::PairShifts( shifts.Matrix() );

    const Int remainder = (numShifts % 2);
    if( remainder == 1 )
    {
        LogicError
        ("Remainder shifts are not yet supported for distributed sweeps");
    }
    auto shiftsEven = shifts(IR(remainder,END),ALL);

    if( winSize >= 4 )
    {
        multibulge::Sweep( H, shiftsEven, Z, ctrlMod );
    }
    else
    {
        // Sweep in pairs
        LogicError("Distributed pair sweeps are not yet supported");
    }
}

template<typename Real>
void SweepHelper
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& shifts,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
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
        util::MakeSubdiagonalReal( H, Z, ctrl );
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
            util::MakeSubdiagonalReal( H, Z, ctrl );
        }
        for( Int i=0; i<numShifts-remainder; ++i )
            single_shift::SweepOpt( H, shiftsEven(i), Z, ctrlMod );
    }
}

template<typename Real>
void SweepHelper
( DistMatrix<Complex<Real>,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Real>,STAR,STAR>& shifts,
  DistMatrix<Complex<Real>,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
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
        util::MakeSubdiagonalReal( H, Z, ctrl );
        LogicError("Single-shift distributed sweeps are not yet supported");
    }

    auto shiftsEven = shifts(IR(remainder,END),ALL);
    if( winSize >= 4 )
    {
        multibulge::Sweep( H, shiftsEven, Z, ctrlMod );
    }
    else
    {
        if( remainder == 0 )
        {
            // We have not already forced the subdiagonal to be real, which is
            // required for single_shift::SweepOpt to properly function
            util::MakeSubdiagonalReal( H, Z, ctrl );
        }
        LogicError("Single-shift distributed sweeps are not yet supported");
    }
}

template<typename F>
void Sweep
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& shifts,
  Matrix<F>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    SweepHelper( H, shifts, Z, ctrl );
}

template<typename F>
void Sweep
( DistMatrix<F,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts,
  DistMatrix<F,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    SweepHelper( H, shifts, Z, ctrl );
}

} // namespace hess_schur

#define PROTO(F) \
  template HessenbergSchurInfo HessenbergSchur \
  ( Matrix<F>& H, \
    Matrix<Complex<Base<F>>>& w, \
    const HessenbergSchurCtrl& ctrl ); \
  template HessenbergSchurInfo HessenbergSchur \
  ( AbstractDistMatrix<F>& H, \
    AbstractDistMatrix<Complex<Base<F>>>& w, \
    const HessenbergSchurCtrl& ctrl ); \
  template HessenbergSchurInfo HessenbergSchur \
  ( Matrix<F>& H, \
    Matrix<Complex<Base<F>>>& w, \
    Matrix<F>& Z, \
    const HessenbergSchurCtrl& ctrl ); \
  template HessenbergSchurInfo HessenbergSchur \
  ( AbstractDistMatrix<F>& H, \
    AbstractDistMatrix<Complex<Base<F>>>& w, \
    AbstractDistMatrix<F>& Z, \
    const HessenbergSchurCtrl& ctrl ); \
  template void hess_schur::Sweep \
  ( Matrix<F>& H, \
    Matrix<Complex<Base<F>>>& shifts, \
    Matrix<F>& Z, \
    const HessenbergSchurCtrl& ctrl ); \
  template void hess_schur::Sweep \
  ( DistMatrix<F,MC,MR,BLOCK>& H, \
    DistMatrix<Complex<Base<F>>,STAR,STAR>& shifts, \
    DistMatrix<F,MC,MR,BLOCK>& Z, \
    const HessenbergSchurCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
