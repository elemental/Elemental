/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP

// Disable this if using C99's <complex.h> is not acceptable/desired
#define EL_REINTERPRET_COMPLEX

#ifdef EL_REINTERPRET_COMPLEX
#include <complex.h>

namespace El {

// TODO(poulson): Move these into the headers
template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
struct FastComplexHelper { };
template<>
struct FastComplexHelper<float> { typedef float _Complex value; };
template<>
struct FastComplexHelper<double> { typedef double _Complex value; };

template<typename Real>
using FastComplex = typename FastComplexHelper<Real>::value;

} // namespace El

#endif // ifdef EL_REINTERPRET_COMPLEX

namespace El {

template<typename F>
struct IsComplexBlasScalar {
  static const bool value = IsComplex<F>::value && IsBlasScalar<F>::value;
};

} // namespace El

#include "./ApplyReflector.hpp"
#include "./VigilantDeflation.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename F>
void ApplyReflectors
( Matrix<F>& H,
  Int winBeg,
  Int winEnd,
  Int chaseBeg,
  Int packetBeg,
  Int transformRowBeg,
  Int transformColEnd,
  Matrix<F>& Z,
  bool wantSchurVecs,
  Matrix<F>& U,
  const Matrix<F>& W,
  Int firstBulge,
  Int numBulges,
  bool accumulate,
  bool progress )
{
    DEBUG_CSE
    const Int lastBulge = firstBulge + numBulges - 1;
    const Int lastBulgeBeg = packetBeg + 3*lastBulge;
    const bool haveSmallBulge = ( lastBulgeBeg == winEnd-3 );
    DEBUG_ONLY(
      if( lastBulgeBeg > winEnd-3 )
          LogicError("Last bulge starts too late");
    )
    const Int numFullBulges = ( haveSmallBulge ? numBulges-1 : numBulges );
    const Int clippedChaseBeg = Max(chaseBeg,winBeg-1);

    // Apply from the left
    // ===================
    for( Int j=Max(chaseBeg,winBeg); j<transformColEnd; ++j )
    {
        // Avoid re-applying freshly generated reflections
        // (which is only relevant when (j-packetBeg) % 3 = 0)
        const Int applyBulgeEnd =
          Min( firstBulge+numFullBulges, (j-packetBeg+2)/3 );
        for( Int bulge=firstBulge; bulge<applyBulgeEnd; ++bulge )
        {
            const Int bulgeBeg = packetBeg + 3*bulge;
            const F* w = &W(0,bulge);
            ApplyLeftReflector
            ( H(bulgeBeg+1,j), H(bulgeBeg+2,j), H(bulgeBeg+3,j), w );
        }
    }
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )
        for( Int j=bulgeBeg+1; j<transformColEnd; ++j )
            ApplyLeftReflector( H(bulgeBeg+1,j), H(bulgeBeg+2,j), w );
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    const Int UHeight = U.Height();
    const Int ZHeight = Z.Height();
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+3); ++i )
            ApplyRightReflector( H(i,bulgeBeg+1), H(i,bulgeBeg+2), w );

        if( accumulate )
        {
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            for( Int i=0; i<UHeight; ++i ) 
                ApplyRightReflector
                ( U(i,bulgeBegRel+1), U(i,bulgeBegRel+2), w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<ZHeight; ++i )
                ApplyRightReflector( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), w );
        }
    }
    for( Int bulge=firstBulge+numFullBulges-1; bulge>=firstBulge; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+4); ++i )
            ApplyRightReflector
            ( H(i,bulgeBeg+1), H(i,bulgeBeg+2), H(i,bulgeBeg+3), w );

        if( accumulate )
        {
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            for( Int i=0; i<UHeight; ++i ) 
                ApplyRightReflector
                ( U(i,bulgeBegRel+1), U(i,bulgeBegRel+2), U(i,bulgeBegRel+3),
                  w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<ZHeight; ++i )
                ApplyRightReflector
                ( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), Z(i,bulgeBeg+3), w );
        }
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    Int firstVigBulge = firstBulge;
    Int numVigBulges = numBulges;
    if( packetBeg + 3*firstBulge == winBeg-1 )
    {
        // The first bulge needs to be introduced
        ++firstVigBulge;
        --numVigBulges;
    }
    VigilantDeflation
    ( H, winBeg, winEnd, packetBeg, firstVigBulge, numVigBulges, progress );

    // Form the last row of the single-step bulge chase
    //
    //       ~ ~ ~                 ~ ~ ~
    //     -----------          -----------
    //    | B B B B x |        | x x x x x |
    //  ~ | B B B B x |      ~ | x B B B B |
    //  ~ | B B B B x | |->  ~ |   B B B B |.
    //  ~ | B B B B x |      ~ |   B B B B |
    //    |       x x |        |   B B B B |
    //     -----------          -----------
    //
    // The last row is introduced from the transformation
    //
    //  H(k+4,k+1:k+3) -= conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2] [1; nu1; nu2]'.
    //
    // For convenience, since H(k+4,k+1:k+2)=0, we start by computing 
    //
    //   innerProd = conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2]
    //             = conj(tau) H(k+4,k+3) nu2
    //
    // TODO(poulson): Avoid memory allocations for several of the temporaries
    // below in the case where we have heap scalars (e.g., BigFloat).
    const Int lastRowBulgeEnd =
      Min( firstBulge+numFullBulges, (winEnd-packetBeg-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=firstBulge; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        const F& tau = w[0];
        const F& nu1 = w[1];
        const F& nu2 = w[2];

        const F innerProd = Conj(tau)*H(bulgeBeg+4,bulgeBeg+3)*nu2;
        H(bulgeBeg+4,bulgeBeg+1) = -innerProd;
        H(bulgeBeg+4,bulgeBeg+2) = -innerProd*Conj(nu1);
        H(bulgeBeg+4,bulgeBeg+3) -= innerProd*Conj(nu2);
    }
}

// Unfortunately, it seems to be the case that it is noticeably faster
// for this routine to manually inline the data access than to use the 
// (presumably inlined) Matrix::operator()(int,int) calls.
#ifdef EL_REINTERPRET_COMPLEX
template<typename F,typename=DisableIf<IsComplexBlasScalar<F>>>
#else
template<typename F>
#endif
void ApplyReflectorsOpt
( Matrix<F>& H,
  Int winBeg,
  Int winEnd,
  Int chaseBeg,
  Int packetBeg,
  Int transformRowBeg,
  Int transformColEnd,
  Matrix<F>& Z,
  bool wantSchurVecs,
  Matrix<F>& U,
  const Matrix<F>& W,
  Int firstBulge,
  Int numBulges,
  bool accumulate,
  bool progress )
{
    DEBUG_CSE
    const Int lastBulge = firstBulge + numBulges - 1;
    const Int lastBulgeBeg = packetBeg + 3*lastBulge;
    const bool haveSmallBulge = ( lastBulgeBeg == winEnd-3 );
    DEBUG_ONLY(
      if( lastBulgeBeg > winEnd-3 )
          LogicError("Last bulge starts too late");
    )
    const Int numFullBulges = ( haveSmallBulge ? numBulges-1 : numBulges );
    const Int clippedChaseBeg = Max(chaseBeg,winBeg-1);

    F* EL_RESTRICT HBuf = H.Buffer();
    F* EL_RESTRICT ZBuf = Z.Buffer();
    F* EL_RESTRICT UBuf = U.Buffer();
    const F* EL_RESTRICT WBuf = W.LockedBuffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();
    const Int ULDim = U.LDim();
    const Int WLDim = W.LDim();

    // Apply from the left
    // ===================
    for( Int j=Max(chaseBeg,winBeg); j<transformColEnd; ++j )
    {
        // Avoid re-applying freshly generated reflections
        // (which is only relevant when (j-packetBeg) % 3 = 0)
        const Int applyBulgeEnd =
          Min( firstBulge+numFullBulges, (j-packetBeg+2)/3 );
        for( Int bulge=firstBulge; bulge<applyBulgeEnd; ++bulge )
        {
            const Int bulgeBeg = packetBeg + 3*bulge;
            const F* EL_RESTRICT w = &WBuf[bulge*WLDim];
            const F omega0 = w[0];
            const F omega1 = w[1];
            const F omega1Conj = Conj(w[1]);
            const F omega2 = w[2];
            const F omega2Conj = Conj(w[2]);
            F* EL_RESTRICT hj = &HBuf[j*HLDim];
            const F innerProd =
              omega0*(hj[bulgeBeg+1] + omega1Conj*hj[bulgeBeg+2] +
                      omega2Conj*hj[bulgeBeg+3]);
            hj[bulgeBeg+1] -= innerProd;
            hj[bulgeBeg+2] -= innerProd*omega1;
            hj[bulgeBeg+3] -= innerProd*omega2;
        }
    }
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* EL_RESTRICT w = &WBuf[bulge*WLDim];
        const F omega0 = w[0];
        const F omega1 = w[1];
        const F omega1Conj = Conj(w[1]);
        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )
        for( Int j=bulgeBeg+1; j<transformColEnd; ++j )
        {
            F* EL_RESTRICT hj = &HBuf[j*HLDim];
            const F innerProd =
              omega0*(hj[bulgeBeg+1] + omega1Conj*hj[bulgeBeg+2]);
            hj[bulgeBeg+1] -= innerProd;
            hj[bulgeBeg+2] -= innerProd*omega1;
        }
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    const Int UHeight = U.Height();
    const Int ZHeight = Z.Height();
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* EL_RESTRICT w = &WBuf[bulge*WLDim];
        const F omega0Conj = Conj(w[0]);
        const F omega1 = w[1];
        const F omega1Conj = Conj(w[1]);

        F* EL_RESTRICT h1 = &HBuf[(bulgeBeg+1)*HLDim];
        F* EL_RESTRICT h2 = &HBuf[(bulgeBeg+2)*HLDim];
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+3); ++i )
        {
            const F innerProd = omega0Conj*(h1[i] + omega1*h2[i]);
            h1[i] -= innerProd;
            h2[i] -= innerProd*omega1Conj;
        }

        if( accumulate )
        {
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            F* EL_RESTRICT u1 = &UBuf[(bulgeBegRel+1)*ULDim];
            F* EL_RESTRICT u2 = &UBuf[(bulgeBegRel+2)*ULDim];
            for( Int i=0; i<UHeight; ++i ) 
            {
                const F innerProd = omega0Conj*(u1[i] + omega1*u2[i]);
                u1[i] -= innerProd;
                u2[i] -= innerProd*omega1Conj;
            }
        }
        else if( wantSchurVecs )
        {
            F* EL_RESTRICT z1 = &ZBuf[(bulgeBeg+1)*ZLDim];
            F* EL_RESTRICT z2 = &ZBuf[(bulgeBeg+2)*ZLDim];
            for( Int i=0; i<ZHeight; ++i )
            {
                const F innerProd = omega0Conj*(z1[i] + omega1*z2[i]);
                z1[i] -= innerProd;
                z2[i] -= innerProd*omega1Conj;
            }
        }
    }
    // Even though we are logically chasing the bulges starting from the
    // bottom-right, it is typically faster to traverse from the top-left
    for( Int bulge=firstBulge; bulge<firstBulge+numFullBulges; ++bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* EL_RESTRICT w = &WBuf[bulge*WLDim];
        const F omega0Conj = Conj(w[0]);
        const F omega1 = w[1];
        const F omega1Conj = Conj(w[1]);
        const F omega2 = w[2];
        const F omega2Conj = Conj(w[2]);

        F* EL_RESTRICT h1 = &HBuf[(bulgeBeg+1)*HLDim];
        F* EL_RESTRICT h2 = &HBuf[(bulgeBeg+2)*HLDim];
        F* EL_RESTRICT h3 = &HBuf[(bulgeBeg+3)*HLDim];
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+4); ++i )
        {
            const F innerProd =
              omega0Conj*(h1[i] + omega1*h2[i] + omega2*h3[i]);
            h1[i] -= innerProd;
            h2[i] -= innerProd*omega1Conj;
            h3[i] -= innerProd*omega2Conj;
        }

        if( accumulate )
        {
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            F* EL_RESTRICT u1 = &UBuf[(bulgeBegRel+1)*ULDim];
            F* EL_RESTRICT u2 = &UBuf[(bulgeBegRel+2)*ULDim];
            F* EL_RESTRICT u3 = &UBuf[(bulgeBegRel+3)*ULDim];
            for( Int i=0; i<UHeight; ++i ) 
            {
                const F innerProd =
                  omega0Conj*(u1[i] + omega1*u2[i] + omega2*u3[i]);
                u1[i] -= innerProd;
                u2[i] -= innerProd*omega1Conj;
                u3[i] -= innerProd*omega2Conj;
            }
        }
        else if( wantSchurVecs )
        {
            F* EL_RESTRICT z1 = &ZBuf[(bulgeBeg+1)*ZLDim];
            F* EL_RESTRICT z2 = &ZBuf[(bulgeBeg+2)*ZLDim];
            F* EL_RESTRICT z3 = &ZBuf[(bulgeBeg+3)*ZLDim];
            for( Int i=0; i<ZHeight; ++i )
            {
                const F innerProd =
                  omega0Conj*(z1[i] + omega1*z2[i] + omega2*z3[i]);
                z1[i] -= innerProd;
                z2[i] -= innerProd*omega1Conj;
                z3[i] -= innerProd*omega2Conj;
            }
        }
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    Int firstVigBulge = firstBulge;
    Int numVigBulges = numBulges;
    if( packetBeg + 3*firstBulge == winBeg-1 )
    {
        // The first bulge still needs to be introduced
        ++firstVigBulge;
        --numVigBulges;
    }
    VigilantDeflation
    ( H, winBeg, winEnd, packetBeg, firstVigBulge, numVigBulges, progress );

    // Form the last row of the single-step bulge chase
    //
    //       ~ ~ ~                 ~ ~ ~
    //     -----------          -----------
    //    | B B B B x |        | x x x x x |
    //  ~ | B B B B x |      ~ | x B B B B |
    //  ~ | B B B B x | |->  ~ |   B B B B |.
    //  ~ | B B B B x |      ~ |   B B B B |
    //    |       x x |        |   B B B B |
    //     -----------          -----------
    //
    // The last row is introduced from the transformation
    //
    //  H(k+4,k+1:k+3) -= conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2] [1; nu1; nu2]'.
    //
    // For convenience, since H(k+4,k+1:k+2)=0, we start by computing 
    //
    //   innerProd = conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2]
    //             = conj(tau) H(k+4,k+3) nu2
    //
    // TODO(poulson): Avoid memory allocations for several of the temporaries
    // below in the case where we have heap scalars (e.g., BigFloat).
    const Int lastRowBulgeEnd =
      Min( firstBulge+numFullBulges, (winEnd-packetBeg-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=firstBulge; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        F* EL_RESTRICT h1 = &HBuf[(bulgeBeg+1)*HLDim]; 
        F* EL_RESTRICT h2 = &HBuf[(bulgeBeg+2)*HLDim];
        F* EL_RESTRICT h3 = &HBuf[(bulgeBeg+3)*HLDim];

        const F* EL_RESTRICT w = &WBuf[bulge*WLDim];
        const F omega0Conj = Conj(w[0]);
        const F omega1Conj = Conj(w[1]);
        const F omega2 = w[2];
        const F omega2Conj = Conj(w[2]);

        const F innerProd = omega0Conj*h3[bulgeBeg+4]*omega2;
        h1[bulgeBeg+4] = -innerProd;
        h2[bulgeBeg+4] = -innerProd*omega1Conj;
        h3[bulgeBeg+4] -= innerProd*omega2Conj;
    }
}

#ifdef EL_REINTERPRET_COMPLEX
// The following performs all computation using the C99 complex classes in
// <complex.h> and seems to yield significantly better performance than the
// El::Complex variant (at least with g++ 5).
template<typename F,
         typename=EnableIf<IsComplexBlasScalar<F>>,
         typename=void>
void ApplyReflectorsOpt
( Matrix<F>& H,
  Int winBeg,
  Int winEnd,
  Int chaseBeg,
  Int packetBeg,
  Int transformRowBeg,
  Int transformColEnd,
  Matrix<F>& Z,
  bool wantSchurVecs,
  Matrix<F>& U,
  const Matrix<F>& W,
  Int firstBulge,
  Int numBulges,
  bool accumulate,
  bool progress )
{
    DEBUG_CSE
    typedef Base<F> Real;
    //typedef std::complex<Real> FastF; using std::conj;
    typedef FastComplex<Real> FastF;

    const Int lastBulge = firstBulge + numBulges - 1;
    const Int lastBulgeBeg = packetBeg + 3*lastBulge;
    const bool haveSmallBulge = ( lastBulgeBeg == winEnd-3 );
    DEBUG_ONLY(
      if( lastBulgeBeg > winEnd-3 )
          LogicError("Last bulge starts too late");
    )
    const Int numFullBulges = ( haveSmallBulge ? numBulges-1 : numBulges );
    const Int clippedChaseBeg = Max(chaseBeg,winBeg-1);

    FastF* EL_RESTRICT HBuf = reinterpret_cast<FastF*>(H.Buffer());
    FastF* EL_RESTRICT ZBuf = reinterpret_cast<FastF*>(Z.Buffer());
    FastF* EL_RESTRICT UBuf = reinterpret_cast<FastF*>(U.Buffer());
    const FastF* EL_RESTRICT WBuf =
      reinterpret_cast<const FastF*>(W.LockedBuffer());
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();
    const Int ULDim = U.LDim();
    const Int WLDim = W.LDim();

    // Apply from the left
    // ===================
    for( Int j=Max(chaseBeg,winBeg); j<transformColEnd; ++j )
    {
        // Avoid re-applying freshly generated reflections
        // (which is only relevant when (j-packetBeg) % 3 = 0)
        const Int applyBulgeEnd =
          Min( firstBulge+numFullBulges, (j-packetBeg+2)/3 );
        FastF* hj = reinterpret_cast<FastF*>(&HBuf[j*HLDim]);
        for( Int bulge=firstBulge; bulge<applyBulgeEnd; ++bulge )
        {
            const Int bulgeBeg = packetBeg + 3*bulge;
            const FastF* EL_RESTRICT w =
              reinterpret_cast<const FastF*>(&WBuf[bulge*WLDim]);
            const FastF innerProd =
              w[0]*(hj[bulgeBeg+1] + conj(w[1])*hj[bulgeBeg+2] +
                    conj(w[2])*hj[bulgeBeg+3]);
            hj[bulgeBeg+1] -= innerProd;
            hj[bulgeBeg+2] -= innerProd*w[1];
            hj[bulgeBeg+3] -= innerProd*w[2];
        }
    }
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const FastF* EL_RESTRICT w =
          reinterpret_cast<const FastF*>(&WBuf[bulge*WLDim]);
        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )
        for( Int j=bulgeBeg+1; j<transformColEnd; ++j )
        {
            FastF* EL_RESTRICT hj = reinterpret_cast<FastF*>(&HBuf[j*HLDim]);
            const FastF innerProd =
              w[0]*(hj[bulgeBeg+1] + conj(w[1])*hj[bulgeBeg+2]);
            hj[bulgeBeg+1] -= innerProd;
            hj[bulgeBeg+2] -= innerProd*w[1];
        }
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    const Int UHeight = U.Height();
    const Int ZHeight = Z.Height();
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const FastF* EL_RESTRICT w =
          reinterpret_cast<const FastF*>(&WBuf[bulge*WLDim]);
        const FastF omega0Conj = conj(w[0]);
        const FastF omega1 = w[1];
        const FastF omega1Conj = conj(w[1]);
        FastF* EL_RESTRICT h1 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+1)*HLDim]);
        FastF* EL_RESTRICT h2 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+2)*HLDim]);
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+3); ++i )
        {
            const FastF innerProd = omega0Conj*(h1[i] + omega1*h2[i]);
            h1[i] -= innerProd;
            h2[i] -= innerProd*omega1Conj;
        }

        if( accumulate )
        {
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            FastF* EL_RESTRICT u1 =
              reinterpret_cast<FastF*>(&UBuf[(bulgeBegRel+1)*ULDim]);
            FastF* EL_RESTRICT u2 =
              reinterpret_cast<FastF*>(&UBuf[(bulgeBegRel+2)*ULDim]);
            for( Int i=0; i<UHeight; ++i ) 
            {
                const FastF innerProd = omega0Conj*(u1[i] + omega1*u2[i]);
                u1[i] -= innerProd;
                u2[i] -= innerProd*omega1Conj;
            }
        }
        else if( wantSchurVecs )
        {
            FastF* EL_RESTRICT z1 =
              reinterpret_cast<FastF*>(&ZBuf[(bulgeBeg+1)*ZLDim]);
            FastF* EL_RESTRICT z2 =
              reinterpret_cast<FastF*>(&ZBuf[(bulgeBeg+2)*ZLDim]);
            for( Int i=0; i<ZHeight; ++i )
            {
                const FastF innerProd = omega0Conj*(z1[i] + omega1*z2[i]);
                z1[i] -= innerProd;
                z2[i] -= innerProd*omega1Conj;
            }
        }
    }
    // Even though we are logically chasing the bulges starting from the
    // bottom-right, it is typically faster to traverse from the top-left
    for( Int bulge=firstBulge; bulge<firstBulge+numFullBulges; ++bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const FastF* EL_RESTRICT w =
          reinterpret_cast<const FastF*>(&WBuf[bulge*WLDim]);
        FastF* EL_RESTRICT h1 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+1)*HLDim]);
        FastF* EL_RESTRICT h2 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+2)*HLDim]);
        FastF* EL_RESTRICT h3 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+3)*HLDim]);
        const FastF omega0Conj = conj(w[0]);
        const FastF omega1 = w[1]; 
        const FastF omega1Conj = conj(w[1]);
        const FastF omega2 = w[2];
        const FastF omega2Conj = conj(w[2]);
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+4); ++i )
        {
            const FastF innerProd =
              omega0Conj*(h1[i] + omega1*h2[i] + omega2*h3[i]);
            h1[i] -= innerProd;
            h2[i] -= innerProd*omega1Conj;
            h3[i] -= innerProd*omega2Conj;
        }

        if( accumulate )
        {
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            FastF* EL_RESTRICT u1 =
              reinterpret_cast<FastF*>(&UBuf[(bulgeBegRel+1)*ULDim]);
            FastF* EL_RESTRICT u2 =
              reinterpret_cast<FastF*>(&UBuf[(bulgeBegRel+2)*ULDim]);
            FastF* EL_RESTRICT u3 =
              reinterpret_cast<FastF*>(&UBuf[(bulgeBegRel+3)*ULDim]);
            for( Int i=0; i<UHeight; ++i ) 
            {
                const FastF innerProd =
                  omega0Conj*(u1[i]+omega1*u2[i]+omega2*u3[i]);
                u1[i] -= innerProd;
                u2[i] -= innerProd*omega1Conj;
                u3[i] -= innerProd*omega2Conj;
            }
        }
        else if( wantSchurVecs )
        {
            FastF* EL_RESTRICT z1 =
              reinterpret_cast<FastF*>(&ZBuf[(bulgeBeg+1)*ZLDim]);
            FastF* EL_RESTRICT z2 =
              reinterpret_cast<FastF*>(&ZBuf[(bulgeBeg+2)*ZLDim]);
            FastF* EL_RESTRICT z3 =
              reinterpret_cast<FastF*>(&ZBuf[(bulgeBeg+3)*ZLDim]);
            for( Int i=0; i<ZHeight; ++i )
            {
                const FastF innerProd =
                  omega0Conj*(z1[i]+omega1*z2[i]+omega2*z3[i]);
                z1[i] -= innerProd;
                z2[i] -= innerProd*omega1Conj;
                z3[i] -= innerProd*omega2Conj;
            }
        }
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    Int firstVigBulge = firstBulge;
    Int numVigBulges = numBulges;
    if( packetBeg + 3*firstBulge == winBeg-1 )
    {
        // The first bulge still needs to be introduced
        ++firstVigBulge;
        --numVigBulges;
    }
    VigilantDeflation
    ( H, winBeg, winEnd, packetBeg, firstVigBulge, numVigBulges, progress );

    // Form the last row of the single-step bulge chase
    //
    //       ~ ~ ~                 ~ ~ ~
    //     -----------          -----------
    //    | B B B B x |        | x x x x x |
    //  ~ | B B B B x |      ~ | x B B B B |
    //  ~ | B B B B x | |->  ~ |   B B B B |.
    //  ~ | B B B B x |      ~ |   B B B B |
    //    |       x x |        |   B B B B |
    //     -----------          -----------
    //
    // The last row is introduced from the transformation
    //
    //  H(k+4,k+1:k+3) -= conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2] [1; nu1; nu2]'.
    //
    // For convenience, since H(k+4,k+1:k+2)=0, we start by computing 
    //
    //   innerProd = conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2]
    //             = conj(tau) H(k+4,k+3) nu2
    //
    // TODO(poulson): Avoid memory allocations for several of the temporaries
    // below in the case where we have heap scalars (e.g., BigFloat).
    const Int lastRowBulgeEnd =
      Min( firstBulge+numFullBulges, (winEnd-packetBeg-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=firstBulge; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        FastF* EL_RESTRICT h1 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+1)*HLDim]);
        FastF* EL_RESTRICT h2 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+2)*HLDim]);
        FastF* EL_RESTRICT h3 =
          reinterpret_cast<FastF*>(&HBuf[(bulgeBeg+3)*HLDim]);
        const FastF* EL_RESTRICT w =
          reinterpret_cast<const FastF*>(&WBuf[bulge*WLDim]);
        const FastF innerProd = conj(w[0])*h3[bulgeBeg+4]*w[2];
        h1[bulgeBeg+4] = -innerProd;
        h2[bulgeBeg+4] = -innerProd*conj(w[1]);
        h3[bulgeBeg+4] -= innerProd*conj(w[2]);
    }
}
#endif // ifdef EL_REINTERPRET_COMPLEX

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP
