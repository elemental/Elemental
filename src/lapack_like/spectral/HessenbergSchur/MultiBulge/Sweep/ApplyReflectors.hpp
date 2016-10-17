/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP

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
  Matrix<F>& W,
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
            const Int UHeight = U.Height();

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
            const Int UHeight = U.Height();
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
        const Int k = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        const F& tau = w[0];
        const F& nu1 = w[1];
        const F& nu2 = w[2];

        const F innerProd = Conj(tau)*H(k+4,k+3)*nu2;
        H(k+4,k+1) = -innerProd;
        H(k+4,k+2) = -innerProd*Conj(nu1);
        H(k+4,k+3) -= innerProd*Conj(nu2);
    }
}

// Unfortunately, it seems to be the case that it is noticeably faster
// for this routine to manually inline the data access than to use the 
// (presumably inlined) Matrix::operator()(int,int) calls.
template<typename F>
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
  Matrix<F>& W,
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

    F* HBuf = H.Buffer();
    F* ZBuf = Z.Buffer();
    F* UBuf = U.Buffer();
    F* WBuf = W.Buffer();
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
            const F* w = &WBuf[bulge*WLDim];
            ApplyLeftReflector
            ( HBuf[(bulgeBeg+1)+j*HLDim],
              HBuf[(bulgeBeg+2)+j*HLDim],
              HBuf[(bulgeBeg+3)+j*HLDim], w );
        }
    }
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &WBuf[bulge*WLDim];
        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )
        for( Int j=bulgeBeg+1; j<transformColEnd; ++j )
            ApplyLeftReflector
            ( HBuf[(bulgeBeg+1)+j*HLDim],
              HBuf[(bulgeBeg+2)+j*HLDim], w );
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    const Int ZHeight = Z.Height();
    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &WBuf[bulge*WLDim];
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+3); ++i )
            ApplyRightReflector
            ( HBuf[i+(bulgeBeg+1)*HLDim],
              HBuf[i+(bulgeBeg+2)*HLDim], w );

        if( accumulate )
        {
            const Int UHeight = U.Height();
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            for( Int i=0; i<UHeight; ++i ) 
                ApplyRightReflector
                ( UBuf[i+(bulgeBegRel+1)*ULDim],
                  UBuf[i+(bulgeBegRel+2)*ULDim], w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<ZHeight; ++i )
                ApplyRightReflector
                ( ZBuf[i+(bulgeBeg+1)*ZLDim],
                  ZBuf[i+(bulgeBeg+2)*ZLDim], w );
        }
    }
    for( Int bulge=firstBulge+numFullBulges-1; bulge>=firstBulge; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        if( bulge >= W.Width() )
            LogicError
            ("W.Width()=",W.Width(),", bulge=",bulge,", firstBulge=",firstBulge,
             ", numFullBulges=",numFullBulges);
        const F* w = &WBuf[bulge*WLDim];
        for( Int i=transformRowBeg; i<Min(winEnd,bulgeBeg+4); ++i )
            ApplyRightReflector
            ( HBuf[i+(bulgeBeg+1)*HLDim],
              HBuf[i+(bulgeBeg+2)*HLDim],
              HBuf[i+(bulgeBeg+3)*HLDim], w );

        if( accumulate )
        {
            const Int UHeight = U.Height();
            const Int bulgeBegRel = (bulgeBeg-clippedChaseBeg) - 1;
            if( bulgeBegRel+3 >= UHeight )
                LogicError
                ("UHeight=",UHeight,", bulgeBeg=",bulgeBeg,
                 ", chaseBeg=",chaseBeg,", winBeg=",winBeg,
                 ", bulgeBegRel=",bulgeBegRel,", firstBulge=",firstBulge,
                 ", numFullBulges=",numFullBulges);
            for( Int i=0; i<UHeight; ++i ) 
                ApplyRightReflector
                ( UBuf[i+(bulgeBegRel+1)*ULDim],
                  UBuf[i+(bulgeBegRel+2)*ULDim],
                  UBuf[i+(bulgeBegRel+3)*ULDim], w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<ZHeight; ++i )
                ApplyRightReflector
                ( ZBuf[i+(bulgeBeg+1)*ZLDim],
                  ZBuf[i+(bulgeBeg+2)*ZLDim],
                  ZBuf[i+(bulgeBeg+3)*ZLDim], w );
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
        const Int k = packetBeg + 3*bulge;
        const F* w = &WBuf[bulge*WLDim];
        const F& tau = w[0];
        const F& nu1 = w[1];
        const F& nu2 = w[2];

        const F innerProd = Conj(tau)*HBuf[(k+4)+(k+3)*HLDim]*nu2;
        HBuf[(k+4)+(k+1)*HLDim] = -innerProd;
        HBuf[(k+4)+(k+2)*HLDim] = -innerProd*Conj(nu1);
        HBuf[(k+4)+(k+3)*HLDim] -= innerProd*Conj(nu2);
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP
