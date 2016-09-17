/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_APPLY_REFLECTORS_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

// TODO(poulson): Avoid temporaries for innerProd and innerProd*nu1 to void
// memory allocations for heap scalars?
template<typename F>
void ApplyLeftReflector( F& eta0, F& eta1, const F* w )
{
    // Update
    //
    //   | eta0 | -= tau |  1  | | 1, conj(nu1) | | eta0 |
    //   | eta1 |        | nu1 |                  | eta1 |
    //
    // where tau is stored in w[0] and nu1 in w[1].
    //
    const F& tau = w[0];
    const F& nu1 = w[1];

    const F innerProd = tau*(eta0+Conj(nu1)*eta1);
    eta0 -= innerProd;
    eta1 -= innerProd*nu1;
}

template<typename F>
void ApplyRightReflector( F& eta0, F& eta1, const F* w )
{
    eta0 = Conj(eta0);
    eta1 = Conj(eta1);
    ApplyLeftReflector( eta0, eta1, w );
    eta0 = Conj(eta0);
    eta1 = Conj(eta1);
}

// TODO(poulson): Avoid temporaries for innerProd, innerProd*nu1, and
// innerProd*nu2 to avoid memory allocations for heap scalars?
template<typename F>
void ApplyLeftReflector( F& eta0, F& eta1, F& eta2, const F* w )
{
    // Update
    //
    //   | eta0 | -= tau |  1  | | 1, conj(nu1), conj(nu2) | | eta0 |
    //   | eta1 |        | nu1 |                             | eta1 |
    //   | eta2 |        | nu2 |                             | eta2 |
    //
    // where tau is stored in w[0], nu1 in w[1], and nu2 in w[2].
    //
    const F& tau = w[0]; 
    const F& nu1 = w[1];
    const F& nu2 = w[2];

    const F innerProd = tau*(eta0+Conj(nu1)*eta1+Conj(nu2)*eta2);
    eta0 -= innerProd;
    eta1 -= innerProd*nu1;
    eta2 -= innerProd*nu2;
}

template<typename F>
void ApplyRightReflector( F& eta0, F& eta1, F& eta2, const F* w )
{
    eta0 = Conj(eta0);
    eta1 = Conj(eta1);
    eta2 = Conj(eta2);
    ApplyLeftReflector( eta0, eta1, eta2, w );
    eta0 = Conj(eta0);
    eta1 = Conj(eta1);
    eta2 = Conj(eta2);
}

// To be performed during the application of the reflections from the right
template<typename F>
void VigilantDeflation
( Matrix<F>& H,
  Int winBeg,
  Int winEnd,
  Int packetBeg,
  Int vigBeg,
  Int vigEnd )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Real zero(0);
    const F complexZero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(H.Height())/ulp);

    // NOTE: LAPACK has a strange increment of vigEnd by one if
    //       packetBeg is equal to winBeg-3 here...
    for( Int bulge=vigEnd-1; bulge>=vigBeg; --bulge )
    {
        const Int k = Min( winEnd-2, packetBeg+3*bulge );
        F& eta00 = H(k  ,k  );
        F& eta01 = H(k  ,k+1);
        F& eta10 = H(k+1,k  );
        F& eta11 = H(k+1,k+1);
        const Real eta00Abs = OneAbs(eta00);
        const Real eta10Abs = OneAbs(eta10);
        const Real eta11Abs = OneAbs(eta11);

        if( eta10 == complexZero )
            continue;
        Real localScale = eta00Abs + eta11Abs;
        if( localScale == zero )
        {
            if( k >= winBeg+3 )
                localScale += OneAbs(H(k,k-3));
            if( k >= winBeg+2 )
                localScale += OneAbs(H(k,k-2));
            if( k >= winBeg+1 )
                localScale += OneAbs(H(k,k-1));
            if( k < winEnd-2 )
                localScale += OneAbs(H(k+2,k+1));
            if( k < winEnd-3 )
                localScale += OneAbs(H(k+3,k+1));
            if( k < winEnd-4 )
                localScale += OneAbs(H(k+4,k+1));
        }
        if( eta10Abs <= Max(smallNum,ulp*localScale) )
        {
            const Real eta01Abs = OneAbs(eta01);
            const Real diagDiffAbs = OneAbs(eta00-eta11);
            Real offMax = Max( eta10Abs, eta01Abs );
            Real offMin = Min( eta10Abs, eta01Abs );
            Real diagMax = Max( eta11Abs, diagDiffAbs );
            Real diagMin = Min( eta11Abs, diagDiffAbs );
            Real scale = diagMax + offMax;
            Real localScale2 = diagMin*(diagMax/scale);
            if( localScale2 == zero ||
                offMin*(offMax/scale) <= Max(smallNum,ulp*localScale2) )
            {
                eta10 = zero;
            }
        }
    }
}

template<typename F>
void ApplyReflectors
( Matrix<F>& H,
  Int winBeg,
  Int winEnd,
  Int slabSize,
  Int ghostCol,
  Int packetBeg,
  Int transformBeg,
  Int transformEnd,
  Matrix<F>& Z,
  bool wantSchurVecs,
  Matrix<F>& U,
  Matrix<F>& W,
  Int fullBeg,
  Int fullEnd,
  bool have3x3,
  bool accumulate )
{
    DEBUG_CSE

    // Apply from the left
    // ===================
    for( Int j=Max(ghostCol,winBeg); j<transformEnd; ++j )
    {
        // Avoid re-applying freshly generated reflections
        // (which is only relevant when (j-packetBeg) % 3 = 0)
        const Int applyBulgeEnd = Min( fullEnd, (j-packetBeg+2)/3 );
        for( Int bulge=fullBeg; bulge<applyBulgeEnd; ++bulge )
        {
            const Int bulgeBeg = packetBeg + 3*bulge;
            const F* w = &W(0,bulge);
            ApplyLeftReflector
            ( H(bulgeBeg+1,j), H(bulgeBeg+2,j), H(bulgeBeg+3,j), w );
        }
    }
    if( have3x3 )
    {
        const Int bulge = fullEnd;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )
        for( Int j=bulgeBeg+1; j<transformEnd; ++j )
            ApplyLeftReflector( H(bulgeBeg+1,j), H(bulgeBeg+2,j), w );
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    const Int nZ = Z.Height();
    // The first relative index of the slab that is in the window
    // (with 4x4 bulges being introduced at winBeg-1)
    const Int slabRelBeg = Max(0,(winBeg-1)-ghostCol);
    if( have3x3 )
    {
        const Int bulge = fullEnd;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+3); ++i )
            ApplyRightReflector( H(i,bulgeBeg+1), H(i,bulgeBeg+2), w );

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
                ApplyRightReflector( U(i,kU), U(i,kU+1), w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
                ApplyRightReflector( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), w );
        }
    }
    for( Int bulge=fullEnd-1; bulge>=fullBeg; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+4); ++i )
            ApplyRightReflector
            ( H(i,bulgeBeg+1), H(i,bulgeBeg+2), H(i,bulgeBeg+3), w );

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
                ApplyRightReflector( U(i,kU), U(i,kU+1), U(i,kU+2), w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
                ApplyRightReflector
                ( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), Z(i,bulgeBeg+3), w );
        }
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    const Int vigBeg =
      ( packetBeg+3*fullBeg == winBeg-1 ? fullBeg+1 : fullBeg );
    const Int vigEnd = ( have3x3 ? fullEnd+1 : fullEnd );
    VigilantDeflation( H, winBeg, winEnd, packetBeg, vigBeg, vigEnd );

    // Form the last row of the result of applying from the right:
    //
    // | X X X X X |     | X X X X X |
    // | X X X X X |     | X X X X X |
    // | X X X X X | |-> |   X X X X |
    // | X X X X X |     |   X X X X |
    // |       X X |     |   X X X X |
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
    const Int lastRowBulgeEnd = Min( fullEnd, (winEnd-packetBeg-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=fullBeg; --bulge )
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
  Int slabSize,
  Int ghostCol,
  Int packetBeg,
  Int transformBeg,
  Int transformEnd,
  Matrix<F>& Z,
  bool wantSchurVecs,
  Matrix<F>& U,
  Matrix<F>& W,
  Int fullBeg,
  Int fullEnd,
  bool have3x3,
  bool accumulate )
{
    DEBUG_CSE
    F* HBuf = H.Buffer();
    F* ZBuf = Z.Buffer();
    F* UBuf = U.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();
    const Int ULDim = U.LDim();

    // Apply from the left
    // ===================
    for( Int j=Max(ghostCol,winBeg); j<transformEnd; ++j )
    {
        // Avoid re-applying freshly generated reflections
        // (which is only relevant when (j-packetBeg) % 3 = 0)
        const Int applyBulgeEnd = Min( fullEnd, (j-packetBeg+2)/3 );
        for( Int bulge=fullBeg; bulge<applyBulgeEnd; ++bulge )
        {
            const Int bulgeBeg = packetBeg + 3*bulge;
            const F* w = &W(0,bulge);
            ApplyLeftReflector
            ( HBuf[(bulgeBeg+1)+j*HLDim],
              HBuf[(bulgeBeg+2)+j*HLDim],
              HBuf[(bulgeBeg+3)+j*HLDim], w );
        }
    }
    if( have3x3 )
    {
        const Int bulge = fullEnd;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )
        for( Int j=bulgeBeg+1; j<transformEnd; ++j )
            ApplyLeftReflector
            ( HBuf[(bulgeBeg+1)+j*HLDim],
              HBuf[(bulgeBeg+2)+j*HLDim], w );
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    const Int nZ = Z.Height();
    // The first relative index of the slab that is in the window
    // (with 4x4 bulges being introduced at winBeg-1)
    const Int slabRelBeg = Max(0,(winBeg-1)-ghostCol);
    if( have3x3 )
    {
        const Int bulge = fullEnd;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+3); ++i )
            ApplyRightReflector
            ( HBuf[i+(bulgeBeg+1)*HLDim],
              HBuf[i+(bulgeBeg+2)*HLDim], w );

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
                ApplyRightReflector
                ( UBuf[i+ kU   *ULDim],
                  UBuf[i+(kU+1)*ULDim], w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
                ApplyRightReflector
                ( ZBuf[i+(bulgeBeg+1)*ZLDim],
                  ZBuf[i+(bulgeBeg+2)*ZLDim], w );
        }
    }
    for( Int bulge=fullEnd-1; bulge>=fullBeg; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+4); ++i )
            ApplyRightReflector
            ( HBuf[i+(bulgeBeg+1)*HLDim],
              HBuf[i+(bulgeBeg+2)*HLDim],
              HBuf[i+(bulgeBeg+3)*HLDim], w );

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
                ApplyRightReflector
                ( UBuf[i+ kU   *ULDim],
                  UBuf[i+(kU+1)*ULDim],
                  UBuf[i+(kU+2)*ULDim], w );
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
                ApplyRightReflector
                ( ZBuf[i+(bulgeBeg+1)*ZLDim],
                  ZBuf[i+(bulgeBeg+2)*ZLDim],
                  ZBuf[i+(bulgeBeg+3)*ZLDim], w );
        }
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    const Int vigBeg =
      ( packetBeg+3*fullBeg == winBeg-1 ? fullBeg+1 : fullBeg );
    const Int vigEnd = ( have3x3 ? fullEnd+1 : fullEnd );
    VigilantDeflation( H, winBeg, winEnd, packetBeg, vigBeg, vigEnd );

    // Form the last row of the result of applying from the right:
    //
    // | X X X X X |     | X X X X X |
    // | X X X X X |     | X X X X X |
    // | X X X X X | |-> |   X X X X |
    // | X X X X X |     |   X X X X |
    // |       X X |     |   X X X X |
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
    const Int lastRowBulgeEnd = Min( fullEnd, (winEnd-packetBeg-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=fullBeg; --bulge )
    {
        const Int k = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);
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

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_HPP
