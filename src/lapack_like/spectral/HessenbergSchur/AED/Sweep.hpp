/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_AED_SWEEP_HPP
#define EL_SCHUR_HESS_AED_SWEEP_HPP

namespace El {
namespace hess_schur {
namespace aed {

template<typename Real>
void ImplicitQQuadraticSeed
( const Matrix<Real>& H,
  const Complex<Real>& shift0, 
  const Complex<Real>& shift1,
        Real* v )
{
    DEBUG_CSE
    const Real zero(0);
    const Int n = H.Height();
    DEBUG_ONLY(
      if( n != 2 && n != 3 )
          LogicError("Expected n to be 2 or 3"); 
      const bool bothReal = ( shift0.imag() == zero && shift1.imag() == zero );
      const bool conjugate = ( shift0.imag() == -shift1.imag() );
      if( !bothReal && !conjugate )
          LogicError("Assumed shifts were either both real or conjugates");
    )
    if( n == 2 )
    {
        const Real& eta00 = H(0,0);
        const Real& eta01 = H(0,1);
        const Real& eta10 = H(1,0);
        const Real& eta11 = H(1,1);

        // It seems arbitrary whether the scale is computed relative
        // to shift0 or shift1, but we follow LAPACK's convention.
        // (While the choice is irrelevant for conjugate shifts, it is not for
        //  real shifts)
        const Real scale = OneAbs(eta00-shift1) + Abs(eta10);
        if( scale == zero )
        {
            v[0] = v[1] = zero;
        }
        else
        {
            // Normalize the first column by the scale
            Real eta10Scale = eta10 / scale;
            v[0] = eta10Scale*eta01 +
                   (eta00-shift0.real())*((eta00-shift1.real())/scale) -
                   shift0.imag()*(shift1.imag()/scale);
            v[1] = eta10Scale*(eta00+eta11-shift0.real()-shift1.real());
        }
    }
    else
    {
        const Real& eta00 = H(0,0);
        const Real& eta01 = H(0,1);
        const Real& eta02 = H(0,2);
        const Real& eta10 = H(1,0);
        const Real& eta11 = H(1,1);
        const Real& eta12 = H(1,2);
        const Real& eta20 = H(2,0);
        const Real& eta21 = H(2,1); 
        const Real& eta22 = H(2,2);

        const Real scale = OneAbs(eta00-shift1) + Abs(eta10) + Abs(eta20);
        if( scale == zero )
        {
            v[0] = v[1] = v[2] = 0;
        }
        else
        {
            // Normalize the first column by the scale
            const Real eta10Scale = eta10 / scale;
            const Real eta20Scale = eta20 / scale;
            v[0] = (eta00-shift0.real())*((eta00-shift1.real())/scale) -
                   shift0.imag()*(shift1.imag()/scale) + eta01*eta10Scale +
                   eta02*eta20Scale;
            v[1] = eta10Scale*(eta00+eta11-shift0.real()-shift1.real()) +
                   eta12*eta20Scale;
            v[2] = eta20Scale*(eta00+eta22-shift0.real()-shift1.real()) +
                   eta21*eta10Scale;
        }
    }
}

template<typename F>
void ImplicitQQuadraticSeed
( const Matrix<F>& H,
  const F& shift0,
  const F& shift1,
        F* v )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Real zero(0);
    const Int n = H.Height();
    DEBUG_ONLY(
      if( n != 2 && n != 3 )
          LogicError("Expected n to be 2 or 3"); 
    )
    if( n == 2 )
    {
        const F& eta00 = H(0,0);
        const F& eta01 = H(0,1);
        const F& eta10 = H(1,0);
        const F& eta11 = H(1,1);

        // It seems arbitrary whether the scale is computed relative
        // to shift0 or shift1, but we follow LAPACK's convention.
        // (While the choice is irrelevant for conjugate shifts, it is not for
        //  real shifts)
        const Real scale = OneAbs(eta00-shift1) + OneAbs(eta10);
        if( scale == zero )
        {
            v[0] = v[1] = zero;
        }
        else
        {
            // Normalize the first column by the scale
            F eta10Scale = eta10 / scale;
            v[0] = eta10Scale*eta01 + (eta00-shift0)*((eta00-shift1)/scale);
            v[1] = eta10Scale*(eta00+eta11-shift0-shift1);
        }
    }
    else
    {
        const F& eta00 = H(0,0);
        const F& eta01 = H(0,1);
        const F& eta02 = H(0,2);
        const F& eta10 = H(1,0);
        const F& eta11 = H(1,1);
        const F& eta12 = H(1,2);
        const F& eta20 = H(2,0);
        const F& eta21 = H(2,1); 
        const F& eta22 = H(2,2);

        const Real scale = OneAbs(eta00-shift1) + OneAbs(eta10) + OneAbs(eta20);
        if( scale == zero )
        {
            v[0] = v[1] = v[2] = 0;
        }
        else
        {
            // Normalize the first column by the scale
            const F eta10Scale = eta10 / scale;
            const F eta20Scale = eta20 / scale;
            v[0] = (eta00-shift0)*((eta00-shift1)/scale) +
                   eta01*eta10Scale + eta02*eta20Scale;
            v[1] = eta10Scale*(eta00+eta11-shift0-shift1) + eta12*eta20Scale;
            v[2] = eta20Scale*(eta00+eta22-shift0-shift1) + eta21*eta10Scale;
        }
    }
}

template<typename F>
void IntroduceBulge
( const Matrix<F>& H,
  const Complex<Base<F>>& shift0,
  const Complex<Base<F>>& shift1,
        F* v )
{
    DEBUG_CSE
    const Int n = H.Height();
    ImplicitQQuadraticSeed( H, shift0, shift1, v ); 
    F beta = v[0];
    v[0] = lapack::Reflector( n, beta, &v[1], 1 );
}

// NOTE: This should only be called for real matrices, where one can assume
// conjugate pairs of shifts
template<typename Real>
void PairShifts( Matrix<Complex<Real>>& shifts )
{
    DEBUG_CSE
    const Int numShifts = shifts.Height();

    Complex<Real> tmp;
    for( Int i=numShifts-1; i>=2; i-=2 )
    {
        if( shifts(i).imag() != -shifts(i-1).imag() )
        {
            tmp = shifts(i);
            shifts(i) = shifts(i-1);
            shifts(i-1) = shifts(i-2);
            shifts(i-2) = tmp;
        }
    }

    DEBUG_ONLY(
      for( Int i=numShifts-1; i>=2; i-=2 )
      {
          if( shifts(i).imag() != -shifts(i-1).imag() )
          {
              RuntimeError("Shifts were not properly paired");
          }
      }
    )
}

template<typename F>
void ComputeReflectors
( Matrix<F>& H,
  Int winBeg,
  Matrix<Complex<Base<F>>>& shifts,
  Matrix<F>& W,
  Int packetBeg,
  Int fullBeg,
  Int fullEnd,
  bool have3x3 )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Real realZero(0);
    const F zero(0);
    const Real ulp = limits::Precision<Real>();

    // Set aside space for a Householder vector for a candidate reinflation of
    // a deflated bulge
    vector<F> wCand(3);

    if( have3x3 )
    {
        const Int bulge = fullEnd;
        const Complex<Real> shift0 = shifts(2*bulge);
        const Complex<Real> shift1 = shifts(2*bulge+1);
        const Int bulgeBeg = packetBeg + 3*bulge;
        F* w = &W(0,bulge);

        if( bulgeBeg == winBeg-1 )
        {
            auto ind1 = IR(bulgeBeg+1,bulgeBeg+3);
            auto H11BR = H(ind1,ind1);
            IntroduceBulge( H11BR, shift0, shift1, w );
        }
        else
        {
            // Find the reflection for chasing the 3x3 bulge
            F beta = H( bulgeBeg+1, bulgeBeg );
            w[1] = H( bulgeBeg+2, bulgeBeg );
            w[0] = lapack::Reflector( 2, beta, &w[1], 1 );
            H( bulgeBeg+1, bulgeBeg ) = beta;
            H( bulgeBeg+2, bulgeBeg ) = realZero;
        }
    }
    for( Int bulge=fullEnd-1; bulge>=fullBeg; --bulge )
    {
        const Complex<Real> shift0 = shifts(2*bulge);
        const Complex<Real> shift1 = shifts(2*bulge+1);
        const Int bulgeBeg = packetBeg + 3*bulge;
        F* w = &W(0,bulge);
        auto ind1 = IR(bulgeBeg+1,bulgeBeg+4);
        auto H11BR = H(ind1,ind1);

        if( bulgeBeg == winBeg-1 )
        {
            IntroduceBulge( H11BR, shift0, shift1, w );
        }
        else
        {
            // Prepare to chase the bulge down a step
            F& eta10 = H(bulgeBeg+1,bulgeBeg);
            F& eta20 = H(bulgeBeg+2,bulgeBeg);
            F& eta30 = H(bulgeBeg+3,bulgeBeg);

            F beta = eta10;
            w[1] = eta20;
            w[2] = eta30;
            w[0] = lapack::Reflector( 3, beta, &w[1], 1 );
                    
            // "Vigilantly" search for a deflation
            // (deflation within the interior is exceedingly rare,
            //  but this check should be essentially free)
            const F& eta31 = H(bulgeBeg+3,bulgeBeg+1);
            const F& eta32 = H(bulgeBeg+3,bulgeBeg+2);
            if( eta30 != zero || eta31 != zero || eta32 == zero )
            {
                // Bulge has not collapsed
                eta10 = beta;
                eta20 = realZero;
                eta30 = realZero;
            }
            else
            {
                // Bulge has collapsed
                IntroduceBulge( H11BR, shift0, shift1, wCand.data() );
                F innerProd = wCand[0]*(eta10+Conj(wCand[1])*eta20);

                const F& eta00 = H(bulgeBeg,  bulgeBeg  );
                const F& eta11 = H(bulgeBeg+1,bulgeBeg+1);
                const F& eta22 = H(bulgeBeg+2,bulgeBeg+2);
                if( OneAbs(eta20-wCand[1]*innerProd) +
                    OneAbs(wCand[2]*innerProd) >=
                    ulp*(OneAbs(eta00)+OneAbs(eta11)+OneAbs(eta22)) )
                {
                    // The proposed bulge was unacceptable;
                    // continue using the collapsed one with regret
                    eta10 = beta;
                    eta20 = realZero;
                    eta30 = realZero;
                }
                else
                {
                    // Accept the proposed replacement bulge
                    eta10 -= innerProd;
                    eta20 = realZero; 
                    eta30 = realZero;
                    w[0] = wCand[0];
                    w[1] = wCand[1];
                    w[2] = wCand[2];
                }
            }
        }
    }
}

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
        {
            ApplyLeftReflector( H(bulgeBeg+1,j), H(bulgeBeg+2,j), w );
        }
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
        {
            ApplyRightReflector( H(i,bulgeBeg+1), H(i,bulgeBeg+2), w );
        }

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
            {
                ApplyRightReflector( U(i,kU), U(i,kU+1), w );
            }
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
            {
                ApplyRightReflector( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), w );
            }
        }
    }
    for( Int bulge=fullEnd-1; bulge>=fullBeg; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);

        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+4); ++i )
        {
            ApplyRightReflector
            ( H(i,bulgeBeg+1), H(i,bulgeBeg+2), H(i,bulgeBeg+3), w );
        }

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
            {
                ApplyRightReflector( U(i,kU), U(i,kU+1), U(i,kU+2), w );
            }
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
            {
                ApplyRightReflector
                ( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), Z(i,bulgeBeg+3), w );
            }
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
        {
            ApplyLeftReflector
            ( HBuf[(bulgeBeg+1)+j*HLDim],
              HBuf[(bulgeBeg+2)+j*HLDim], w );
        }
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
        {
            ApplyRightReflector
            ( HBuf[i+(bulgeBeg+1)*HLDim],
              HBuf[i+(bulgeBeg+2)*HLDim], w );
        }

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
            {
                ApplyRightReflector
                ( UBuf[i+ kU   *ULDim],
                  UBuf[i+(kU+1)*ULDim], w );
            }
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
            {
                ApplyRightReflector
                ( ZBuf[i+(bulgeBeg+1)*ZLDim],
                  ZBuf[i+(bulgeBeg+2)*ZLDim], w );
            }
        }
    }
    for( Int bulge=fullEnd-1; bulge>=fullBeg; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const F* w = &W(0,bulge);

        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+4); ++i )
        {
            ApplyRightReflector
            ( HBuf[i+(bulgeBeg+1)*HLDim],
              HBuf[i+(bulgeBeg+2)*HLDim],
              HBuf[i+(bulgeBeg+3)*HLDim], w );
        }

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
            {
                ApplyRightReflector
                ( UBuf[i+ kU   *ULDim],
                  UBuf[i+(kU+1)*ULDim],
                  UBuf[i+(kU+2)*ULDim], w );
            }
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
            {
                ApplyRightReflector
                ( ZBuf[i+(bulgeBeg+1)*ZLDim],
                  ZBuf[i+(bulgeBeg+2)*ZLDim],
                  ZBuf[i+(bulgeBeg+3)*ZLDim], w );
            }
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

template<typename F>
void Sweep
( Matrix<F>& H,
  Matrix<Complex<Base<F>>>& shifts,
  Matrix<F>& Z,
  Matrix<F>& U,
  Matrix<F>& W,
  Matrix<F>& WAccum,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Real realZero(0);
    const Int n = H.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );

    const Int numShifts = shifts.Height();
    DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts..."); 
      if( numShifts % 2 != 0 )
          LogicError("Expected an even number of sweeps");
    )
    const Int numBulges = numShifts / 2;

    if( !IsComplex<F>::value )
        PairShifts( shifts );

    // TODO: Decide if this is strictly necessary
    H(winBeg+2,winBeg) = realZero;

    // Set aside space for storing either the three nonzero entries of the first
    // column of a quadratic polynomial for introducing each bulge or the scalar
    // 'tau' and non-unit subvector v such that the 2x2 or 3x3 Householder
    // reflection I - tau*[1;v] [1;v]' can safely be stored within a vector of
    // length 3 (following the LAPACK convention that 'tau' will lie in the
    // first entry and 'v' will begin in the second entry).
    W.Resize( 3, numBulges );

    // Each time a new 4x4 bulge is introduced into the upper-left corner
    // of the matrix, we think of the bulge as having started at index
    // (winBeg-1), as, after the introduction of the 3x3 Householder
    // similarity using the implicit Q theorem on (H-sigma_0 I)*(H-sigma_1 I),
    // the bulge has starting index winBeg.
    //
    // Initialize with the last bulge about to be introduced in the upper-left
    const Int ghostBeg = (winBeg-1) - 3*(numBulges-1);

    // The last bulge must be at least 3x3 in order to involve a 2x2 
    // reflection, so it must start before winEnd-2
    const Int ghostEnd = winEnd-2;

    // Each movement of a packet involves shifting the first bulge one 
    // column right of the starting position of the last bulge. This involves
    // shifting every bulge right 3*(numBulges-1) + 1 columns.
    const Int ghostStride = 3*(numBulges-1) + 1;

    // The total effected region of each movement is therefore (at most)
    // the span of numBulges 3x3 Householder reflections and the translation
    // distance of the packet (ghostStride)
    const Int slabSize = 3*numBulges + ghostStride;

    for( Int ghostCol=ghostBeg; ghostCol<ghostEnd; ghostCol+=ghostStride )
    {
        // Note that this slab endpoint may be past winEnd
        const Int slabEnd = ghostCol + slabSize;
        if( ctrl.accumulateReflections )
            Identity( U, slabSize-1, slabSize-1 );

        const Int packetEnd = Min(ghostCol+ghostStride,ghostEnd);
        for( Int packetBeg=ghostCol; packetBeg<packetEnd; ++packetBeg )
        {
            const Int fullBeg = Max( 0, ((winBeg-1)-packetBeg+2)/3 );
            const Int fullEnd = Min( numBulges, (winEnd-packetBeg-1)/3 );

            // 2x2 reflectors can only occur if a bulge occupies the 3x3 in the
            // bottom-right corner (which begins at index winEnd-3)
            const bool have3x3 =
              ( fullEnd < numBulges && packetBeg+3*fullEnd == winEnd-3 );

            ComputeReflectors
            ( H, winBeg, shifts, W, packetBeg, fullBeg, fullEnd, have3x3 );

            Int transformBeg;
            if( ctrl.accumulateReflections )
                transformBeg = Max( winBeg, ghostCol );
            else if( ctrl.fullTriangle )
                transformBeg = 0;
            else
                transformBeg = winBeg;

            Int transformEnd;
            if( ctrl.accumulateReflections )
                transformEnd = Min( slabEnd, winEnd );
            else if( ctrl.fullTriangle )
                transformEnd = n;
            else
                transformEnd = winEnd;

            ApplyReflectorsOpt
            ( H, winBeg, winEnd,
              slabSize, ghostCol, packetBeg, transformBeg, transformEnd,
              Z, ctrl.wantSchurVecs, U, W,
              fullBeg, fullEnd, have3x3, ctrl.accumulateReflections );
        }

        if( ctrl.accumulateReflections )
        {
            const Int transformBeg = ( ctrl.fullTriangle ? 0 : winBeg );
            const Int transformEnd = ( ctrl.fullTriangle ? n : winEnd );
 
            const Int slabRelBeg = Max(0,(winBeg-1)-ghostCol);
            const Int nU = (slabSize-1) - Max(0,slabEnd-winEnd) - slabRelBeg;

            auto contractInd = IR(0,nU) + slabRelBeg;
            auto UAccum = U( contractInd, contractInd );

            // Horizontal far-from-diagonal application
            const Int rightIndBeg = Min(ghostCol+slabSize,winEnd);
            const Int rightIndEnd = transformEnd;
            const auto rightInd = IR(rightIndBeg,rightIndEnd);
            auto horzInd = IR(0,nU) + (ghostCol+slabRelBeg+1);
            auto HHorzFar = H( horzInd, rightInd );
            Gemm( ADJOINT, NORMAL, F(1), UAccum, HHorzFar, WAccum );
            HHorzFar = WAccum;

            // Vertical far-from-diagonal application
            auto vertInd = IR(transformBeg,Max(winBeg,ghostCol));
            auto HVertFar = H( vertInd, horzInd );
            Gemm( NORMAL, NORMAL, F(1), HVertFar, UAccum, WAccum );
            HVertFar = WAccum;

            if( ctrl.wantSchurVecs )
            {
                auto ZSub = Z( ALL, horzInd );
                Gemm( NORMAL, NORMAL, F(1), ZSub, UAccum, WAccum );
                ZSub = WAccum;
            }
        }
    }
}

} // namespace aed
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_AED_SWEEP_HPP
