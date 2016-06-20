/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace hess_qr {
namespace small_bulge_sweep {

template<typename Real>
void ImplicitQQuadraticSeed
( const Matrix<Real>& H,
  const Real& shift0Real, const Real& shift0Imag, 
  const Real& shift1Real, const Real& shift1Imag,
        Real* v )
{
    DEBUG_CSE
    const Real zero(0);
    const Int n = H.Height();
    DEBUG_ONLY(
      if( n != 2 && n != 3 )
          LogicError("Expected n to be 2 or 3"); 
      const bool bothReal = ( shift0Imag == zero && shift1Imag == zero );
      const bool conjugate = ( shift0Imag == -shift1Imag );
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
        const Real scale = Abs(eta00-shift1Real) + Abs(shift1Imag) + Abs(eta10);
        if( scale == zero )
        {
            v[0] = v[1] = zero;
        }
        else
        {
            // Normalize the first column by the scale
            Real eta10Scale = eta10 / scale;
            v[0] = eta10Scale*eta01 +
                   (eta00-shift0Real)*((eta00-shift1Real)/scale) -
                   shift0Imag*(shift1Imag/scale);
            v[1] = eta10Scale*(eta00+eta11-shift0Real-shift1Real);
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

        const Real scale =
          Abs(eta00-shift1Real) + Abs(shift1Imag) + Abs(eta10) + Abs(eta20);
        if( scale == zero )
        {
            v[0] = v[1] = v[2] = 0;
        }
        else
        {
            // Normalize the first column by the scale
            const Real eta10Scale = eta10 / scale;
            const Real eta20Scale = eta20 / scale;
            v[0] = (eta00-shift0Real)*((eta00-shift1Real)/scale) -
                   shift0Imag*(shift1Imag/scale) + eta01*eta10Scale +
                   eta02*eta20Scale;
            v[1] = eta10Scale*(eta00+eta11-shift0Real-shift1Real) +
                   eta12*eta20Scale;
            v[2] = eta20Scale*(eta00+eta22-shift0Real-shift1Real) +
                   eta21*eta10Scale;
        }
    }
}

template<typename Real>
void IntroduceBulge
( const Matrix<Real>& H,
  const Real& shift0Real, const Real& shift0Imag, 
  const Real& shift1Real, const Real& shift1Imag,
        Real* v )
{
    DEBUG_CSE
    const Int n = H.Height();
    ImplicitQQuadraticSeed
    ( H,
      shift0Real, shift0Imag,
      shift1Real, shift1Imag,
      v ); 
    Real beta = v[0];
    v[0] = lapack::Reflector( n, beta, &v[1], 1 );
}

template<typename Real>
void PairShifts( Matrix<Real>& wReal, Matrix<Real>& wImag )
{
    DEBUG_CSE
    const Int numShifts = wReal.Height();
    const Real zero(0);

    Print( wReal, "wReal" );
    Print( wImag, "wImag" );
      
    Real tmp;
    for( Int i=numShifts-1; i>=2; i-=2 )
    {
        if( wImag(i) != -wImag(i-1) )
        {
            tmp = wReal(i);
            wReal(i) = wReal(i-1);
            wReal(i-1) = wReal(i-2);
            wReal(i-2) = tmp;

            tmp = wImag(i);
            wImag(i) = wImag(i-1);
            wImag(i-1) = wImag(i-2);
            wImag(i-2) = tmp;
        }
    }

    DEBUG_ONLY(
      for( Int i=numShifts-1; i>=2; i-=2 )
      {
          if( wImag(i) != -wImag(i-1) )
          {
              Print( wReal, "wRealPair" );
              Print( wImag, "wImagPair" );
              RuntimeError("Shifts were not properly paired");
          }
      }
    )
}

template<typename Real>
void ComputeReflectors
( Matrix<Real>& H,
  Int winBeg,
  Matrix<Real>& realShifts,
  Matrix<Real>& imagShifts,
  Matrix<Real>& W,
  Int packetBeg,
  Int fullBeg,
  Int fullEnd,
  bool have3x3 )
{
    DEBUG_CSE
    const Real zero(0);
    const Real ulp = limits::Precision<Real>();

    // Set aside space for a Householder vector for a candidate reinflation of
    // a deflated bulge
    vector<Real> wCand(3);

    if( have3x3 )
    {
        const Int bulge = fullEnd;
        const Real realShift0 = realShifts(2*bulge);
        const Real imagShift0 = imagShifts(2*bulge);
        const Real realShift1 = realShifts(2*bulge+1);
        const Real imagShift1 = imagShifts(2*bulge+1);
        const Int bulgeBeg = packetBeg + 3*bulge;
        Real* w = &W(0,bulge);

        if( bulgeBeg == winBeg-1 )
        {
            IntroduceBulge
            ( H(IR(bulgeBeg+1,bulgeBeg+3),IR(bulgeBeg+1,bulgeBeg+3)),
              realShift0, imagShift0,
              realShift1, imagShift1,
              w ); 
        }
        else
        {
            // Find the reflection for chasing the 3x3 bulge
            Real beta = H( bulgeBeg+1, bulgeBeg );
            w[1] = H( bulgeBeg+2, bulgeBeg );
            w[0] = lapack::Reflector( 2, beta, &w[1], 1 );
            H( bulgeBeg+1, bulgeBeg ) = beta;
            H( bulgeBeg+2, bulgeBeg ) = zero;
        }
    }
    for( Int bulge=fullEnd-1; bulge>=fullBeg; --bulge )
    {
        const Real realShift0 = realShifts(2*bulge);
        const Real imagShift0 = imagShifts(2*bulge);
        const Real realShift1 = realShifts(2*bulge+1);
        const Real imagShift1 = imagShifts(2*bulge+1);
        const Int bulgeBeg = packetBeg + 3*bulge;
        Real* w = &W(0,bulge);
        auto ind1 = IR(bulgeBeg+1,bulgeBeg+4);
        auto H11BR = H(ind1,ind1);

        if( bulgeBeg == winBeg-1 )
        {
            IntroduceBulge
            ( H11BR,
              realShift0, imagShift0,
              realShift1, imagShift1,
              w ); 
        }
        else
        {
            // Prepare to chase the bulge down a step
            Real& eta10 = H(bulgeBeg+1,bulgeBeg);
            Real& eta20 = H(bulgeBeg+2,bulgeBeg);
            Real& eta30 = H(bulgeBeg+3,bulgeBeg);

            Real beta = eta10;
            w[1] = eta20;
            w[2] = eta30;
            w[0] = lapack::Reflector( 3, beta, &w[1], 1 );
                    
            // "Vigilantly" search for a deflation
            // (deflation within the interior is exceedingly rare,
            //  but this check should be essentially free)
            const Real& eta31 = H(bulgeBeg+3,bulgeBeg+1);
            const Real& eta32 = H(bulgeBeg+3,bulgeBeg+2);
            if( eta30 != zero || eta31 != zero || eta32 == zero )
            {
                // Bulge has not collapsed
                eta10 = beta;
                eta20 = zero;
                eta30 = zero;
            }
            else
            {
                // Bulge has collapsed
                IntroduceBulge
                ( H11BR,
                  realShift0, imagShift0,
                  realShift1, imagShift1,
                  wCand.data() ); 
                Real innerProd = wCand[0]*(eta10+wCand[1]*eta20);

                const Real& eta00 = H(bulgeBeg,  bulgeBeg  );
                const Real& eta11 = H(bulgeBeg+1,bulgeBeg+1);
                const Real& eta22 = H(bulgeBeg+2,bulgeBeg+2);
                if( Abs(eta20-innerProd*wCand[1]) + Abs(innerProd*wCand[2]) >=
                    ulp*(Abs(eta00)+Abs(eta11)+Abs(eta22)) )
                {
                    // The proposed bulge was unacceptable;
                    // continue using the collapsed one with regret
                    eta10 = beta;
                    eta20 = zero;
                    eta30 = zero;
                }
                else
                {
                    // Accept the proposed replacement bulge
                    eta10 -= innerProd;
                    eta20 = zero; 
                    eta30 = zero;
                    w[0] = wCand[0];
                    w[1] = wCand[1];
                    w[2] = wCand[2];
                }
            }
        }
    }
}

template<typename Real>
void ApplyReflector( Real& eta1, Real& eta2, const Real* w )
{
    // Apply I - tau*[1;nu0]*[1;nu0]'
    const Real& tau = w[0];
    const Real& nu0 = w[1];

    const Real innerProd = tau*(eta1+eta2*nu0);
    eta1 -= innerProd;
    eta2 -= innerProd*nu0;
}

template<typename Real>
void ApplyReflector( Real& eta1, Real& eta2, Real& eta3, const Real* w )
{
    // Apply I - tau*[1;v]*[1;v]', where v=[nu0;nu1]
    const Real& tau = w[0]; 
    const Real& nu0 = w[1];
    const Real& nu1 = w[2];

    const Real innerProd = tau*(eta1+eta2*nu0+eta3*nu1);
    eta1 -= innerProd;
    eta2 -= innerProd*nu0;
    eta3 -= innerProd*nu1;
}

// To be performed during the application of the reflections from the right
template<typename Real>
void VigilantDeflation
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  Int packetBeg,
  Int vigBeg,
  Int vigEnd )
{
    DEBUG_CSE
    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(H.Height())/ulp);

    // NOTE: LAPACK has a strange increment of vigEnd by one if
    //       packetBeg is equal to winBeg-3 here...
    for( Int bulge=vigEnd-1; bulge>=vigBeg; --bulge )
    {
        const Int k = Min( winEnd-2, packetBeg+3*bulge );
        Real& eta00 = H(k  ,k  );
        Real& eta01 = H(k  ,k+1);
        Real& eta10 = H(k+1,k  );
        Real& eta11 = H(k+1,k+1);
        const Real eta00Abs = Abs(eta00);
        const Real eta10Abs = Abs(eta10);
        const Real eta11Abs = Abs(eta11);

        if( eta10 == zero )
            continue;
        Real localScale = eta00Abs + eta11Abs;
        if( localScale == zero )
        {
            if( k >= winBeg+3 )
                localScale += Abs(H(k,k-3));
            if( k >= winBeg+2 )
                localScale += Abs(H(k,k-2));
            if( k >= winBeg+1 )
                localScale += Abs(H(k,k-1));
            if( k < winEnd-2 )
                localScale += Abs(H(k+2,k+1));
            if( k < winEnd-3 )
                localScale += Abs(H(k+3,k+1));
            if( k < winEnd-4 )
                localScale += Abs(H(k+4,k+1));
        }
        if( eta10Abs <= Max(smallNum,ulp*localScale) )
        {
            const Real eta01Abs = Abs(eta01);
            const Real diagDiffAbs = Abs(eta00-eta11);
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

template<typename Real>
void ApplyReflectors
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  Int slabSize,
  Int ghostCol,
  Int packetBeg,
  Int transformBeg,
  Int transformEnd,
  Matrix<Real>& Z,
  bool wantSchurVecs,
  Matrix<Real>& U,
  Matrix<Real>& W,
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
            const Real* w = &W(0,bulge);
            ApplyReflector
            ( H(bulgeBeg+1,j), H(bulgeBeg+2,j), H(bulgeBeg+3,j), w );
        }
    }
    if( have3x3 )
    {
        const Int bulge = fullEnd;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const Real* w = &W(0,bulge);

        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )

        for( Int j=bulgeBeg+1; j<transformEnd; ++j )
        {
            ApplyReflector( H(bulgeBeg+1,j), H(bulgeBeg+2,j), w );
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
        const Real* w = &W(0,bulge);

        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+3); ++i )
        {
            ApplyReflector( H(i,bulgeBeg+1), H(i,bulgeBeg+2), w );
        }

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
            {
                ApplyReflector( U(i,kU), U(i,kU+1), w );
            }
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
            {
                ApplyReflector( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), w );
            }
        }
    }
    for( Int bulge=fullEnd-1; bulge>=fullBeg; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const Real* w = &W(0,bulge);

        for( Int i=transformBeg; i<Min(winEnd,bulgeBeg+4); ++i )
        {
            ApplyReflector
            ( H(i,bulgeBeg+1), H(i,bulgeBeg+2), H(i,bulgeBeg+3), w );
        }

        if( accumulate )
        {
            const Int kU = bulgeBeg - ghostCol;
            for( Int i=slabRelBeg; i<slabSize-1; ++i ) 
            {
                ApplyReflector( U(i,kU), U(i,kU+1), U(i,kU+2), w );
            }
        }
        else if( wantSchurVecs )
        {
            for( Int i=0; i<nZ; ++i )
            {
                ApplyReflector
                ( Z(i,bulgeBeg+1), Z(i,bulgeBeg+2), Z(i,bulgeBeg+3), w );
            }
        }
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    const Int vigBeg =
      ( packetBeg+3*fullBeg == winBeg-1 ? fullBeg+1 : fullBeg );
    const Int vigEnd = ( have3x3 ? fullEnd+1 : fullEnd );
    small_bulge_sweep::VigilantDeflation
    ( H, winBeg, winEnd, packetBeg, vigBeg, vigEnd );

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
    //   H(k+4,k+1:k+3) -= tau H(k+4,k+1:k+3) [1; nu0; nu1] [1; nu0; nu1]'.
    //
    // For convenience, since H(k+4,k+1:k+2)=0, we start by computing 
    //
    //   innerProd = tau H(k+4,k+1:k+3) [1; nu0; nu1]
    //             = tau H(k+4,k+3) nu1
    //
    const Int lastRowBulgeEnd = Min( fullEnd, (winEnd-packetBeg-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=fullBeg; --bulge )
    {
        const Int k = packetBeg + 3*bulge;
        const Real* w = &W(0,bulge);
        const Real& tau = w[0];
        const Real& nu0 = w[1];
        const Real& nu1 = w[2];

        const Real innerProd = tau*H(k+4,k+3)*nu1;
        H(k+4,k+1) = -innerProd;
        H(k+4,k+2) = -innerProd*nu0;
        H(k+4,k+3) -= innerProd*nu1;
    }
}

} // namespace small_bulge_sweep

template<typename Real>
void SmallBulgeSweep
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  Matrix<Real>& realShifts,
  Matrix<Real>& imagShifts,
  bool fullTriangle,
  Matrix<Real>& Z,
  bool wantSchurVecs,
  Matrix<Real>& U,
  Matrix<Real>& W,
  Matrix<Real>& WAccum,
  bool accumulate=true )
{
    DEBUG_CSE
    const Real zero(0);
    const Int n = H.Height();
    const Int nZ = Z.Height();

    const Int numShifts = realShifts.Height();
    DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts..."); 
      if( numShifts % 2 != 0 )
          LogicError("Expected an even number of sweeps");
    )
    const Int numBulges = numShifts / 2;

    small_bulge_sweep::PairShifts( realShifts, imagShifts );

    // TODO: Decide if this is strictly necessary
    H(winBeg+2,winBeg) = zero;

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
        if( accumulate )
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

            small_bulge_sweep::ComputeReflectors
            ( H, winBeg, realShifts, imagShifts, W,
              packetBeg, fullBeg, fullEnd, have3x3 );

            Int transformBeg;
            if( accumulate )
                transformBeg = Max( winBeg, ghostCol );
            else if( fullTriangle )
                transformBeg = 0;
            else
                transformBeg = winBeg;

            Int transformEnd;
            if( accumulate )
                transformEnd = Min( slabEnd, winEnd );
            else if( fullTriangle )
                transformEnd = n;
            else
                transformEnd = winEnd;

            small_bulge_sweep::ApplyReflectors
            ( H, winBeg, winEnd,
              slabSize, ghostCol, packetBeg, transformBeg, transformEnd,
              Z, wantSchurVecs, U, W,
              fullBeg, fullEnd, have3x3, accumulate );
        }
        if( accumulate )
        {
            const Int transformBeg = ( fullTriangle ? 0 : winBeg );
            const Int transformEnd = ( fullTriangle ? n : winEnd );
 
            const Int slabRelBeg = Max(0,(winBeg-1)-ghostCol);
            const Int slabWinEnd = Min(slabEnd,winEnd);

            const Int nU = (slabSize-1) - Max(0,slabEnd-winEnd) - slabRelBeg;

            auto contractInd = IR(0,nU) + slabRelBeg;
            auto UAccum = U( contractInd, contractInd );

            // Horizontal far-from-diagonal application
            const Int rightIndBeg = Min(ghostCol+slabSize,winEnd);
            const Int rightIndEnd = transformEnd;
            const auto rightInd = IR(rightIndBeg,rightIndEnd);
            auto horzInd = IR(0,nU) + (ghostCol+slabRelBeg+1);
            auto HHorzFar = H( horzInd, rightInd );
            Gemm( ADJOINT, NORMAL, Real(1), UAccum, HHorzFar, WAccum );
            HHorzFar = WAccum;

            // Vertical far-from-diagonal application
            auto vertInd = IR(transformBeg,Max(winBeg,ghostCol));
            auto HVertFar = H( vertInd, horzInd );
            Gemm( NORMAL, NORMAL, Real(1), HVertFar, UAccum, WAccum );
            HVertFar = WAccum;

            if( wantSchurVecs )
            {
                auto ZSub = Z( ALL, horzInd );
                Gemm( NORMAL, NORMAL, Real(1), ZSub, UAccum, WAccum );
                ZSub = WAccum;
            }
        }
    }
}

} // namespace hess_qr
