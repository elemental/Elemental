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
    DEBUG_ONLY(CSE cse("small_bulge_sweep::ImplicitQQuadraticSeed"))
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
    DEBUG_ONLY(CSE cse("small_bulge_sweep::IntroduceBulge"))
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
void PairShifts( vector<Real>& realShifts, vector<Real>& imagShifts )
{
    DEBUG_ONLY(CSE cse("small_bulge_sweep::PairShifts"))
    const Int numShifts = realShifts.size();
    const Real zero(0);

    // TODO: Make this routine simpler

    // Count the number of real shifts
    Int numRealShifts = 0;
    for( Int i=0; i<numShifts; ++i )
    {
        if( imagShifts[i] == zero )
            ++numRealShifts;
    }
    DEBUG_ONLY(
      if( numRealShifts % 2 != 0 )
          LogicError("Expected an even number of real shifts");
    )
    // Group the real shifts together up front, followed by conjugate pairs
    vector<Real> realShiftsSort(numShifts), imagShiftsSort(numShifts);
    Int realOff=0, complexOff=numRealShifts;
    for( Int i=0; i<numShifts; ++i )
    {
        if( imagShifts[i] == zero ) 
        {
            realShiftsSort[realOff] = realShifts[i];
            imagShiftsSort[realOff] = zero;
            ++realOff;
        }
        else
        {
            realShiftsSort[complexOff] = realShifts[i]; 
            imagShiftsSort[complexOff] = imagShifts[i];
            ++complexOff;
        }
    }
    DEBUG_ONLY(
      for( Int i=numRealShifts; i<numShifts; i+=2 )
      {
          if( imagShiftsSort[i] != -imagShiftsSort[i+1] )
              LogicError("Complex shifts did not come in conjugate pairs");
      }
    )

    realShifts = realShiftsSort;
    imagShifts = imagShiftsSort;
}

template<typename Real>
void ComputeReflectors
( Matrix<Real>& H,
  Int winBeg,
  vector<Real>& realShifts,
  vector<Real>& imagShifts,
  Matrix<Real>& W,
  Int packetBeg,
  Int firstFull,
  Int lastFull,
  bool have2x2 )
{
    DEBUG_ONLY(CSE cse("small_bulge_sweep::ComputeReflectors"))
    const Real zero(0);
    const Real ulp = limits::Precision<Real>();

    // Set aside space for a Householder vector for a candidate reinflation of
    // a deflated bulge
    vector<Real> wCand(3);

    if( have2x2 )
    {
        const Int bulge = lastFull+1;
        const Real realShift0 = realShifts[2*bulge];
        const Real imagShift0 = imagShifts[2*bulge];
        const Real realShift1 = realShifts[2*bulge+1];
        const Real imagShift1 = imagShifts[2*bulge+1];
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
            // Find the reflection for chasing the 2x2 bulge
            Real beta = H( bulgeBeg+1, bulgeBeg );
            w[1] = H( bulgeBeg+2, bulgeBeg );
            w[0] = lapack::Reflector( 2, beta, &w[1], 1 );
            H( bulgeBeg+1, bulgeBeg ) = beta;
            H( bulgeBeg+2, bulgeBeg ) = zero;
        }
    }
    for( Int bulge=lastFull; bulge>=firstFull; --bulge )
    {
        const Real realShift0 = realShifts[2*bulge];
        const Real imagShift0 = imagShifts[2*bulge];
        const Real realShift1 = realShifts[2*bulge+1];
        const Real imagShift1 = imagShifts[2*bulge+1];
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
            Real& eta00 = H(bulgeBeg,  bulgeBeg  );
            Real& eta10 = H(bulgeBeg+1,bulgeBeg  );
            Real& eta11 = H(bulgeBeg+1,bulgeBeg+1);
            Real& eta20 = H(bulgeBeg+2,bulgeBeg  );
            Real& eta22 = H(bulgeBeg+2,bulgeBeg+2);
            Real& eta30 = H(bulgeBeg+3,bulgeBeg  );
            Real& eta31 = H(bulgeBeg+3,bulgeBeg+1);
            Real& eta32 = H(bulgeBeg+3,bulgeBeg+2);

            Real beta = eta10;
            w[1] = eta20;
            w[2] = eta30;
            w[0] = lapack::Reflector( 3, beta, &w[1], 1 );
                    
            // "Vigilantly" search for a deflation
            // (deflation within the interior is exceedingly rare,
            //  but this check should be essentially free)
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
  Int firstFull,
  Int lastFull,
  bool have2x2,
  bool accumulate )
{
    DEBUG_ONLY(CSE cse("small_bulge_sweep::ApplyReflectors"))

    // Apply from the left
    // ===================
    if( have2x2 )
    {
        const Int bulge = lastFull+1;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const Real* w = &W(0,bulge);

        DEBUG_ONLY(
          if( bulgeBeg+1 < winBeg )
              LogicError("bulgeBeg=",bulgeBeg,", winBeg=",winBeg);
        )
        for( Int j=bulgeBeg+1; j<transformEnd; ++j )
        {
            Real& eta1 = H(bulgeBeg+1,j);
            Real& eta2 = H(bulgeBeg+2,j);

            // Apply I - tau*[1;v]*[1;v]', where tau is stored in
            // w[0] and v is stored in w[1]
            Real innerProd = w[0]*(eta1+w[1]*eta2);
            eta1 -=      innerProd;
            eta2 -= w[1]*innerProd;
        }
    }
    for( Int j=Max(packetBeg,winBeg); j<transformEnd; ++j )
    {
        const Int lastBulge = Min( lastFull, (j-packetBeg-1)/3 );
        for( Int bulge=firstFull; bulge<=lastBulge; ++bulge )
        {
            const Int bulgeBeg = packetBeg + 3*bulge;
            const Real* w = &W(0,bulge);
            Real& eta1 = H(bulgeBeg+1,j); 
            Real& eta2 = H(bulgeBeg+2,j);
            Real& eta3 = H(bulgeBeg+3,j);

            // Apply I - tau*[1;v]*[1;v]', where tau is stored in
            // w[0] and v is stored in w[1] and w[2]
            Real innerProd = w[0]*(eta1+w[1]*eta2+w[2]*eta3);
            eta1 -=      innerProd;
            eta2 -= w[1]*innerProd;
            eta3 -= w[2]*innerProd;
        }
    }

    // Apply from the right
    // ====================
    const Int nZ = Z.Height();
    // The first relative index of the slab that is in the window
    // (with 4x4 bulges being introduced at winBeg-1)
    const Int slabRelBeg = Max(0,(winBeg-1)-ghostCol);
    if( have2x2 )
    {
        const Int bulge = lastFull+1;
        const Int bulgeBeg = packetBeg + 3*bulge;
        const Real* w = &W(0,bulge);
        auto HR = H(ALL,IR(bulgeBeg,END));

        for( Int j=transformBeg; j<Min(winEnd,bulgeBeg+3); ++j )
        {
            Real& eta1 = HR(j,1);
            Real& eta2 = HR(j,2);

            Real innerProd = w[0]*(eta1+eta2*w[1]);
            eta1 -= innerProd;
            eta2 -= innerProd*w[1];
        }

        if( accumulate )
        {
            Int kU = bulgeBeg - ghostCol;
            auto UR = U(ALL,IR(kU,END));
            for( Int j=slabRelBeg; j<slabSize; ++j ) 
            {
                Real& ups1 = UR(j,0);
                Real& ups2 = UR(j,1);

                Real innerProd = w[0]*(ups1+ups2*w[1]);
                ups1 -= innerProd;
                ups2 -= innerProd*w[1];
            }
        }
        else if( wantSchurVecs )
        {
            auto ZR = Z(ALL,IR(bulgeBeg,END));
            for( Int j=0; j<nZ; ++j )
            {
                Real& zeta1 = ZR(j,1);
                Real& zeta2 = ZR(j,2);

                Real innerProd = w[0]*(zeta1+zeta2*w[1]);
                zeta1 -= innerProd;
                zeta2 -= innerProd*w[1];
            }
        }
    }
    for( Int bulge=lastFull; bulge>=firstFull; --bulge )
    {
        const Int bulgeBeg = packetBeg + 3*bulge;
        const Real* w = &W(0,bulge);
        auto HR = H(ALL,IR(bulgeBeg,END));

        for( Int j=transformBeg; j<Min(winEnd,bulgeBeg+4); ++j )
        {
            Real& eta1 = HR(j,1);
            Real& eta2 = HR(j,2);
            Real& eta3 = HR(j,3);

            Real innerProd = w[0]*(eta1+eta2*w[1]+eta3*w[2]);
            eta1 -= innerProd;
            eta2 -= innerProd*w[1];
            eta3 -= innerProd*w[2];
        }

        if( accumulate )
        {
            Int kU = bulgeBeg - ghostCol;
            auto UR = U(ALL,IR(kU,END));
            for( Int j=slabRelBeg; j<slabSize; ++j ) 
            {
                Real& ups1 = UR(j,0);
                Real& ups2 = UR(j,1);
                Real& ups3 = UR(j,2);

                Real innerProd = w[0]*(ups1+ups2*w[1]+ups3*w[2]);
                ups1 -= innerProd;
                ups2 -= innerProd*w[1];
                ups3 -= innerProd*w[2];
            }
        }
        else if( wantSchurVecs )
        {
            auto ZR = Z(ALL,IR(bulgeBeg,END));
            for( Int j=0; j<nZ; ++j )
            {
                Real& zeta1 = ZR(j,1);
                Real& zeta2 = ZR(j,2);
                Real& zeta3 = ZR(j,3);

                Real innerProd = w[0]*(zeta1+zeta2*w[1]+zeta3*w[2]);
                zeta1 -= innerProd;
                zeta2 -= innerProd*w[1];
                zeta3 -= innerProd*w[2];
            }
        }
    }
}

template<typename Real>
void VigilantDeflation
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  Int packetBeg,
  Int vigFirst,
  Int vigLast )
{
    DEBUG_ONLY(CSE cse("small_bulge_sweep::VigilantDeflation"))
    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(H.Height())/ulp);

    // NOTE: LAPACK has a strange increment of vigLast by one if
    //       packetBeg is equal to winBeg-3 here...
    for( Int bulge=vigLast; bulge>=vigFirst; --bulge )
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

} // namespace small_bulge_sweep

template<typename Real>
void SmallBulgeSweep
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  vector<Real>& realShifts,
  vector<Real>& imagShifts,
  bool fullTriangle,
  Matrix<Real>& Z,
  bool wantSchurVecs,
  Matrix<Real>& U,
  Matrix<Real>& W,
  Matrix<Real>& WAccum,
  bool accumulate=true )
{
    const Real zero(0);
    const Int n = H.Height();
    const Int nZ = Z.Height();

    const Int numShifts = realShifts.size();
    DEBUG_ONLY(
      if( numShifts < 2 )
          LogicError("Expected at least one pair of shifts..."); 
      if( numShifts % 2 != 0 )
          LogicError("Expected an even number of sweeps");
    )
    const Int numBulges = numShifts / 2;

    PairShifts( realShifts, imagShifts );

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
            const Int firstFull = Max( 0, ((winBeg-1)-packetBeg+2)/3 );
            const Int numFull = Min( numBulges, (winEnd-packetBeg-1)/3 );
            const Int lastFull = firstFull+numFull-1;

            // 2x2 reflectors can only occur if a bulge occupies the 3x3 in the
            // bottom-right corner (which begins at index winEnd-3)
            const bool have2x2 =
              ( (lastFull+1) < numBulges &&
                packetBeg+3*(lastFull+1) == winEnd-3 );

            small_bulge_sweep::ComputeReflectors
            ( H, winBeg, realShifts, imagShifts, W,
              packetBeg, firstFull, lastFull, have2x2 );

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
              firstFull, lastFull, have2x2, accumulate );

            // Vigilant deflation check using Ahues/Tisseur criteria
            const Int vigFirst =
              ( packetBeg+3*firstFull == winBeg-1 ? firstFull+1 : firstFull );
            const Int vigLast = ( have2x2 ? lastFull+1 : lastFull );
            small_bulge_sweep::VigilantDeflation
            ( H, winBeg, winEnd, packetBeg, vigFirst, vigLast );

            // Fill in the last row of the translated full bulges
            //
            // | X X X X X |     | X X X X X |
            // | X X X X X |     | X X X X X |
            // | X X X X X | |-> |   X X X X |
            // | X X X X X |     |   X X X X |
            // |       X X |     |   X X X X |
            //
            // where the last row is only introduced from the transformation
            //
            //   H(k+4,k+1:k+3) -= tau H(k+4,k+1:k+3) [1; v0; v1] [1; v0; v1]'.
            //
            // For convenience, since H(k+4,k+1:k+2)=0, we start by computing 
            //
            //   innerProd = tau H(k+4,k+1:k+3) [1; v0; v1]
            //             = tau H(k+4,k+3) v1
            //
            for( Int bulge=lastFull; bulge>=firstFull; --bulge )
            {
                const Int k = packetBeg + 3*bulge;
                const Real* w = &W(0,bulge);

                const Real innerProd = w[0]*w[2]*H(k+4,k+3);
                H(k+4,k+1) = -innerProd;
                H(k+4,k+2) = -innerProd*w[1];
                H(k+4,k+3) -= innerProd*w[2];
            }
        }
        if( accumulate )
        {
            // TODO: Debug this section's index ranges
            const Int transformBeg = ( fullTriangle ? 0 : winBeg );
            const Int transformEnd = ( fullTriangle ? n : winEnd );
 
            const Int slabRelBeg = Max(0,(winBeg-1)-ghostCol);
            const Int slabWinEnd = Min(slabEnd,winEnd);

            const Int nU = (slabSize-1) - Max(0,slabEnd-winEnd) - slabRelBeg;

            auto contractInd = IR(0,nU) + slabRelBeg;
            auto UAccum = U( contractInd, contractInd );

            // Horizontal far-from-diagonal application
            const Int horzSize = transformEnd-transformBeg;
            auto horzInd = IR(0,nU) + (ghostCol+slabRelBeg);
            auto HHorzFar = H( horzInd, IR(0,horzSize)+slabWinEnd );
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
