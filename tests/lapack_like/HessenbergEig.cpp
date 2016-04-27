/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

namespace hess_qr {

// Use Ahues and Tissuer's (LAPACK Working Note 122, 1997) refinement of
// Wilkinson's criteria for determining if a subdiagonal entry of a Hessenberg
// matrix is negligible
template<typename F>
Int DetectSmallSubdiagonal( const Matrix<F>& H )
{
    typedef Base<F> Real;

    const Int n = H.Height();
    const F* HBuf = H.LockedBuffer();
    const Int HLDim = H.LDim();

    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real safeMax = Real(1)/safeMin;
    const Real smallNum = safeMin*(Real(n)/ulp);

    // Search up the subdiagonal
    for( Int k=n-1; k>0; --k )
    {
        const F eta_km1_km1 = HBuf[(k-1)+(k-1)*HLDim];
        const F eta_km1_k = HBuf[(k-1)+k*HLDim];
        const Real eta_k_km1 = RealPart(HBuf[k+(k-1)*HLDim]);
        const F eta_k_k = HBuf[k+k*HLDim];
        if( OneAbs(eta_k_km1) <= smallNum )
        {
            return k;
        }

        Real localScale = OneAbs(eta_km1_km1) + OneAbs(eta_k_k);
        if( localScale == Real(0) )
        {
            // Search outward a bit to get a sense of the matrix's local scale
            if( k-2 >= 0 )
                localScale += Abs(RealPart(HBuf[(k-1)+(k-2)*HLDim]));
            if( k+1 <= n-1 )
                localScale += Abs(RealPart(HBuf[(k+1)+ k   *HLDim]));
        }
        
        if( Abs(eta_k_km1) <= ulp*localScale )
        {
            const Real maxOff = Max( Abs(eta_k_km1), OneAbs(eta_km1_k) );
            const Real minOff = Min( Abs(eta_k_km1), OneAbs(eta_km1_k) ); 

            const Real diagDiff = OneAbs(eta_km1_km1-eta_k_k);
            const Real maxDiag = Max( OneAbs(eta_k_k), diagDiff );
            const Real minDiag = Min( OneAbs(eta_k_k), diagDiff );

            const Real sigma = maxDiag + maxOff;
            if( minOff*(maxOff/sigma) <=
                Max(smallNum,ulp*(minDiag*(maxDiag/sigma))) )
            {
                return k;
            }
        }
    }
    return 0;
}

template<typename Real>
Complex<Real> WilkinsonShift( const Matrix<Complex<Real>>& H )
{
    typedef Complex<Real> F;
    const Int n = H.Height();
    DEBUG_ONLY(
      if( n < 2 )
          LogicError("WilkinsonShift requires n >= 2");
    )
    const F* HBuf = H.LockedBuffer(n-2,n-2);
    const Int HLDim = H.LDim();
    const Real zero = Real(0);

    const F eta00 = HBuf[0+0*HLDim];
    const F eta01 = HBuf[0+1*HLDim];
    const F eta10 = HBuf[1+0*HLDim];
    const F eta11 = HBuf[1+1*HLDim];
    // NOTE:
    // eta10 should be real, but it is possibly negative, and so we will
    // interpret it as a complex number so that the square-root is well-defined
    DEBUG_ONLY(
      if( ImagPart(eta10) != zero )
          LogicError("Subdiagonal assumed real");
    )

    F shift = eta11;
    const F gamma = Sqrt(eta01)*Sqrt(eta10);
    Real sigma = OneAbs(gamma);
    if( sigma != zero )
    {
        const F xi = (eta00-shift)/Real(2);   
        const Real xiAbs = OneAbs(xi);
        sigma = Max( sigma, xiAbs );
        const F xiSigma = xi/sigma;
        const F gammaSigma = gamma/sigma;
        F zeta = sigma*Sqrt(xiSigma*xiSigma+gammaSigma*gammaSigma);
        if( xiAbs > zero )
        {
            const F xiUnit = xi/xiAbs;
            if( RealPart(xiUnit)*RealPart(zeta) +
                ImagPart(xiUnit)*ImagPart(zeta) < zero )
            {
                zeta = -zeta;
            }
        }
        shift -= gamma*SafeDiv(gamma,xi+zeta);
    }
    return shift;
}

template<typename Real>
std::tuple<Int,Complex<Real>,Complex<Real>>
ChooseSingleShiftStart( const Matrix<Complex<Real>>& H, Complex<Real> shift )
{
    typedef Complex<Real> F;
    const Real ulp = limits::Precision<Real>();

    const Int n = H.Height();
    const F* HBuf = H.LockedBuffer();
    const Int HLDim = H.LDim();
    for( Int k=n-2; k>0; --k )
    {
        const Real eta10 = RealPart(HBuf[k+(k-1)*HLDim]);
        const F eta11 = HBuf[k+k*HLDim];
        const F eta22 = HBuf[(k+1)+(k+1)*HLDim];
        Real eta21 = RealPart(HBuf[(k+1)+k*HLDim]);

        F eta11Shift = eta11 - shift;
        const Real sigma = OneAbs(eta11Shift) + Abs(eta21);    
        eta11Shift /= sigma;
        eta21 /= sigma; 
        if( Abs(eta10)*Abs(eta21) <=
            ulp*(OneAbs(eta11Shift)*(OneAbs(eta11)+OneAbs(eta22))) ) 
        {
            return std::make_tuple(k,eta11Shift,Complex<Real>(eta21));
        }
    }
    // Start at the top
    const F eta11 = HBuf[0+0*HLDim];
    const F eta22 = HBuf[1+1*HLDim]; 
    Real eta21 = RealPart(HBuf[1+0*HLDim]);

    const F eta11Shift = eta11 - shift;
    const Real sigma = OneAbs(eta11Shift) + Abs(eta21); 
    return std::make_tuple(0,eta11Shift/sigma,Complex<Real>(eta21/sigma)); 
}

template<typename Real>
void PrepareDoubleShift
( Real eta00, Real eta01,
  Real eta10, Real eta11,
  Real& rho0Real, Real& rho0Imag,
  Real& rho1Real, Real& rho1Imag )
{
    const Real zero(0);
    const Real scale = Abs(eta00) + Abs(eta01) + Abs(eta10) + Abs(eta11);
    if( scale == zero )
    {
        rho0Real = rho0Imag = rho1Real = rho1Imag = zero;
    }
    else
    {
        eta00 /= scale;
        eta01 /= scale;
        eta10 /= scale;
        eta11 /= scale;
        Real halfTrace = (eta00+eta11) / 2;
        Real det = (eta00-halfTrace)*(eta11-halfTrace) - eta01*eta10;
        Real absDisc = Sqrt(Abs(det));
        if( det >= zero )
        {
            rho0Real = halfTrace*scale;
            rho0Imag = absDisc*scale;

            rho1Real =  rho0Real;
            rho1Imag = -rho0Imag;
        }
        else
        {
            if( Abs(halfTrace+absDisc-eta11) <= Abs(halfTrace-absDisc-eta11) )
            {
                rho0Real = (halfTrace+absDisc)*scale;
                rho1Real = rho0Real;
            }
            else
            {
                rho1Real = (halfTrace-absDisc)*scale;
                rho0Real = rho1Real;
            }
            rho0Imag = rho1Imag = zero;
        }
    }
}

template<typename Real>
Int ChooseDoubleShiftStart
( const Matrix<Real>& H, 
  const Real& rho0Real, const Real& rho0Imag,
  const Real& rho1Real, const Real& rho1Imag,
        Real* v /* vector of length 3 */ )
{
    const Real ulp = limits::Precision<Real>();
    const Int n = H.Height();
    const Real* HBuf = H.LockedBuffer();
    const Int HLDim = H.LDim();

    for( Int k=n-3; k>=0; --k )
    {
        const Real eta11 = HBuf[ k   + k   *HLDim];
        const Real eta12 = HBuf[ k   +(k+1)*HLDim];
        const Real eta21 = HBuf[(k+1)+ k   *HLDim];
        const Real eta22 = HBuf[(k+1)+(k+1)*HLDim];
        const Real eta32 = HBuf[(k+2)+(k+1)*HLDim];
  
        Real scale = Abs(eta11-rho1Real) + Abs(rho1Imag) + Abs(eta21);
        Real eta21Scale = eta21 / scale;

        v[0] =
          eta21Scale*eta12 +
          (eta11-rho0Real)*((eta11-rho1Real)/scale) -
          rho0Imag*(rho1Imag/scale);
        v[1] = eta21Scale*(eta11+eta22-rho0Real-rho1Real);
        v[2] = eta21Scale*eta32;

        scale = Abs(v[0]) + Abs(v[1]) + Abs(v[2]);
        v[0] /= scale;
        v[1] /= scale;
        v[2] /= scale;
        if( k == 0 )
        {
            return k;
        }
        
        const Real eta00 = HBuf[(k-1)+(k-1)*HLDim];
        const Real eta10 = HBuf[ k   +(k-1)*HLDim];
        if( Abs(eta10)*(Abs(v[1])+Abs(v[2])) <= 
            ulp*Abs(v[0])*(Abs(eta00)+Abs(eta11)+Abs(eta22)) )
        {
            return k;
        }
    }
    // This should never be reached but is to avoid compiler warnings
    return 0;
}

template<typename Real>
void SingleShiftStep
(       Matrix<Complex<Real>>& H,
        Int windowBeg,
        Int windowEnd,
        Int shiftStart,
        Complex<Real> nu0,
        Complex<Real> nu1,
  bool fullTriangle,
        Matrix<Complex<Real>>& Z,
  bool wantSchurVecs,
  Int schurVecBeg,
  Int schurVecEnd )
{
    typedef Complex<Real> F;
    const Int n = H.Height();
    const Int nZ = schurVecEnd-schurVecBeg;
    F* HBuf = H.Buffer();
    F* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    DEBUG_ONLY(
      if( windowBeg > shiftStart || shiftStart >= windowEnd )
          LogicError
          ("shiftStart was ",shiftStart,", but window was [",windowBeg,",",
           windowEnd,")");
    )

    Int transformBeg = ( fullTriangle ? 0 : windowBeg ); 
    Int transformEnd = ( fullTriangle ? n : windowEnd );

    for( Int k=shiftStart; k<windowEnd-1; ++k )
    {
        if( k > shiftStart )
        {
            nu0 = HBuf[k+(k-1)*HLDim];    
            nu1 = RealPart(HBuf[(k+1)+(k-1)*HLDim]);
        }
        // TODO: Assert nu1 is real
        F tau0 = lapack::Reflector( 2, nu0, &nu1, 1 );
        if( k > shiftStart )
        {
            HBuf[ k   +(k-1)*HLDim] = nu0;
            HBuf[(k+1)+(k-1)*HLDim] = 0;
        }
        // The formulas within lapack::Reflector trivially imply that 
        // tau0*Conj(nu1) will be real if nu1 was real on entry to
        // lapack::Reflector (an equivalent claim is made within ZLAHQR)
        F tau1 = RealPart(tau0*Conj(nu1));

        // Apply the Householder reflector from the left
        for( Int j=k; j<transformEnd; ++j )
        {
            F innerProd =
              tau0*HBuf[ k   +j*HLDim] +
              tau1*HBuf[(k+1)+j*HLDim];
            HBuf[ k   +j*HLDim] -= innerProd;
            HBuf[(k+1)+j*HLDim] -= innerProd*nu1;
        }

        // Apply the Householder reflector from the right
        const Int rightApplyEnd = Min(k+3,windowEnd);
        for( Int j=transformBeg; j<rightApplyEnd; ++j )
        {
            F innerProd =
              Conj(tau0)*HBuf[j+ k   *HLDim] +
                   tau1 *HBuf[j+(k+1)*HLDim]; 
            HBuf[j+ k   *HLDim] -= innerProd;
            HBuf[j+(k+1)*HLDim] -= innerProd*Conj(nu1);
        }

        if( wantSchurVecs )
        {
            // Accumulate the Schur vectors
            for( Int j=schurVecBeg; j<schurVecEnd; ++j )
            {
                F innerProd =
                  Conj(tau0)*ZBuf[j+ k   *ZLDim] +
                       tau1 *ZBuf[j+(k+1)*ZLDim];
                ZBuf[j+ k   *ZLDim] -= innerProd;
                ZBuf[j+(k+1)*ZLDim] -= innerProd*Conj(nu1);
            }
        }

        if( k == shiftStart && shiftStart > windowBeg )
        {
            // Make H[shiftStart,shiftStart-1] real by scaling by phase
            // TODO: Investigate more carefully
            F subdiagVal = Real(1) - Conj(tau0);
            F phase = subdiagVal / Abs(subdiagVal);
            HBuf[(shiftStart+1)+shiftStart*HLDim] *= Conj(phase);
            if( shiftStart+2 < windowEnd )
                HBuf[(shiftStart+2)+(shiftStart+1)*HLDim] *= phase;
            for( Int j=shiftStart; j<windowEnd; ++j )
            {
                if( j != shiftStart+1 )
                {
                    if( j+1 < transformEnd )
                    {
                        blas::Scal
                        ( transformEnd-(j+1), phase,
                          &HBuf[j+(j+1)*HLDim], HLDim );
                    }
                    blas::Scal
                    ( j-transformBeg, Conj(phase),
                      &HBuf[transformBeg+j*HLDim], 1 );
                    if( wantSchurVecs )
                    {
                        blas::Scal
                        ( nZ, Conj(phase),
                          &ZBuf[schurVecBeg+j*ZLDim], 1 );
                    }
                }
            }
        }
    }
    // Make H(windowEnd-1,windowEnd-2) real by scaling by phase
    F subdiagVal = HBuf[(windowEnd-1)+(windowEnd-2)*HLDim];
    if( ImagPart(subdiagVal) != Real(0) )
    {
        Real subdiagAbs = Abs(subdiagVal);
        HBuf[(windowEnd-1)+(windowEnd-2)*HLDim] = subdiagAbs;
        F phase = subdiagVal / subdiagAbs;
        if( windowEnd < transformEnd ) 
        {
            blas::Scal
            ( transformEnd-windowEnd, Conj(phase),
              &HBuf[(windowEnd-1)+windowEnd*HLDim], HLDim );
        }
        blas::Scal
        ( (windowEnd-1)-transformBeg, phase,
          &HBuf[transformBeg+(windowEnd-1)*HLDim], 1 );
        if( wantSchurVecs )
        {
            blas::Scal( nZ, phase, &ZBuf[schurVecBeg+(windowEnd-1)*ZLDim], 1 );
        }
    }
}

template<typename Real>
void DoubleShiftStep
(       Matrix<Real>& H,
        Int windowBeg,
        Int windowEnd,
        Int shiftStart,
        Real* v, // vector of length three
  bool fullTriangle,
        Matrix<Real>& Z,
  bool wantSchurVecs,
  Int schurVecBeg,
  Int schurVecEnd )
{
    const Real zero(0), one(1);
    const Int n = H.Height();
    const Int nZ = schurVecEnd-schurVecBeg;
    Real* HBuf = H.Buffer();
    Real* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    DEBUG_ONLY(
      if( windowBeg > shiftStart || shiftStart >= windowEnd-1 )
          LogicError
          ("shiftStart was ",shiftStart,", but window was [",windowBeg,",",
           windowEnd,")");
    )

    Int transformBeg = ( fullTriangle ? 0 : windowBeg ); 
    Int transformEnd = ( fullTriangle ? n : windowEnd );

    for( Int k=shiftStart; k<windowEnd-1; ++k )
    {
        const Int numReflect = Min( 3, windowEnd-k );
        if( k > shiftStart )
        {
            MemCopy( v, &HBuf[k+(k-1)*HLDim], numReflect );
        }
        Real tau0 = lapack::Reflector( numReflect, v[0], &v[1], 1 );
        if( k > shiftStart )
        {
            HBuf[ k   +(k-1)*HLDim] = v[0];
            HBuf[(k+1)+(k-1)*HLDim] = zero;
            if( k < windowEnd-2 )
                HBuf[(k+2)+(k-1)*HLDim] = zero;
        }
        else if( shiftStart > windowBeg )
        {
            // The following is supposedly more reliable than
            // H(k,k-1) = -H(k,k-1) when v(1) and v(2) underflow
            // (cf. LAPACK's {s,d}lahqr)
            HBuf[k+(k-1)*HLDim] *= (one-tau0);
        }
        Real tau1 = tau0*v[1];
        if( numReflect == 3 )
        {
            Real tau2 = tau0*v[2];

            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd =
                       HBuf[ k   +j*HLDim] +
                  v[1]*HBuf[(k+1)+j*HLDim] +
                  v[2]*HBuf[(k+2)+j*HLDim]; 
                HBuf[ k   +j*HLDim] -= innerProd*tau0;
                HBuf[(k+1)+j*HLDim] -= innerProd*tau1;
                HBuf[(k+2)+j*HLDim] -= innerProd*tau2;
            }
           
            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+4,windowEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd =
                       HBuf[j+ k   *HLDim] +
                  v[1]*HBuf[j+(k+1)*HLDim] +
                  v[2]*HBuf[j+(k+2)*HLDim];
                HBuf[j+ k   *HLDim] -= innerProd*tau0;
                HBuf[j+(k+1)*HLDim] -= innerProd*tau1;
                HBuf[j+(k+2)*HLDim] -= innerProd*tau2;
            }

            if( wantSchurVecs )
            {
                for( Int j=schurVecBeg; j<schurVecEnd; ++j )
                {
                    Real innerProd =
                           ZBuf[j+ k   *ZLDim] +
                      v[1]*ZBuf[j+(k+1)*ZLDim] +
                      v[2]*ZBuf[j+(k+2)*ZLDim];
                    ZBuf[j+ k   *ZLDim] -= innerProd*tau0;
                    ZBuf[j+(k+1)*ZLDim] -= innerProd*tau1;
                    ZBuf[j+(k+2)*ZLDim] -= innerProd*tau2;
                }
            }
        }
        else if( numReflect == 2 )
        {
            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd =
                       HBuf[ k   +j*HLDim] +
                  v[1]*HBuf[(k+1)+j*HLDim];
                HBuf[ k   +j*HLDim] -= innerProd*tau0;
                HBuf[(k+1)+j*HLDim] -= innerProd*tau1;
            }

            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+3,windowEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd =
                       HBuf[j+ k   *HLDim] +
                  v[1]*HBuf[j+(k+1)*HLDim]; 
                HBuf[j+ k   *HLDim] -= innerProd*tau0;
                HBuf[j+(k+1)*HLDim] -= innerProd*tau1;
            }

            if( wantSchurVecs )
            {
                // Accumulate the Schur vectors
                for( Int j=schurVecBeg; j<schurVecEnd; ++j )
                {
                    Real innerProd =
                           ZBuf[j+ k   *ZLDim] +
                      v[1]*ZBuf[j+(k+1)*ZLDim];
                    ZBuf[j+ k   *ZLDim] -= innerProd*tau0;
                    ZBuf[j+(k+1)*ZLDim] -= innerProd*tau1;
                }
            }
        }
    }
}


template<typename Real>
void Windowed
(       Matrix<Real>& H,
        Int windowBeg,
        Int windowEnd,
        Matrix<Real>& wReal,
        Matrix<Real>& wImag,
  bool fullTriangle,
        Matrix<Real>& Z,
  bool wantSchurVecs,
  Int schurVecBeg,
  Int schurVecEnd )
{
    const Real zero(0), threeFourths=Real(3)/Real(4);
    const Int maxIter=30;    
    // Cf. LAPACK for these somewhat arbitrary constants
    const Real exceptScale0=Real(3)/Real(4),
               exceptScale1=Real(-4375)/Real(10000);
    const Int n = H.Height();
    const Int windowSize = windowEnd - windowBeg;
    const Int nZ = schurVecEnd - schurVecBeg;
    Real* HBuf = H.Buffer(); 
    Real* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    wReal.Resize( n, 1 );
    wImag.Resize( n, 1 );
    Real* wRealBuf = wReal.Buffer();
    Real* wImagBuf = wImag.Buffer();

    if( windowSize == 0 )
    {
        return;
    }
    if( windowSize == 1 )
    {
        wRealBuf[windowBeg] = HBuf[windowBeg+windowBeg*HLDim];
        wImagBuf[windowBeg] = zero;
        return;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=windowBeg; j<windowEnd-3; ++j ) 
    {
        HBuf[(j+2)+j*HLDim] = 0;
        HBuf[(j+3)+j*HLDim] = 0;
    }
    if( windowBeg <= windowEnd-3 )
        HBuf[(windowEnd-1)+(windowEnd-3)*HLDim] = 0;
    
    // Attempt to converge the eigenvalues one or two at a time
    Real v[3]; // for computing Householder reflectors
    while( windowBeg < windowEnd )
    {
        Int iterBeg = windowBeg;
        Int iter;
        for( iter=0; iter<maxIter; ++iter )
        {
            {
                auto winInd = IR(iterBeg,windowEnd);
                iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
            }
            if( iterBeg > windowBeg )
            {
                HBuf[iterBeg+(iterBeg-1)*HLDim] = zero;
            }
            if( iterBeg == windowEnd-1 )
            {
                wRealBuf[iterBeg] = HBuf[iterBeg+iterBeg*HLDim];
                wImagBuf[iterBeg] = zero;
                --windowEnd;
                break;
            }
            else if( iterBeg == windowEnd-2 )
            {
                Real c, s;
                Real& alpha00 = HBuf[(windowEnd-2)+(windowEnd-2)*HLDim];
                Real& alpha01 = HBuf[(windowEnd-2)+(windowEnd-1)*HLDim];
                Real& alpha10 = HBuf[(windowEnd-1)+(windowEnd-2)*HLDim];
                Real& alpha11 = HBuf[(windowEnd-1)+(windowEnd-1)*HLDim];
                lapack::TwoByTwoSchur
                ( alpha00, alpha01,
                  alpha10, alpha11,
                  c, s,
                  wRealBuf[iterBeg],   wImagBuf[iterBeg], 
                  wRealBuf[iterBeg+1], wImagBuf[iterBeg+1] );
                if( fullTriangle )
                {
                    if( n > windowEnd )
                        blas::Rot
                        ( n-windowEnd,
                          &HBuf[(windowEnd-2)+windowEnd*HLDim], HLDim,
                          &HBuf[(windowEnd-1)+windowEnd*HLDim], HLDim,
                          c, s );
                    blas::Rot
                    ( windowEnd-2,
                      &HBuf[0+(windowEnd-2)*HLDim], 1,
                      &HBuf[0+(windowEnd-1)*HLDim], 1,
                      c, s );
                }
                if( wantSchurVecs )
                {
                    blas::Rot
                    ( nZ,
                      &ZBuf[schurVecBeg+(windowEnd-2)*ZLDim], 1,
                      &ZBuf[schurVecBeg+(windowEnd-1)*ZLDim], 1,
                      c, s );
                }
                windowEnd -= 2;
                break;
            }

            // Pick either the Francis shifts or exceptional shifts
            Real eta00, eta01, eta10, eta11;
            if( iter == maxIter/3 )
            {
                const Real* HSub = &HBuf[iterBeg+iterBeg*HLDim];
                const Real scale = Abs(HSub[1+0*HLDim]) + Abs(HSub[2+1*HLDim]);
                eta00 = exceptScale0*scale + HSub[0+0*HLDim];
                eta01 = exceptScale1*scale;
                eta10 = scale;
                eta11 = eta00; 
            } 
            else if( iter == 2*maxIter/3 )
            {
                const Real* HSub = &HBuf[(windowEnd-3)+(windowEnd-3)*HLDim];
                const Real scale = Abs(HSub[2+1*HLDim]) + Abs(HSub[1+0*HLDim]);
                eta00 = exceptScale0*scale + HSub[2+2*HLDim];
                eta11 = exceptScale1*scale;
                eta10 = scale;
                eta11 = eta00;
            }
            else
            {
                const Real* HSub = &HBuf[(windowEnd-2)+(windowEnd-2)*HLDim];
                eta00 = HSub[0+0*HLDim];
                eta01 = HSub[0+1*HLDim];
                eta10 = HSub[1+0*HLDim];
                eta11 = HSub[1+1*HLDim];
            }
            Real rho0Real, rho0Imag, rho1Real, rho1Imag;
            PrepareDoubleShift
            ( eta00, eta01,
              eta10, eta11,
              rho0Real, rho0Imag,
              rho1Real, rho1Imag );

            auto subInd = IR(iterBeg,windowEnd);
            Int shiftStart = iterBeg +
              ChooseDoubleShiftStart
              ( H(subInd,subInd), rho0Real, rho0Imag, rho1Real, rho1Imag, v );
            DoubleShiftStep
            ( H, iterBeg, windowEnd, shiftStart, v, fullTriangle,
              Z, wantSchurVecs, schurVecBeg, schurVecEnd );
        }
        if( iter == maxIter )
            RuntimeError("QR iteration did not converge");
    }
}

template<typename Real>
void Windowed
(       Matrix<Complex<Real>>& H,
        Int windowBeg,
        Int windowEnd,
        Matrix<Complex<Real>>& w,
  bool fullTriangle,
        Matrix<Complex<Real>>& Z,
  bool wantSchurVecs,
  Int schurVecBeg,
  Int schurVecEnd )
{
    typedef Complex<Real> F;
    const Real zero(0), threeFourths=Real(3)/Real(4);
    const Int maxIter=30;    
    const Int n = H.Height();
    const Int windowSize = windowEnd - windowBeg;
    const Int nZ = schurVecEnd - schurVecBeg;
    F* HBuf = H.Buffer(); 
    F* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    w.Resize( n, 1 );
    F* wBuf = w.Buffer();

    if( windowSize == 0 )
    {
        return;
    }
    if( windowSize == 1 )
    {
        wBuf[windowBeg] = HBuf[windowBeg+windowBeg*HLDim];
        return;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=windowBeg; j<windowEnd-3; ++j ) 
    {
        HBuf[(j+2)+j*HLDim] = 0;
        HBuf[(j+3)+j*HLDim] = 0;
    }
    if( windowBeg <= windowEnd-3 )
        HBuf[(windowEnd-1)+(windowEnd-3)*HLDim] = 0;
    
    // Rotate the matrix so that the subdiagonals are real
    const Int scaleBeg = ( fullTriangle ? 0 : windowBeg );
    const Int scaleEnd = ( fullTriangle ? n : windowEnd );
    for( Int i=windowBeg+1; i<windowEnd; ++i )
    {
        if( ImagPart(HBuf[i+(i-1)*HLDim]) != zero )
        {
            const F eta_i_im1 = HBuf[i+(i-1)*HLDim];
            F phase = eta_i_im1 / OneAbs(eta_i_im1);
            phase = Conj(phase) / Abs(phase);
            HBuf[i+(i-1)*HLDim] = Abs(eta_i_im1);
            blas::Scal( scaleEnd-i, phase, &HBuf[i+i*HLDim], HLDim );
            blas::Scal
            ( Min(scaleEnd,i+2)-scaleBeg, Conj(phase),
              &HBuf[scaleBeg+i*HLDim], 1 );
            if( wantSchurVecs )
                blas::Scal( nZ, Conj(phase), &ZBuf[schurVecBeg+i*ZLDim], 1 );
        }
    }

    // Attempt to converge the eigenvalues one at a time
    while( windowBeg < windowEnd )
    {
        Int iterBeg = windowBeg;
        Int iter;
        for( iter=0; iter<maxIter; ++iter )
        {
            {
                auto winInd = IR(iterBeg,windowEnd);
                iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
            }
            if( iterBeg > windowBeg )
            {
                HBuf[iterBeg+(iterBeg-1)*HLDim] = zero;
            }
            if( iterBeg == windowEnd-1 )
            {
                wBuf[iterBeg] = HBuf[iterBeg+iterBeg*HLDim];
                --windowEnd;
                break;
            }

            // Pick either the Wilkinson shift or an exceptional shift
            F shift;
            auto subInd = IR(iterBeg,windowEnd);
            if( iter == maxIter/3 )
            {
                const F diagVal = HBuf[iterBeg+iterBeg*HLDim];
                const Real subdiagVal =
                  RealPart(HBuf[(iterBeg+1)+iterBeg*HLDim]);
                shift = diagVal + threeFourths*Abs(subdiagVal);
            } 
            else if( iter == 2*maxIter/3 )
            {
                const F diagVal = HBuf[(windowEnd-1)+(windowEnd-1)*HLDim];
                const Real subdiagVal =
                  RealPart(HBuf[(windowEnd-1)+(windowEnd-2)*HLDim]);
                shift = diagVal + threeFourths*Abs(subdiagVal);
            }
            else
            {
                shift = WilkinsonShift( H(subInd,subInd) );
            }

            auto qrTuple = ChooseSingleShiftStart( H(subInd,subInd), shift );
            const Int shiftStart = iterBeg+std::get<0>(qrTuple);
            const Complex<Real> nu0 = std::get<1>(qrTuple);
            const Complex<Real> nu1 = std::get<2>(qrTuple);
            SingleShiftStep
            ( H, iterBeg, windowEnd, shiftStart, nu0, nu1, fullTriangle,
              Z, wantSchurVecs, schurVecBeg, schurVecEnd );
        }
        if( iter == maxIter )
            RuntimeError("QR iteration did not converge");
    }
}

} // namespace hess_qr

template<typename Real>
void HessenbergQR
(       Matrix<Real>& H,
        Matrix<Real>& wReal,
        Matrix<Real>& wImag,
  bool fullTriangle=true )
{
    Int windowBeg=0, windowEnd=H.Height();
    Int schurVecBeg=0, schurVecEnd=0;
    bool wantSchurVecs=false;
    Matrix<Real> Z;
    hess_qr::Windowed
    ( H, windowBeg, windowEnd, wReal, wImag, fullTriangle, 
      Z, wantSchurVecs, schurVecBeg, schurVecEnd );
}

template<typename Real>
void HessenbergQR
(       Matrix<Real>& H,
        Matrix<Real>& wReal,
        Matrix<Real>& wImag,
        Matrix<Real>& Z,
  bool fullTriangle=true )
{
    Int windowBeg=0, windowEnd=H.Height();
    Int schurVecBeg=0, schurVecEnd=Z.Height();
    bool wantSchurVecs=true;
    hess_qr::Windowed
    ( H, windowBeg, windowEnd, wReal, wImag, fullTriangle, 
      Z, wantSchurVecs, schurVecBeg, schurVecEnd );
}

template<typename Real>
void HessenbergQR
(       Matrix<Complex<Real>>& H,
        Matrix<Complex<Real>>& w,
  bool fullTriangle=true )
{
    Int windowBeg=0, windowEnd=H.Height();
    Int schurVecBeg=0, schurVecEnd=0;
    bool wantSchurVecs=false;
    Matrix<Complex<Real>> Z;
    hess_qr::Windowed
    ( H, windowBeg, windowEnd, w, fullTriangle, 
      Z, wantSchurVecs, schurVecBeg, schurVecEnd );
}

template<typename Real>
void HessenbergQR
(       Matrix<Complex<Real>>& H,
        Matrix<Complex<Real>>& w,
        Matrix<Complex<Real>>& Z,
  bool fullTriangle=true )
{
    Int windowBeg=0, windowEnd=H.Height();
    Int schurVecBeg=0, schurVecEnd=Z.Height();
    bool wantSchurVecs=true;
    hess_qr::Windowed
    ( H, windowBeg, windowEnd, w, fullTriangle, 
      Z, wantSchurVecs, schurVecBeg, schurVecEnd );
}

template<typename Real>
void TestAhuesTisseur( bool print )
{
    typedef Complex<Real> F;
    const Int n = 3;
    Output("Testing Ahues/Tisseur with ",TypeName<F>());

    Matrix<F> H;
    Zeros( H, n, n );
    H.Set( 0, 0, F(1.) );
    H.Set( 0, 1, F(1.1e5) );
    H.Set( 0, 2, F(0.) );
    H.Set( 1, 0, F(1.1e-8) );
    H.Set( 1, 1, F(1.+1.e-2) );
    H.Set( 1, 2, F(1.1e5) );
    H.Set( 2, 0, F(0.) );
    H.Set( 2, 1, F(1.1e-8) );
    H.Set( 2, 2, F(1.+2.*1.e-2) );
    const Real HFrob = FrobeniusNorm( H );
    Output("|| H ||_F = ",HFrob);
    if( print )
        Print( H, "H" );

    Matrix<F> T, w, Z;
    Identity( Z, n, n );
    T = H;
    Timer timer;
    timer.Start();
    HessenbergQR( T, w, Z );
    Output("  HessenbergQR: ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "w" );
        Print( Z, "Z" );
        Print( T, "T" );
    }

    Matrix<F> R;
    Gemm( NORMAL, NORMAL, F(1), Z, T, R );
    Gemm( NORMAL, NORMAL, F(1), H, Z, F(-1), R );
    const Real errFrob = FrobeniusNorm( R ); 
    Output("|| H Z - Z T ||_F / || H ||_F = ",errFrob/HFrob);
    if( print )
        Print( R );
}

template<typename Real>
void TestAhuesTisseurQuasi( bool print )
{
    typedef Real F;
    const Int n = 3;
    Output("Testing Ahues/Tisseur with ",TypeName<F>());

    Matrix<F> H;
    Zeros( H, n, n );
    H.Set( 0, 0, F(1.) );
    H.Set( 0, 1, F(1.1e5) );
    H.Set( 0, 2, F(0.) );
    H.Set( 1, 0, F(1.1e-8) );
    H.Set( 1, 1, F(1.+1.e-2) );
    H.Set( 1, 2, F(1.1e5) );
    H.Set( 2, 0, F(0.) );
    H.Set( 2, 1, F(1.1e-8) );
    H.Set( 2, 2, F(1.+2.*1.e-2) );
    const Real HFrob = FrobeniusNorm( H );
    Output("|| H ||_F = ",HFrob);
    if( print )
        Print( H, "H" );

    Matrix<F> T, wReal, wImag, Z;
    Identity( Z, n, n );
    T = H;
    Timer timer;
    timer.Start();
    HessenbergQR( T, wReal, wImag, Z );
    Output("  HessenbergQR: ",timer.Stop()," seconds");
    if( print )
    {
        Print( wReal, "wReal" );
        Print( wImag, "wImag" );
        Print( Z, "Z" );
        Print( T, "T" );
    }

    Matrix<F> R;
    Gemm( NORMAL, NORMAL, F(1), Z, T, R );
    Gemm( NORMAL, NORMAL, F(1), H, Z, F(-1), R );
    const Real errFrob = FrobeniusNorm( R ); 
    Output("|| H Z - Z T ||_F / || H ||_F = ",errFrob/HFrob);
    if( print )
        Print( R );
}

template<typename Real>
void TestRandom( Int n, bool print )
{
    typedef Complex<Real> F;
    Output("Testing uniform Hessenberg with ",TypeName<F>());

    Matrix<F> H;
    Uniform( H, n, n );
    MakeTrapezoidal( UPPER, H, -1 );
    const Real HFrob = FrobeniusNorm( H );
    Output("|| H ||_F = ",HFrob);
    if( print )
        Print( H, "H" );

    Matrix<F> T, w, Z;
    Identity( Z, H.Height(), H.Height() );
    T = H;
    Timer timer;
    timer.Start();
    HessenbergQR( T, w, Z );
    Output("  HessenbergQR: ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "w" );
        Print( Z, "Z" );
        Print( T, "T" );
    }

    Matrix<F> R;
    Gemm( NORMAL, NORMAL, F(1), Z, T, R );
    Gemm( NORMAL, NORMAL, F(1), H, Z, F(-1), R );
    const Real errFrob = FrobeniusNorm( R ); 
    Output("|| H Z - Z T ||_F / || H ||_F = ",errFrob/HFrob);
    if( print )
        Print( R );
}

template<typename Real>
void TestRandomQuasi( Int n, bool print )
{
    typedef Real F;
    Output("Testing uniform Hessenberg with ",TypeName<F>());

    Matrix<F> H;
    Uniform( H, n, n );
    MakeTrapezoidal( UPPER, H, -1 );
    const Real HFrob = FrobeniusNorm( H );
    Output("|| H ||_F = ",HFrob);
    if( print )
        Print( H, "H" );

    Matrix<F> T, wReal, wImag, Z;
    Identity( Z, H.Height(), H.Height() );
    T = H;
    Timer timer;
    timer.Start();
    HessenbergQR( T, wReal, wImag, Z );
    Output("  HessenbergQR: ",timer.Stop()," seconds");
    if( print )
    {
        Print( wReal, "wReal" );
        Print( wImag, "wImag" );
        Print( Z, "Z" );
        Print( T, "T" );
    }

    Matrix<F> R;
    Gemm( NORMAL, NORMAL, F(1), Z, T, R );
    Gemm( NORMAL, NORMAL, F(1), H, Z, F(-1), R );
    const Real errFrob = FrobeniusNorm( R ); 
    Output("|| H Z - Z T ||_F / || H ||_F = ",errFrob/HFrob);
    if( print )
        Print( R );
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","random matrix size",80);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        TestAhuesTisseurQuasi<float>( print );
        TestAhuesTisseurQuasi<double>( print );
#ifdef EL_HAVE_QUAD
        TestAhuesTisseurQuasi<Quad>( print );
#endif
#ifdef EL_HAVE_QD
        TestAhuesTisseurQuasi<DoubleDouble>( print );
        TestAhuesTisseurQuasi<QuadDouble>( print );
#endif
#ifdef EL_HAVE_MPC
        TestAhuesTisseurQuasi<BigFloat>( print );
#endif

        TestAhuesTisseur<float>( print );
        TestAhuesTisseur<double>( print );
#ifdef EL_HAVE_QUAD
        TestAhuesTisseur<Quad>( print );
#endif

        TestRandomQuasi<float>( n, print );
        TestRandomQuasi<double>( n, print );
#ifdef EL_HAVE_QUAD
        TestRandomQuasi<Quad>( n, print );
#endif
#ifdef EL_HAVE_QD
        TestRandomQuasi<DoubleDouble>( n, print );
        TestRandomQuasi<QuadDouble>( n, print );
#endif
#ifdef EL_HAVE_QD
        TestRandomQuasi<BigFloat>( n, print );
#endif

        TestRandom<float>( n, print );
        TestRandom<double>( n, print );
#ifdef EL_HAVE_QUAD
        TestRandom<Quad>( n, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
