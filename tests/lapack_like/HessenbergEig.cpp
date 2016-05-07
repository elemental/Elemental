/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
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
    const Int offset = n-2;
    const F* HBuf = H.LockedBuffer(offset,offset);
    const Int HLDim = H.LDim();
    const Real zero = Real(0);

    const F& eta00 = HBuf[0+0*HLDim];
    const F& eta01 = HBuf[0+1*HLDim];
    const F& eta10 = HBuf[1+0*HLDim];
    const F& eta11 = HBuf[1+1*HLDim];
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
ChooseSingleShiftStart
( const Matrix<Complex<Real>>& H,
        Complex<Real> shift )
{
    typedef Complex<Real> F;
    const Real ulp = limits::Precision<Real>();

    const Int n = H.Height();
    const F* HBuf = H.LockedBuffer();
    const Int HLDim = H.LDim();
    for( Int k=n-2; k>0; --k )
    {
        const Real eta10 = RealPart(HBuf[k+(k-1)*HLDim]);
        const F& eta11 = HBuf[k+k*HLDim];
        const F& eta22 = HBuf[(k+1)+(k+1)*HLDim];
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
    const F& eta11 = HBuf[0+0*HLDim];
    const F& eta22 = HBuf[1+1*HLDim]; 
    Real eta21 = RealPart(HBuf[1+0*HLDim]);

    const F eta11Shift = eta11 - shift;
    const Real sigma = OneAbs(eta11Shift) + Abs(eta21); 
    return std::make_tuple(0,eta11Shift/sigma,Complex<Real>(eta21/sigma)); 
}

template<typename Real>
void PrepareDoubleShift
( Real eta00, Real eta01,
  Real eta10, Real eta11,
  Real& shift0Real, Real& shift0Imag,
  Real& shift1Real, Real& shift1Imag )
{
    const Real zero(0);
    const Real scale = Abs(eta00) + Abs(eta01) + Abs(eta10) + Abs(eta11);
    if( scale == zero )
    {
        shift0Real = shift0Imag = shift1Real = shift1Imag = zero;
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
            shift0Real = halfTrace*scale;
            shift0Imag = absDisc*scale;

            shift1Real =  shift0Real;
            shift1Imag = -shift0Imag;
        }
        else
        {
            if( Abs(halfTrace+absDisc-eta11) <= Abs(halfTrace-absDisc-eta11) )
            {
                shift0Real = (halfTrace+absDisc)*scale;
                shift1Real = shift0Real;
            }
            else
            {
                shift1Real = (halfTrace-absDisc)*scale;
                shift0Real = shift1Real;
            }
            shift0Imag = shift1Imag = zero;
        }
    }
}

template<typename Real>
Int ChooseDoubleShiftStart
( const Matrix<Real>& H, 
  const Real& shift0Real, const Real& shift0Imag,
  const Real& shift1Real, const Real& shift1Imag,
        vector<Real>& v )
{
    const Real ulp = limits::Precision<Real>();
    const Int n = H.Height();
    const Real* HBuf = H.LockedBuffer();
    const Int HLDim = H.LDim();

    v.resize( 3 );
    for( Int k=n-3; k>=0; --k )
    {
        const Real& eta11 = HBuf[ k   + k   *HLDim];
        const Real& eta12 = HBuf[ k   +(k+1)*HLDim];
        const Real& eta21 = HBuf[(k+1)+ k   *HLDim];
        const Real& eta22 = HBuf[(k+1)+(k+1)*HLDim];
        const Real& eta32 = HBuf[(k+2)+(k+1)*HLDim];
  
        Real scale = Abs(eta11-shift1Real) + Abs(shift1Imag) + Abs(eta21);
        Real eta21Scale = eta21 / scale;

        v[0] =
          eta21Scale*eta12 +
          (eta11-shift0Real)*((eta11-shift1Real)/scale) -
          shift0Imag*(shift1Imag/scale);
        v[1] = eta21Scale*(eta11+eta22-shift0Real-shift1Real);
        v[2] = eta21Scale*eta32;

        scale = Abs(v[0]) + Abs(v[1]) + Abs(v[2]);
        v[0] /= scale;
        v[1] /= scale;
        v[2] /= scale;
        if( k == 0 )
        {
            return k;
        }
        
        const Real& eta00 = HBuf[(k-1)+(k-1)*HLDim];
        const Real& eta10 = HBuf[ k   +(k-1)*HLDim];
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
void SingleShiftSweep
(       Matrix<Complex<Real>>& H,
        Int winBeg,
        Int winEnd,
        Complex<Real> shift,
  bool fullTriangle,
        Matrix<Complex<Real>>& Z,
  bool wantSchurVecs )
{
    typedef Complex<Real> F;
    const Int n = H.Height();
    const Int nZ = Z.Height();
    F* HBuf = H.Buffer();
    F* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    const Int transformBeg = ( fullTriangle ? 0 : winBeg ); 
    const Int transformEnd = ( fullTriangle ? n : winEnd );

    auto subInd = IR(winBeg,winEnd);
    auto qrTuple = ChooseSingleShiftStart( H(subInd,subInd), shift );
    const Int shiftStart = winBeg + std::get<0>(qrTuple);
    Complex<Real> nu0 = std::get<1>(qrTuple);
    Complex<Real> nu1 = std::get<2>(qrTuple);

    for( Int k=shiftStart; k<winEnd-1; ++k )
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
        const Int rightApplyEnd = Min(k+3,winEnd);
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
            for( Int j=0; j<nZ; ++j )
            {
                F innerProd =
                  Conj(tau0)*ZBuf[j+ k   *ZLDim] +
                       tau1 *ZBuf[j+(k+1)*ZLDim];
                ZBuf[j+ k   *ZLDim] -= innerProd;
                ZBuf[j+(k+1)*ZLDim] -= innerProd*Conj(nu1);
            }
        }

        if( k == shiftStart && shiftStart > winBeg )
        {
            // Make H[shiftStart,shiftStart-1] real by scaling by phase
            // TODO: Investigate more carefully
            F subdiagVal = Real(1) - Conj(tau0);
            F phase = subdiagVal / Abs(subdiagVal);
            HBuf[(shiftStart+1)+shiftStart*HLDim] *= Conj(phase);
            if( shiftStart+2 < winEnd )
                HBuf[(shiftStart+2)+(shiftStart+1)*HLDim] *= phase;
            for( Int j=shiftStart; j<winEnd; ++j )
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
                        blas::Scal( nZ, Conj(phase), &ZBuf[j*ZLDim], 1 );
                    }
                }
            }
        }
    }
    // Make H(winEnd-1,winEnd-2) real by scaling by phase
    F subdiagVal = HBuf[(winEnd-1)+(winEnd-2)*HLDim];
    if( ImagPart(subdiagVal) != Real(0) )
    {
        Real subdiagAbs = Abs(subdiagVal);
        HBuf[(winEnd-1)+(winEnd-2)*HLDim] = subdiagAbs;
        F phase = subdiagVal / subdiagAbs;
        if( winEnd < transformEnd ) 
        {
            blas::Scal
            ( transformEnd-winEnd, Conj(phase),
              &HBuf[(winEnd-1)+winEnd*HLDim], HLDim );
        }
        blas::Scal
        ( (winEnd-1)-transformBeg, phase,
          &HBuf[transformBeg+(winEnd-1)*HLDim], 1 );
        if( wantSchurVecs )
        {
            blas::Scal( nZ, phase, &ZBuf[(winEnd-1)*ZLDim], 1 );
        }
    }
}

template<typename Real>
void DoubleShiftSweep
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  Real shift0Real, Real shift0Imag, 
  Real shift1Real, Real shift1Imag,
  bool fullTriangle,
  Matrix<Real>& Z,
  bool wantSchurVecs )
{
    const Real zero(0), one(1);
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Real* HBuf = H.Buffer();
    Real* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    const Int transformBeg = ( fullTriangle ? 0 : winBeg ); 
    const Int transformEnd = ( fullTriangle ? n : winEnd );

    vector<Real> v(3);
    auto subInd = IR(winBeg,winEnd);
    Int shiftStart = winBeg +
      ChooseDoubleShiftStart
      ( H(subInd,subInd),
        shift0Real, shift0Imag,
        shift1Real, shift1Imag,
        v );

    for( Int k=shiftStart; k<winEnd-1; ++k )
    {
        const Int numReflect = Min( 3, winEnd-k );
        if( k > shiftStart )
        {
            MemCopy( v.data(), &HBuf[k+(k-1)*HLDim], numReflect );
        }
        Real tau0 = lapack::Reflector( numReflect, v[0], &v[1], 1 );
        if( k > shiftStart )
        {
            HBuf[ k   +(k-1)*HLDim] = v[0];
            HBuf[(k+1)+(k-1)*HLDim] = zero;
            if( k < winEnd-2 )
                HBuf[(k+2)+(k-1)*HLDim] = zero;
        }
        else if( shiftStart > winBeg )
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
            const Int rightApplyEnd = Min(k+4,winEnd);
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
                for( Int j=0; j<nZ; ++j )
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
            const Int rightApplyEnd = Min(k+3,winEnd);
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
                for( Int j=0; j<nZ; ++j )
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
void ImplicitQQuadraticSeed
( const Int n,
  const Real* HBuf, const Int HLDim,
  const Real& shift0Real, const Real& shift0Imag, 
  const Real& shift1Real, const Real& shift1Imag,
        Real* v )
{
    const Real zero(0);
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
        const Real& eta00 = HBuf[0+0*HLDim];
        const Real& eta01 = HBuf[0+1*HLDim];
        const Real& eta10 = HBuf[1+0*HLDim];
        const Real& eta11 = HBuf[1+1*HLDim];

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
        const Real& eta00 = HBuf[0+0*HLDim];
        const Real& eta01 = HBuf[0+1*HLDim];
        const Real& eta02 = HBuf[0+2*HLDim];
        const Real& eta10 = HBuf[1+0*HLDim];
        const Real& eta11 = HBuf[1+1*HLDim];
        const Real& eta12 = HBuf[1+2*HLDim];
        const Real& eta20 = HBuf[2+0*HLDim];
        const Real& eta21 = HBuf[2+1*HLDim]; 
        const Real& eta22 = HBuf[2+2*HLDim];

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
void PairShifts( vector<Real>& realShifts, vector<Real>& imagShifts )
{
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
    const Real zero(0), one(1);
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Real* HBuf = H.Buffer();
    Real* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

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
    HBuf[(winBeg+2)+winBeg] = zero;

    // Set aside space for storing either the three nonzero entries of the first
    // column of a quadratic polynomial for introducing each bulge or the scalar
    // 'tau' and non-unit subvector v such that the 2x2 or 3x3 Householder
    // reflection I - tau*[1;v] [1;v]' can safely be stored within a vector of
    // length 3 (following the LAPACK convention that 'tau' will lie in the
    // first entry and 'v' will begin in the second entry).
    W.Resize( 3, numBulges );
    Real* WBuf = W.Buffer();
    const Int WLDim = W.LDim();

    // Set aside space for a Householder vector for a candidate reinflation of
    // a deflated bulge
    vector<Real> wCand(3);

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

            if( have2x2 )
            {
                const Int bulge = lastFull+1;
                const Int k = packetBeg + 3*bulge;
                Real* w = &WBuf[bulge*WLDim];

                if( k == winBeg-1 )
                {
                    // Introduce a 2x2 bulge
                    ImplicitQQuadraticSeed
                    ( 2, &HBuf[(k+1)+(k+1)*HLDim], HLDim,
                      realShifts[2*bulge],   imagShifts[2*bulge],
                      realShifts[2*bulge+1], imagShifts[2*bulge+1],
                      w ); 
                    Real beta = w[0];
                    w[0] = lapack::Reflector( 2, beta, &w[1], 1 );
                }
                else
                {
                    // Find the reflection for chasing the 2x2 bulge
                    Real beta = HBuf[(k+1)+k*HLDim];
                    w[1] = HBuf[(k+2)+k*HLDim]; 
                    w[0] = lapack::Reflector( 2, beta, &w[1], 1 );
                    HBuf[(k+1)+k*HLDim] = beta;
                    HBuf[(k+2)+k*HLDim] = zero;
                }
            }
            for( Int bulge=lastFull; bulge>=firstFull; --bulge )
            {
                const Int k = packetBeg + 3*bulge;
                Real* w = &WBuf[bulge*WLDim];

                if( k == winBeg-1 )
                {
                    // Introduce the bulge into the upper-left corner
                    ImplicitQQuadraticSeed
                    ( 3, &HBuf[(k+1)+(k+1)*HLDim], HLDim,
                      realShifts[2*bulge],   imagShifts[2*bulge],
                      realShifts[2*bulge+1], imagShifts[2*bulge+1],
                      w ); 

                    Real alpha = w[0];
                    w[0] = lapack::Reflector( 3, alpha, &w[1], 1 );
                }
                else
                {
                    // Prepare to chase the bulge down a step
                    Real& eta_k_k     = HBuf[ k   + k   *HLDim];
                    Real& eta_kp1_k   = HBuf[(k+1)+ k   *HLDim];
                    Real& eta_kp1_kp1 = HBuf[(k+1)+(k+1)*HLDim];
                    Real& eta_kp2_k   = HBuf[(k+2)+ k   *HLDim];
                    Real& eta_kp2_kp2 = HBuf[(k+2)+(k+2)*HLDim];
                    Real& eta_kp3_k   = HBuf[(k+3)+ k   *HLDim];
                    Real& eta_kp3_kp1 = HBuf[(k+3)+(k+1)*HLDim];
                    Real& eta_kp3_kp2 = HBuf[(k+3)+(k+2)*HLDim];

                    Real beta = eta_kp1_k;
                    w[1]  = eta_kp2_k;
                    w[2]  = eta_kp3_k;
                    w[0] = lapack::Reflector( 3, beta, &w[1], 1 );
                    
                    // "Vigilantly" search for a deflation
                    // (deflation within the interior is exceedingly rare,
                    //  but this check should be essentially free)
                    if( eta_kp3_k   != zero ||
                        eta_kp3_kp1 != zero ||
                        eta_kp3_kp2 == zero )
                    {
                        // Bulge has not collapsed
                        eta_kp1_k = beta;
                        eta_kp2_k = zero;
                        eta_kp3_k = zero;
                    }
                    else
                    {
                        // Bulge has collapsed
                        ImplicitQQuadraticSeed
                        ( 3, &HBuf[(k+1)+(k+1)*HLDim], HLDim,
                          realShifts[2*bulge],   imagShifts[2*bulge],
                          realShifts[2*bulge+1], imagShifts[2*bulge+1],
                          wCand.data() ); 
                        Real alpha = wCand[0];
                        wCand[0] = lapack::Reflector( 3, alpha, &wCand[1], 1 );
                        Real innerProd =
                          wCand[0]*(eta_kp1_k+wCand[1]*eta_kp2_k);
                        if( Abs(eta_kp2_k-innerProd*wCand[1]) +
                            Abs(innerProd*wCand[2]) >=
                            ulp*(Abs(eta_k_k)+
                                 Abs(eta_kp1_kp1)+
                                 Abs(eta_kp2_kp2)) )
                        {
                            // The proposed bulge was unacceptable;
                            // continue using the collapsed one with regret
                            eta_kp1_k = beta;
                            eta_kp2_k = zero;
                            eta_kp3_k = zero;
                        }
                        else
                        {
                            // Accept the proposed replacement bulge
                            eta_kp1_k -= innerProd;
                            eta_kp2_k = zero; 
                            eta_kp3_k = zero;
                            w[0] = wCand[0];
                            w[1] = wCand[1];
                            w[2] = wCand[2];
                        }
                    }
                }
            }

            // Apply reflections from the left
            Int transformEnd;
            if( accumulate )
                transformEnd = Min( slabEnd, winEnd );
            else if( fullTriangle )
                transformEnd = n;
            else
                transformEnd = winEnd;
            if( have2x2 )
            {
                const Int bulge = lastFull+1;
                const Int k = packetBeg + 3*bulge;
                const Real* w = &WBuf[bulge*WLDim];

                for( Int j=Max(k+1,winBeg); j<transformEnd; ++j )
                {
                    Real& eta_kp1_j = HBuf[(k+1)+j*HLDim];
                    Real& eta_kp2_j = HBuf[(k+2)+j*HLDim];

                    // Apply I - tau*[1;v]*[1;v]', where tau is stored in
                    // w[0] and v is stored in w[1]
                    Real innerProd = w[0]*(eta_kp1_j+w[1]*eta_kp2_j);
                    eta_kp1_j -=      innerProd;
                    eta_kp2_j -= w[1]*innerProd;
                }
            }
            for( Int j=Max(packetBeg,winBeg); j<transformEnd; ++j )
            {
                const Int lastBulge = Min( lastFull, (j-packetBeg-1)/3 );
                for( Int bulge=lastFull; bulge>=firstFull; --bulge )
                {
                    const Int k = packetBeg + 3*bulge;
                    const Real* w = &WBuf[bulge*WLDim];
                    Real& eta_kp1_j = HBuf[(k+1)+j*HLDim]; 
                    Real& eta_kp2_j = HBuf[(k+2)+j*HLDim];
                    Real& eta_kp3_j = HBuf[(k+3)+j*HLDim];

                    // Apply I - tau*[1;v]*[1;v]', where tau is stored in
                    // w[0] and v is stored in w[1] and w[2]
                    Real innerProd =
                      w[0]*(eta_kp1_j+w[1]*eta_kp2_j+w[2]*eta_kp3_j);
                    eta_kp1_j -=      innerProd;
                    eta_kp2_j -= w[1]*innerProd;
                    eta_kp3_j -= w[2]*innerProd;
                }
            }

            // Apply reflections from the right
            Int transformBeg;
            if( accumulate )
                transformBeg = Max( winBeg, ghostBeg );
            else if( fullTriangle )
                transformBeg = 0;
            else
                transformBeg = winBeg;
            if( have2x2 )
            {
                const Int bulge = lastFull+1;
                const Int k = packetBeg + 3*bulge;
                const Real* w = &WBuf[bulge*WLDim];

                for( Int j=transformBeg; j<Min(winEnd,k+3); ++j )
                {
                    Real& eta_j_kp1 = HBuf[j+(k+1)*HLDim];
                    Real& eta_j_kp2 = HBuf[j+(k+2)*HLDim];

                    Real innerProd = w[0]*(eta_j_kp1+eta_j_kp2*w[1]);
                    eta_j_kp1 -= innerProd;
                    eta_j_kp2 -= innerProd*w[1];
                }

                if( accumulate )
                {
                    Real* UBuf = U.Buffer();
                    const Int ULDim = U.LDim();

                    Int kU = k - ghostBeg - 1; // TODO: Check for off-by-one
                    for( Int j=Max(0,(winBeg-1)-ghostBeg); j<slabSize; ++j ) 
                    {
                        Real& ups_j_kUp1 = UBuf[j+(kU+1)*ULDim];
                        Real& ups_j_kUp2 = UBuf[j+(kU+2)*ULDim];

                        Real innerProd = w[0]*(ups_j_kUp1+ups_j_kUp2*w[1]);
                        ups_j_kUp1 -= innerProd;
                        ups_j_kUp2 -= innerProd*w[1];
                    }
                }
                else if( wantSchurVecs )
                {
                    for( Int j=0; j<nZ; ++j )
                    {
                        Real& zeta_j_kp1 = ZBuf[j+(k+1)*ZLDim];
                        Real& zeta_j_kp2 = ZBuf[j+(k+2)*ZLDim];

                        Real innerProd = w[0]*(zeta_j_kp1+zeta_j_kp2*w[1]);
                        zeta_j_kp1 -= innerProd;
                        zeta_j_kp2 -= innerProd*w[1];
                    }
                }
            }
            for( Int bulge=lastFull; bulge>=firstFull; --bulge )
            {
                const Int k = packetBeg + 3*bulge;
                const Real* w = &WBuf[bulge*WLDim];

                for( Int j=transformBeg; j<Min(winEnd,k+4); ++j )
                {
                    Real& eta_j_kp1 = HBuf[j+(k+1)*HLDim];
                    Real& eta_j_kp2 = HBuf[j+(k+2)*HLDim];
                    Real& eta_j_kp3 = HBuf[j+(k+3)*HLDim];

                    Real innerProd =
                      w[0]*(eta_j_kp1+eta_j_kp2*w[1]+eta_j_kp3*w[2]);
                    eta_j_kp1 -= innerProd;
                    eta_j_kp2 -= innerProd*w[1];
                    eta_j_kp3 -= innerProd*w[2];
                }

                if( accumulate )
                {
                    Real* UBuf = U.Buffer();
                    const Int ULDim = U.LDim();

                    Int kU = k - ghostBeg - 1; // TODO: Check for off-by-one
                    for( Int j=Max(0,(winBeg-1)-ghostBeg); j<slabSize; ++j ) 
                    {
                        Real& ups_j_kUp1 = UBuf[j+(kU+1)*ULDim];
                        Real& ups_j_kUp2 = UBuf[j+(kU+2)*ULDim];
                        Real& ups_j_kUp3 = UBuf[j+(kU+3)*ULDim];

                        Real innerProd =
                          w[0]*(ups_j_kUp1+ups_j_kUp2*w[1]+ups_j_kUp3*w[2]);
                        ups_j_kUp1 -= innerProd;
                        ups_j_kUp2 -= innerProd*w[1];
                        ups_j_kUp3 -= innerProd*w[2];
                    }
                }
                else if( wantSchurVecs )
                {
                    for( Int j=0; j<nZ; ++j )
                    {
                        Real& zeta_j_kp1 = ZBuf[j+(k+1)*ZLDim];
                        Real& zeta_j_kp2 = ZBuf[j+(k+2)*ZLDim];
                        Real& zeta_j_kp3 = ZBuf[j+(k+3)*ZLDim];

                        Real innerProd =
                          w[0]*(zeta_j_kp1+zeta_j_kp2*w[1]+zeta_j_kp3*w[2]);
                        zeta_j_kp1 -= innerProd;
                        zeta_j_kp2 -= innerProd*w[1];
                        zeta_j_kp3 -= innerProd*w[2];
                    }
                }
            }

            // Vigilant deflation check using Ahues/Tisseur criteria
            const Int vigFirst =
              ( packetBeg+3*firstFull == winBeg-1 ? firstFull+1 : firstFull );
            const Int vigLast = ( have2x2 ? lastFull+1 : lastFull );
            // NOTE: LAPACK has a strange increment of vigLast by one if
            //       packetBeg is equal to winBeg-3 here...
            for( Int bulge=vigLast; bulge>=vigFirst; --bulge )
            {
                const Int k = Min( winEnd-2, packetBeg+3*bulge );
                Real& eta00 = HBuf[ k   + k   *HLDim];
                Real& eta01 = HBuf[ k   +(k+1)*HLDim];
                Real& eta10 = HBuf[(k+1)+ k   *HLDim];
                Real& eta11 = HBuf[(k+1)+(k+1)*HLDim];

                if( eta10 == zero )
                    continue;
                Real localScale = Abs(eta00) + Abs(eta11);
                if( localScale == zero )
                {
                    if( k >= winBeg+3 )
                        localScale += Abs(HBuf[k+(k-3)*HLDim]);
                    if( k >= winBeg+2 )
                        localScale += Abs(HBuf[k+(k-2)*HLDim]);
                    if( k >= winBeg+1 )
                        localScale += Abs(HBuf[k+(k-1)*HLDim]);
                    if( k < winEnd-2 )
                        localScale += Abs(HBuf[(k+2)+(k+1)*HLDim]);
                    if( k < winEnd-3 )
                        localScale += Abs(HBuf[(k+3)+(k+1)*HLDim]);
                    if( k < winEnd-4 )
                        localScale += Abs(HBuf[(k+4)+(k+1)*HLDim]);
                }
                if( Abs(eta10) <= Max(smallNum,ulp*localScale) )
                {
                    Real offMax = Max( Abs(eta10), Abs(eta01) );
                    Real offMin = Min( Abs(eta10), Abs(eta01) ); 
                    Real diagMax = Max( Abs(eta11), Abs(eta00-eta11) );
                    Real diagMin = Min( Abs(eta11), Abs(eta00-eta11) );
                    Real scale = diagMax + offMax;
                    Real localScale2 = diagMin*(diagMax/scale);
                    if( localScale2 == zero ||
                        offMin*(offMax/scale) <= Max(smallNum,ulp*localScale2) )
                    {
                        eta10 = zero;
                    }
                }
            }

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
                const Real* w = &WBuf[bulge*WLDim];
                Real& eta_kp4_kp1 = HBuf[(k+4)+(k+1)*HLDim];
                Real& eta_kp4_kp2 = HBuf[(k+4)+(k+2)*HLDim]; 
                Real& eta_kp4_kp3 = HBuf[(k+4)+(k+3)*HLDim];

                const Real innerProd = w[0]*w[2]*eta_kp4_kp3;
                eta_kp4_kp1 = -innerProd;
                eta_kp4_kp2 = -innerProd*w[1];
                eta_kp4_kp3 -= innerProd*w[2];
            }
        }
        if( accumulate )
        {
            Real* UBuf = U.Buffer();
            const Int ULDim = U.LDim();
            const Int transformBeg = ( fullTriangle ? 0 : winBeg );
            const Int transformEnd = ( fullTriangle ? n : winEnd );
 
            const Int k1 = Max( 0, winBeg-(ghostBeg+1) );
            const Int nU = ((slabSize-1)-Max(0,slabEnd-winEnd)) - k1;
            
            const Int j = Min(slabEnd,winEnd);

            // Horizontal far-from-diagonal application
            const Int horzSize = transformEnd-transformBeg;
            Zeros( WAccum, nU, horzSize );
            blas::Gemm
            ( 'C', 'N', nU, horzSize, nU, 
              Real(1), &UBuf[k1+k1*ULDim],           ULDim,
                       &HBuf[(ghostBeg+k1)+j*HLDim], HLDim,
              Real(0), WAccum.Buffer(), WAccum.LDim() );
            lapack::Copy
            ( 'A', nU, horzSize,
              WAccum.Buffer(), WAccum.LDim(),
              &HBuf[(ghostBeg+k1)+j*HLDim], HLDim );

            // Vertical far-from-diagonal application
            const Int vertSize = Max( winBeg, ghostBeg ) - transformBeg; 
            Zeros( WAccum, vertSize, nU );
            blas::Gemm
            ( 'N', 'N', vertSize, nU, nU, 
              Real(1), &HBuf[transformBeg+(ghostBeg+k1)*HLDim], HLDim,
                       &UBuf[k1+k1*ULDim],                      ULDim,
              Real(0), WAccum.Buffer(), WAccum.LDim() );
            lapack::Copy
            ( 'A', vertSize, nU,
              WAccum.Buffer(), WAccum.LDim(),
              &HBuf[transformBeg+(ghostBeg+k1)*HLDim], HLDim );

            if( wantSchurVecs )
            {
                Zeros( WAccum, nZ, nU );
                blas::Gemm
                ( 'N', 'N', nZ, nU, nU,
                  Real(1), &ZBuf[(ghostBeg+k1)*ZLDim], ZLDim,
                           &UBuf[k1+k1*ULDim],         ULDim,
                  Real(0), WAccum.Buffer(), WAccum.LDim() );
                blas::Copy
                ( 'A', nZ, nU,
                  WAccum.Buffer(), WAccum.LDim(),
                  &ZBuf[(ghostBeg+k1)*ZLDim], ZLDim );
            }
        }
    }
}

template<typename Real>
Int WindowedSingle
( Matrix<Real>& H,
  Int winBeg,
  Int winEnd,
  Matrix<Real>& wReal,
  Matrix<Real>& wImag,
  bool fullTriangle,
  Matrix<Real>& Z,
  bool wantSchurVecs,
  bool demandConverged )
{
    const Real zero(0), threeFourths=Real(3)/Real(4);
    const Int maxIter=30;    
    // Cf. LAPACK for these somewhat arbitrary constants
    const Real exceptScale0=Real(3)/Real(4),
               exceptScale1=Real(-4375)/Real(10000);
    const Int n = H.Height();
    const Int windowSize = winEnd - winBeg;
    const Int nZ = Z.Height();
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
        return winBeg;
    }
    if( windowSize == 1 )
    {
        wRealBuf[winBeg] = HBuf[winBeg+winBeg*HLDim];
        wImagBuf[winBeg] = zero;
        return winBeg;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=winBeg; j<winEnd-3; ++j ) 
    {
        HBuf[(j+2)+j*HLDim] = 0;
        HBuf[(j+3)+j*HLDim] = 0;
    }
    if( winBeg <= winEnd-3 )
        HBuf[(winEnd-1)+(winEnd-3)*HLDim] = 0;
    
    // Attempt to converge the eigenvalues one or two at a time
    vector<Real> v; // for computing Householder reflectors
    while( winBeg < winEnd )
    {
        Int iterBeg = winBeg;
        Int iter;
        for( iter=0; iter<maxIter; ++iter )
        {
            {
                auto winInd = IR(iterBeg,winEnd);
                iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
            }
            if( iterBeg > winBeg )
            {
                HBuf[iterBeg+(iterBeg-1)*HLDim] = zero;
            }
            if( iterBeg == winEnd-1 )
            {
                wRealBuf[iterBeg] = HBuf[iterBeg+iterBeg*HLDim];
                wImagBuf[iterBeg] = zero;
                --winEnd;
                break;
            }
            else if( iterBeg == winEnd-2 )
            {
                Real c, s;
                Real& alpha00 = HBuf[(winEnd-2)+(winEnd-2)*HLDim];
                Real& alpha01 = HBuf[(winEnd-2)+(winEnd-1)*HLDim];
                Real& alpha10 = HBuf[(winEnd-1)+(winEnd-2)*HLDim];
                Real& alpha11 = HBuf[(winEnd-1)+(winEnd-1)*HLDim];
                lapack::TwoByTwoSchur
                ( alpha00, alpha01,
                  alpha10, alpha11,
                  c, s,
                  wRealBuf[iterBeg],   wImagBuf[iterBeg], 
                  wRealBuf[iterBeg+1], wImagBuf[iterBeg+1] );
                if( fullTriangle )
                {
                    if( n > winEnd )
                        blas::Rot
                        ( n-winEnd,
                          &HBuf[(winEnd-2)+winEnd*HLDim], HLDim,
                          &HBuf[(winEnd-1)+winEnd*HLDim], HLDim,
                          c, s );
                    blas::Rot
                    ( winEnd-2,
                      &HBuf[0+(winEnd-2)*HLDim], 1,
                      &HBuf[0+(winEnd-1)*HLDim], 1,
                      c, s );
                }
                if( wantSchurVecs )
                {
                    blas::Rot
                    ( nZ,
                      &ZBuf[(winEnd-2)*ZLDim], 1,
                      &ZBuf[(winEnd-1)*ZLDim], 1,
                      c, s );
                }
                winEnd -= 2;
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
                const Real* HSub = &HBuf[(winEnd-3)+(winEnd-3)*HLDim];
                const Real scale = Abs(HSub[2+1*HLDim]) + Abs(HSub[1+0*HLDim]);
                eta00 = exceptScale0*scale + HSub[2+2*HLDim];
                eta11 = exceptScale1*scale;
                eta10 = scale;
                eta11 = eta00;
            }
            else
            {
                const Real* HSub = &HBuf[(winEnd-2)+(winEnd-2)*HLDim];
                eta00 = HSub[0+0*HLDim];
                eta01 = HSub[0+1*HLDim];
                eta10 = HSub[1+0*HLDim];
                eta11 = HSub[1+1*HLDim];
            }
            Real shift0Real, shift0Imag, shift1Real, shift1Imag;
            PrepareDoubleShift
            ( eta00, eta01,
              eta10, eta11,
              shift0Real, shift0Imag,
              shift1Real, shift1Imag );

            auto subInd = IR(iterBeg,winEnd);
            DoubleShiftSweep
            ( H,
              iterBeg, winEnd,
              shift0Real, shift0Imag,
              shift1Real, shift1Imag,
              fullTriangle,
              Z,
              wantSchurVecs );
        }
        if( iter == maxIter && demandConverged )
            RuntimeError("QR iteration did not converge");
    }
    return winEnd;
}

template<typename Real>
Int WindowedSingle
(       Matrix<Complex<Real>>& H,
        Int winBeg,
        Int winEnd,
        Matrix<Complex<Real>>& w,
  bool fullTriangle,
        Matrix<Complex<Real>>& Z,
  bool wantSchurVecs,
  bool demandConverged )
{
    typedef Complex<Real> F;
    const Real zero(0), threeFourths=Real(3)/Real(4);
    const Int maxIter=30;    
    const Int n = H.Height();
    const Int windowSize = winEnd - winBeg;
    const Int nZ = Z.Height();
    F* HBuf = H.Buffer(); 
    F* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    w.Resize( n, 1 );
    F* wBuf = w.Buffer();

    if( windowSize == 0 )
    {
        return winBeg;
    }
    if( windowSize == 1 )
    {
        wBuf[winBeg] = HBuf[winBeg+winBeg*HLDim];
        return winBeg;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=winBeg; j<winEnd-3; ++j ) 
    {
        HBuf[(j+2)+j*HLDim] = 0;
        HBuf[(j+3)+j*HLDim] = 0;
    }
    if( winBeg <= winEnd-3 )
        HBuf[(winEnd-1)+(winEnd-3)*HLDim] = 0;
    
    // Rotate the matrix so that the subdiagonals are real
    const Int scaleBeg = ( fullTriangle ? 0 : winBeg );
    const Int scaleEnd = ( fullTriangle ? n : winEnd );
    for( Int i=winBeg+1; i<winEnd; ++i )
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
                blas::Scal( nZ, Conj(phase), &ZBuf[i*ZLDim], 1 );
        }
    }

    // Attempt to converge the eigenvalues one at a time
    while( winBeg < winEnd )
    {
        Int iterBeg = winBeg;
        Int iter;
        for( iter=0; iter<maxIter; ++iter )
        {
            {
                auto winInd = IR(iterBeg,winEnd);
                iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
            }
            if( iterBeg > winBeg )
            {
                HBuf[iterBeg+(iterBeg-1)*HLDim] = zero;
            }
            if( iterBeg == winEnd-1 )
            {
                wBuf[iterBeg] = HBuf[iterBeg+iterBeg*HLDim];
                --winEnd;
                break;
            }

            // Pick either the Wilkinson shift or an exceptional shift
            F shift;
            auto subInd = IR(iterBeg,winEnd);
            if( iter == maxIter/3 )
            {
                const F diagVal = HBuf[iterBeg+iterBeg*HLDim];
                const Real subdiagVal =
                  RealPart(HBuf[(iterBeg+1)+iterBeg*HLDim]);
                shift = diagVal + threeFourths*Abs(subdiagVal);
            } 
            else if( iter == 2*maxIter/3 )
            {
                const F diagVal = HBuf[(winEnd-1)+(winEnd-1)*HLDim];
                const Real subdiagVal =
                  RealPart(HBuf[(winEnd-1)+(winEnd-2)*HLDim]);
                shift = diagVal + threeFourths*Abs(subdiagVal);
            }
            else
            {
                shift = WilkinsonShift( H(subInd,subInd) );
            }

            SingleShiftSweep
            ( H,
              iterBeg, winEnd,
              shift,
              fullTriangle,
              Z,
              wantSchurVecs );
        }
        if( iter == maxIter && demandConverged )
            RuntimeError("QR iteration did not converge");
    }
    return winEnd;
}

// The best references are
//
//   Karen Braman, Ralph Byers, and Roy Mathias,
//   "The multishift QR algorithm. Part II: Aggressive Early Deflation",
//   SIAM J. Matrix Anal. Appl., Vol. 23, No. 4, pp. 948--973, 2002
//
// and the LAPACK implementation DLAQR2, which has several distinct differences
// from the suggestions of Braman et al., such as:
//
//   1) Solely using "nearby-diagonal deflation" instead of Braman et al.'s 
//      suggestion of also allowing for "window-Schur deflation".
//
//   2) Using the largest (in magnitude) eigenvalue of a 2x2 Schur block to 
//      determine whether it qualifies for "nearby-diagonal deflation" rather
//      that using the square-root of the absolute value of the determinant
//      (which would correspond to the geometric mean of the eigenvalue
//       magnitudes). 
//
// In both respects, the LAPACK implementation is significantly more
// conservative than the original suggestions of Braman et al.
//
template<typename Real>
Int SpikeDeflation
(       Matrix<Real>& T,
        Matrix<Real>& V,
  const Real& eta,
        Int winBeg=0 )
{
    DEBUG_ONLY(CSE cse("hess_qr::SpikeDeflation"))

    const Int n = T.Height();
    Real* TBuf = T.Buffer();
    Real* VBuf = V.Buffer();
    const Int TLDim = T.LDim();
    const Int VLDim = V.LDim();

    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    // TODO: Take this as an argument to SpikeDeflation to avoid reallocation?
    vector<Real> work(n);

    Int winEnd = n;
    while( winBeg < winEnd )
    {
        const bool twoByTwo =
          ( winEnd==1 ? false : TBuf[(winEnd-1)+(winEnd-2)*TLDim] != zero );

        if( twoByTwo )
        {
            // Follow LAPACK's suit (rather than Braman et al.) and use the 
            // eigenvalue of the 2x2 with largest magnitude in order to 
            // determine if this entry of the spike qualifies for 
            // "nearby-diagonal" deflation. Recall that the 2x2 block is assumed
            // to be in standard form,
            //
            //    | alpha, gamma |, where beta*gamma < 0,
            //    | beta,  alpha |
            //
            // so that the eigenvalues are alpha +- sqrt(beta*gamma) and the
            // spectral radius is |alpha| + sqrt(beta*gamma).
            //
            const Real& alpha = TBuf[(winEnd-2)+(winEnd-2)*TLDim];
            const Real& beta  = TBuf[(winEnd-1)+(winEnd-2)*TLDim]; 
            const Real& gamma = TBuf[(winEnd-2)+(winEnd-1)*TLDim];
            const Real spectralRadius =
                Abs(alpha) + Sqrt(Abs(beta))*Sqrt(Abs(gamma));
            const Real scale =
              ( spectralRadius > 0 ? spectralRadius : Abs(eta) );

            // The relevant two entries of spike V^T [eta; zeros(n-1,1)]
            const Real sigma0 = eta*VBuf[(winEnd-2)*VLDim];
            const Real sigma1 = eta*VBuf[(winEnd-1)*VLDim];

            if( Max( Abs(sigma0), Abs(sigma1) ) <= Max( smallNum, ulp*scale ) )
            {
                // The two-by-two block satisfies the "nearby-diagonal" test
                winEnd -= 2;
            }
            else
            {
                // Move this undeflatable 2x2 block to the top of the window 
                lapack::SchurExchange
                ( n, T, TLDim, V, VLDim, winEnd-2, winBeg, work.data() );
                winBeg += 2;
            }
        }
        else
        {
            const Real spectralRadius = Abs(TBuf[(winEnd-1)+(winEnd-1)*TLDim]);
            const Real scale =
              ( spectralRadius > 0 ? spectralRadius : Abs(eta) );

            // The relevant entry of the spike V^T [eta; zeros(n-1,1)]
            const Real sigma = eta*VBuf[(winEnd-1)*VLDim];

            if( Abs(sigma) <= Max( smallNum, ulp*scale ) )
            {
                // The one-by-one block satisfies the "nearby-diagonal" test
                winEnd -= 1;
            }
            else
            {
                // Move the undeflatable 1x1 block to the top of the window
                lapack::SchurExchange
                ( n, T, TLDim, V, VLDim, winEnd-1, winBeg, work.data() );
                winBeg += 1;
            }
        }
    }
    // Return the number of unconverged/undeflated eigenvalues
    return winBeg;
}

} // namespace hess_qr

template<typename Real>
void HessenbergQR
( Matrix<Real>& H,
  Matrix<Real>& wReal,
  Matrix<Real>& wImag,
  bool fullTriangle=true )
{
    Int winBeg=0, winEnd=H.Height();
    bool wantSchurVecs=false;
    bool demandConverged=true;
    Matrix<Real> Z;
    hess_qr::WindowedSingle
    ( H,
      winBeg, winEnd,
      wReal, wImag,
      fullTriangle, 
      Z,
      wantSchurVecs,
      demandConverged );
}

template<typename Real>
void HessenbergQR
( Matrix<Real>& H,
  Matrix<Real>& wReal,
  Matrix<Real>& wImag,
  Matrix<Real>& Z,
  bool fullTriangle=true )
{
    Int winBeg=0, winEnd=H.Height();
    bool wantSchurVecs=true;
    bool demandConverged=true;
    hess_qr::WindowedSingle
    ( H,
      winBeg, winEnd,
      wReal, wImag,
      fullTriangle, 
      Z,
      wantSchurVecs,
      demandConverged );
}

template<typename Real>
void HessenbergQR
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& w,
  bool fullTriangle=true )
{
    Int winBeg=0, winEnd=H.Height();
    bool wantSchurVecs=false;
    bool demandConverged=true;
    Matrix<Complex<Real>> Z;
    hess_qr::WindowedSingle
    ( H,
      winBeg, winEnd,
      w,
      fullTriangle, 
      Z,
      wantSchurVecs,
      demandConverged );
}

template<typename Real>
void HessenbergQR
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& Z,
  bool fullTriangle=true )
{
    Int winBeg=0, winEnd=H.Height();
    bool wantSchurVecs=true;
    bool demandConverged=true;
    hess_qr::WindowedSingle
    ( H,
      winBeg, winEnd,
      w,
      fullTriangle, 
      Z,
      wantSchurVecs,
      demandConverged );
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
