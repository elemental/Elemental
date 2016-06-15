/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace hess_qr {

// Use Ahues and Tissuer's (LAPACK Working Note 122, 1997) refinement of
// Wilkinson's criteria for determining if a subdiagonal entry of a Hessenberg
// matrix is negligible
template<typename F>
Int DetectSmallSubdiagonal( const Matrix<F>& H )
{
    typedef Base<F> Real;

    const Int n = H.Height();
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    // Search up the subdiagonal
    for( Int k=n-2; k>=0; --k )
    {
        const F eta00 = H(k,k);
        const F eta01 = H(k,k+1);
        const Real eta10 = RealPart( H(k+1,k) );
        const F eta11 = H(k+1,k+1);
        if( OneAbs(eta10) <= smallNum )
        {
            return k+1;
        }

        Real localScale = OneAbs(eta00) + OneAbs(eta11);
        if( localScale == Real(0) )
        {
            // Search outward a bit to get a sense of the matrix's local scale
            if( k-1 >= 0 )
                localScale += Abs( RealPart( H(k,k-1) ) );
            if( k+2 <= n-1 )
                localScale += Abs( RealPart( H(k+2,k+1) ) );
        }
        
        if( Abs(eta10) <= ulp*localScale )
        {
            const Real maxOff = Max( Abs(eta10), OneAbs(eta01) );
            const Real minOff = Min( Abs(eta10), OneAbs(eta01) ); 

            const Real diagDiff = OneAbs(eta00-eta11);
            const Real maxDiag = Max( OneAbs(eta11), diagDiff );
            const Real minDiag = Min( OneAbs(eta11), diagDiff );

            const Real sigma = maxDiag + maxOff;
            if( minOff*(maxOff/sigma) <=
                Max(smallNum,ulp*(minDiag*(maxDiag/sigma))) )
            {
                return k+1;
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
    const Real zero = Real(0);

    const F& eta00 = H(offset,offset);
    const F& eta01 = H(offset,offset+1);
    const F& eta10 = H(offset+1,offset);
    const F& eta11 = H(offset+1,offset+1);
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
    for( Int k=n-3; k>=0; --k )
    {
        const Real eta10 = RealPart( H(k+1,k) );
        const F& eta11 = H(k+1,k+1);
        const F& eta22 = H(k+2,k+2);
        Real eta21 = RealPart( H(k+2,k+1) );

        F eta11Shift = eta11 - shift;
        const Real sigma = OneAbs(eta11Shift) + Abs(eta21);    
        eta11Shift /= sigma;
        eta21 /= sigma; 
        if( Abs(eta10)*Abs(eta21) <=
            ulp*(OneAbs(eta11Shift)*(OneAbs(eta11)+OneAbs(eta22))) ) 
        {
            return std::make_tuple(k+1,eta11Shift,Complex<Real>(eta21));
        }
    }
    // Start at the top
    const F& eta11 = H(0,0);
    Real eta21 = RealPart( H(1,0) );

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

    v.resize( 3 );
    for( Int k=n-4; k+1>=0; --k )
    {
        const Real& eta11 = H(k+1,k+1);
        const Real& eta12 = H(k+1,k+2);
        const Real& eta21 = H(k+2,k+1);
        const Real& eta22 = H(k+2,k+2);
        const Real& eta32 = H(k+3,k+2);
  
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
        if( k == -1 )
        {
            return 0;
        }
        
        const Real& eta00 = H(k,  k);
        const Real& eta10 = H(k+1,k);
        if( Abs(eta10)*(Abs(v[1])+Abs(v[2])) <= 
            ulp*Abs(v[0])*(Abs(eta00)+Abs(eta11)+Abs(eta22)) )
        {
            return k+1;
        }
    }
    // This should never be reached but is to avoid compiler warnings
    return 0;
}

template<typename Real>
void SingleShiftSweep
( Matrix<Complex<Real>>& H,
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
            nu0 = H(k,k-1);
            nu1 = RealPart( H(k+1,k-1) );
        }
        // TODO: Assert nu1 is real
        F tau0 = lapack::Reflector( 2, nu0, &nu1, 1 );
        if( k > shiftStart )
        {
            H(k,  k-1) = nu0;
            H(k+1,k-1) = 0;
        }
        // The formulas within lapack::Reflector trivially imply that 
        // tau0*Conj(nu1) will be real if nu1 was real on entry to
        // lapack::Reflector (an equivalent claim is made within ZLAHQR)
        F tau1 = RealPart(tau0*Conj(nu1));

        // Apply the Householder reflector from the left
        for( Int j=k; j<transformEnd; ++j )
        {
            F innerProd = tau0*H(k,j) + tau1*H(k+1,j);
            H(k,  j) -= innerProd;
            H(k+1,j) -= innerProd*nu1;
        }

        // Apply the Householder reflector from the right
        const Int rightApplyEnd = Min(k+3,winEnd);
        for( Int j=transformBeg; j<rightApplyEnd; ++j )
        {
            F innerProd = Conj(tau0)*H(j,k) + tau1*H(j,k+1);
            H(j,k  ) -= innerProd;
            H(j,k+1) -= innerProd*Conj(nu1);
        }

        if( wantSchurVecs )
        {
            // Accumulate the Schur vectors
            for( Int j=0; j<nZ; ++j )
            {
                F innerProd = Conj(tau0)*Z(j,k) + tau1*Z(j,k+1);
                Z(j,k  ) -= innerProd;
                Z(j,k+1) -= innerProd*Conj(nu1);
            }
        }

        if( k == shiftStart && shiftStart > winBeg )
        {
            // Make H[shiftStart,shiftStart-1] real by scaling by phase
            // TODO: Investigate more carefully
            F subdiagVal = Real(1) - Conj(tau0);
            F phase = subdiagVal / Abs(subdiagVal);
            H( shiftStart+1, shiftStart ) *= Conj(phase);
            if( shiftStart+2 < winEnd )
                H( shiftStart+2, shiftStart+1 ) *= phase;
            for( Int j=shiftStart; j<winEnd; ++j )
            {
                if( j != shiftStart+1 )
                {
                    if( j+1 < transformEnd )
                    {
                        blas::Scal
                        ( transformEnd-(j+1), phase, &H(j,j+1), H.LDim() );
                    }
                    blas::Scal
                    ( j-transformBeg, Conj(phase), &H(transformBeg,j), 1 );
                    if( wantSchurVecs )
                    {
                        blas::Scal( nZ, Conj(phase), &Z(0,j), 1 );
                    }
                }
            }
        }
    }
    // Make H(winEnd-1,winEnd-2) real by scaling by phase
    F subdiagVal = H( winEnd-1, winEnd-2 );
    if( ImagPart(subdiagVal) != Real(0) )
    {
        Real subdiagAbs = Abs(subdiagVal);
        H( winEnd-1, winEnd-2 ) = subdiagAbs;
        F phase = subdiagVal / subdiagAbs;
        if( winEnd < transformEnd ) 
        {
            blas::Scal
            ( transformEnd-winEnd, Conj(phase),
              &H(winEnd-1,winEnd), H.LDim() );
        }
        blas::Scal
        ( (winEnd-1)-transformBeg, phase, &H(transformBeg,winEnd-1), 1 );
        if( wantSchurVecs )
        {
            blas::Scal( nZ, phase, &Z(0,winEnd-1), 1 );
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
            MemCopy( v.data(), &H(k,k-1), numReflect );
        }
        Real tau0 = lapack::Reflector( numReflect, v[0], &v[1], 1 );
        if( k > shiftStart )
        {
            H(k,  k-1) = v[0];
            H(k+1,k-1) = zero;
            if( k < winEnd-2 )
                H(k+2,k-1) = zero;
        }
        else if( shiftStart > winBeg )
        {
            // The following is supposedly more reliable than
            // H(k,k-1) = -H(k,k-1) when v(1) and v(2) underflow
            // (cf. LAPACK's {s,d}lahqr)
            H(k,k-1) *= (one-tau0);
        }
        Real tau1 = tau0*v[1];
        if( numReflect == 3 )
        {
            Real tau2 = tau0*v[2];

            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd = H(k,j) + v[1]*H(k+1,j) + v[2]*H(k+2,j);
                H(k,  j) -= innerProd*tau0;
                H(k+1,j) -= innerProd*tau1;
                H(k+2,j) -= innerProd*tau2;
            }
           
            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+4,winEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd = H(j,k) + v[1]*H(j,k+1) + v[2]*H(j,k+2);
                H(j,k  ) -= innerProd*tau0;
                H(j,k+1) -= innerProd*tau1;
                H(j,k+2) -= innerProd*tau2;
            }

            if( wantSchurVecs )
            {
                for( Int j=0; j<nZ; ++j )
                {
                    Real innerProd = Z(j,k) + v[1]*Z(j,k+1) + v[2]*Z(j,k+2);
                    Z(j,k  ) -= innerProd*tau0;
                    Z(j,k+1) -= innerProd*tau1;
                    Z(j,k+2) -= innerProd*tau2;
                }
            }
        }
        else if( numReflect == 2 )
        {
            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd = H(k,j) + v[1]*H(k+1,j);
                H(k,  j) -= innerProd*tau0;
                H(k+1,j) -= innerProd*tau1;
            }

            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+3,winEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd = H(j,k) + v[1]*H(j,k+1);
                H(j,k  ) -= innerProd*tau0;
                H(j,k+1) -= innerProd*tau1;
            }

            if( wantSchurVecs )
            {
                // Accumulate the Schur vectors
                for( Int j=0; j<nZ; ++j )
                {
                    Real innerProd = Z(j,k) + v[1]*Z(j,k+1);
                    Z(j,k  ) -= innerProd*tau0;
                    Z(j,k+1) -= innerProd*tau1;
                }
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
    const Real zero(0);
    const Int maxIter=30;    
    // Cf. LAPACK for these somewhat arbitrary constants
    const Real exceptScale0=Real(3)/Real(4),
               exceptScale1=Real(-4375)/Real(10000);
    const Int n = H.Height();
    const Int windowSize = winEnd - winBeg;
    const Int nZ = Z.Height();

    wReal.Resize( n, 1 );
    wImag.Resize( n, 1 );

    if( windowSize == 0 )
    {
        return winBeg;
    }
    if( windowSize == 1 )
    {
        wReal(winBeg) = H(winBeg,winBeg);
        wImag(winBeg) = zero;
        return winBeg;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=winBeg; j<winEnd-3; ++j ) 
    {
        H(j+2,j) = zero;
        H(j+3,j) = zero;
    }
    if( winBeg <= winEnd-3 )
        H(winEnd-1,winEnd-3) = zero;
    
    // Attempt to converge the eigenvalues one or two at a time
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
                H(iterBeg,iterBeg-1) = zero;
            }
            if( iterBeg == winEnd-1 )
            {
                wReal(iterBeg) = H(iterBeg,iterBeg);
                wImag(iterBeg) = zero;
                --winEnd;
                break;
            }
            else if( iterBeg == winEnd-2 )
            {
                Real c, s;
                lapack::TwoByTwoSchur
                ( H(winEnd-2,winEnd-2), H(winEnd-2,winEnd-1),
                  H(winEnd-1,winEnd-2), H(winEnd-1,winEnd-1),
                  c, s,
                  wReal(iterBeg),   wImag(iterBeg), 
                  wReal(iterBeg+1), wImag(iterBeg+1) );
                if( fullTriangle )
                {
                    if( n > winEnd )
                        blas::Rot
                        ( n-winEnd,
                          &H(winEnd-2,winEnd), H.LDim(),
                          &H(winEnd-1,winEnd), H.LDim(),
                          c, s );
                    blas::Rot
                    ( winEnd-2,
                      &H(0,winEnd-2), 1,
                      &H(0,winEnd-1), 1,
                      c, s );
                }
                if( wantSchurVecs )
                {
                    blas::Rot
                    ( nZ,
                      &Z(0,winEnd-2), 1,
                      &Z(0,winEnd-1), 1,
                      c, s );
                }
                winEnd -= 2;
                break;
            }

            // Pick either the Francis shifts or exceptional shifts
            Real eta00, eta01, eta10, eta11;
            if( iter == maxIter/3 )
            {
                const Real scale =
                  Abs(H(iterBeg+1,iterBeg)) + Abs(H(iterBeg+2,iterBeg+1));
                eta00 = exceptScale0*scale + H(iterBeg,iterBeg);
                eta01 = exceptScale1*scale;
                eta10 = scale;
                eta11 = eta00; 
            } 
            else if( iter == 2*maxIter/3 )
            {
                const Real scale = 
                  Abs(H(winEnd-1,winEnd-2)) + Abs(H(winEnd-2,winEnd-3));
                eta00 = exceptScale0*scale + H(winEnd-1,winEnd-1);
                eta01 = exceptScale1*scale;
                eta10 = scale;
                eta11 = eta00;
            }
            else
            {
                eta00 = H(winEnd-2,winEnd-2);
                eta01 = H(winEnd-2,winEnd-1);
                eta10 = H(winEnd-1,winEnd-2);
                eta11 = H(winEnd-1,winEnd-1);
            }
            Real shift0Real, shift0Imag, shift1Real, shift1Imag;
            PrepareDoubleShift
            ( eta00, eta01,
              eta10, eta11,
              shift0Real, shift0Imag,
              shift1Real, shift1Imag );

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

    w.Resize( n, 1 );

    if( windowSize == 0 )
    {
        return winBeg;
    }
    if( windowSize == 1 )
    {
        w(winBeg) = H(winBeg,winBeg);
        return winBeg;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=winBeg; j<winEnd-3; ++j ) 
    {
        H(j+2,j) = zero;
        H(j+3,j) = zero;
    }
    if( winBeg <= winEnd-3 )
        H(winEnd-1,winEnd-3) = zero;
    
    // Rotate the matrix so that the subdiagonals are real
    const Int scaleBeg = ( fullTriangle ? 0 : winBeg );
    const Int scaleEnd = ( fullTriangle ? n : winEnd );
    for( Int i=winBeg+1; i<winEnd; ++i )
    {
        F& eta = H(i,i-1);
        if( ImagPart(eta) != zero )
        {
            F phase = eta / OneAbs(eta);
            phase = Conj(phase) / Abs(phase);
            eta = Abs(eta);
            blas::Scal( scaleEnd-i, phase, &H(i,i), H.LDim() );
            blas::Scal
            ( Min(scaleEnd,i+2)-scaleBeg, Conj(phase), &H(scaleBeg,i), 1 );
            if( wantSchurVecs )
                blas::Scal( nZ, Conj(phase), &Z(0,i), 1 );
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
                H(iterBeg,iterBeg-1) = zero;
            }
            if( iterBeg == winEnd-1 )
            {
                w(iterBeg) = H(iterBeg,iterBeg);
                --winEnd;
                break;
            }

            // Pick either the Wilkinson shift or an exceptional shift
            F shift;
            if( iter == maxIter/3 )
            {
                const F diagVal = H(iterBeg,iterBeg);
                const Real subdiagVal = RealPart(H(iterBeg+1,iterBeg));
                shift = diagVal + threeFourths*Abs(subdiagVal);
            } 
            else if( iter == 2*maxIter/3 )
            {
                const F diagVal = H(winEnd-1,winEnd-1);
                const Real subdiagVal = RealPart(H(winEnd-1,winEnd-2));
                shift = diagVal + threeFourths*Abs(subdiagVal);
            }
            else
            {
                auto subInd = IR(iterBeg,winEnd);
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

} // namespace hess_qr
