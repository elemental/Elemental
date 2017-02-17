/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_SINGLE_SHIFT_SWEEP_HPP
#define EL_HESS_SCHUR_SINGLE_SHIFT_SWEEP_HPP

namespace El {
namespace hess_schur {

// Use Ahues and Tissuer's (LAPACK Working Note 122, 1997) refinement of
// Wilkinson's criteria for determining if a subdiagonal entry of a Hessenberg
// matrix is negligible
template<typename Field>
Int DetectSmallSubdiagonal( const Matrix<Field>& H )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    const Int n = H.Height();
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    // Search up the subdiagonal
    for( Int k=n-2; k>=0; --k )
    {
        const Field eta00 = H(k,k);
        const Field eta01 = H(k,k+1);
        const Field eta10 = H(k+1,k);
        const Field eta11 = H(k+1,k+1);
        if( OneAbs(eta10) <= smallNum )
        {
            return k+1;
        }

        Real localScale = OneAbs(eta00) + OneAbs(eta11);
        if( localScale == Real(0) )
        {
            // Search outward a bit to get a sense of the matrix's local scale
            if( k-1 >= 0 )
                localScale += OneAbs( H(k,k-1) );
            if( k+2 <= n-1 )
                localScale += OneAbs( H(k+2,k+1) );
        }

        if( Abs(eta10) <= ulp*localScale )
        {
            const Real maxOff = Max( OneAbs(eta10), OneAbs(eta01) );
            const Real minOff = Min( OneAbs(eta10), OneAbs(eta01) );

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

// This is an adaptation of the above that only takes in the relevant
// tridiagonal information.
template<typename Field>
Int DetectSmallSubdiagonal
( const Matrix<Field>& hMain,
  const Matrix<Field>& hSub,
  const Matrix<Field>& hSuper )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    const Int n = hMain.Height();
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    // Search up the subdiagonal
    for( Int k=n-2; k>=0; --k )
    {
        const Field eta00 = hMain(k);
        const Field eta01 = hSuper(k);
        const Field eta10 = hSub(k);
        const Field eta11 = hMain(k+1);
        if( OneAbs(eta10) <= smallNum )
        {
            return k+1;
        }

        Real localScale = OneAbs(eta00) + OneAbs(eta11);
        if( localScale == Real(0) )
        {
            // Search outward a bit to get a sense of the matrix's local scale
            if( k-1 >= 0 )
                localScale += OneAbs( hSub(k-1) );
            if( k+2 <= n-1 )
                localScale += OneAbs( hSub(k+1) );
        }

        if( Abs(eta10) <= ulp*localScale )
        {
            const Real maxOff = Max( OneAbs(eta10), OneAbs(eta01) );
            const Real minOff = Min( OneAbs(eta10), OneAbs(eta01) );

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
    EL_DEBUG_CSE
    typedef Complex<Real> Field;
    const Int n = H.Height();
    EL_DEBUG_ONLY(
      if( n < 2 )
          LogicError("WilkinsonShift requires n >= 2");
    )
    const Int offset = n-2;
    const Real zero = Real(0);

    const Field& eta00 = H(offset,offset);
    const Field& eta01 = H(offset,offset+1);
    const Field& eta10 = H(offset+1,offset);
    const Field& eta11 = H(offset+1,offset+1);
    // NOTE:
    // eta10 should be real, but it is possibly negative, and so we will
    // interpret it as a complex number so that the square-root is well-defined
    EL_DEBUG_ONLY(
      if( ImagPart(eta10) != zero )
          LogicError("Subdiagonal assumed real");
    )

    Field shift = eta11;
    const Field gamma = Sqrt(eta01)*Sqrt(eta10);
    Real sigma = OneAbs(gamma);
    if( sigma != zero )
    {
        const Field xi = (eta00-shift)/Real(2);
        const Real xiAbs = OneAbs(xi);
        sigma = Max( sigma, xiAbs );
        const Field xiSigma = xi/sigma;
        const Field gammaSigma = gamma/sigma;
        Field zeta = sigma*Sqrt(xiSigma*xiSigma+gammaSigma*gammaSigma);
        if( xiAbs > zero )
        {
            const Field xiUnit = xi/xiAbs;
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

namespace single_shift {

template<typename Real>
std::tuple<Int,Complex<Real>,Complex<Real>>
ChooseStart
( const Matrix<Complex<Real>>& H,
        Complex<Real> shift )
{
    EL_DEBUG_CSE
    typedef Complex<Real> Field;
    const Real ulp = limits::Precision<Real>();

    const Int n = H.Height();
    for( Int k=n-3; k>=0; --k )
    {
        const Real eta10 = RealPart( H(k+1,k) );
        const Field& eta11 = H(k+1,k+1);
        const Field& eta22 = H(k+2,k+2);
        Real eta21 = RealPart( H(k+2,k+1) );

        Field eta11Shift = eta11 - shift;
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
    const Field& eta11 = H(0,0);
    Real eta21 = RealPart( H(1,0) );

    const Field eta11Shift = eta11 - shift;
    const Real sigma = OneAbs(eta11Shift) + Abs(eta21);
    return std::make_tuple(0,eta11Shift/sigma,Complex<Real>(eta21/sigma));
}

template<typename Real>
void Sweep
( Matrix<Complex<Real>>& H,
  Complex<Real> shift,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> Field;
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );

    const Int transformBeg = ( ctrl.fullTriangle ? 0 : winBeg );
    const Int transformEnd = ( ctrl.fullTriangle ? n : winEnd );

    // TODO(poulson): Assert that H(k+1,k-1) is real for all k

    auto subInd = IR(winBeg,winEnd);
    auto qrTuple = ChooseStart( H(subInd,subInd), shift );
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
        Field tau0 = lapack::Reflector( 2, nu0, &nu1, 1 );
        if( k > shiftStart )
        {
            H(k,  k-1) = nu0;
            H(k+1,k-1) = 0;
        }
        // The formulas within lapack::Reflector trivially imply that
        // tau0*Conj(nu1) will be real if nu1 was real on entry to
        // lapack::Reflector (an equivalent claim is made within ZLAHQR)
        Real tau1 = RealPart(tau0*Conj(nu1));

        // Apply the Householder reflector from the left
        for( Int j=k; j<transformEnd; ++j )
        {
            Field innerProd = tau0*H(k,j) + tau1*H(k+1,j);
            H(k,  j) -= innerProd;
            H(k+1,j) -= innerProd*nu1;
        }

        // Apply the Householder reflector from the right
        const Int rightApplyEnd = Min(k+3,winEnd);
        for( Int j=transformBeg; j<rightApplyEnd; ++j )
        {
            Field innerProd = Conj(tau0)*H(j,k) + tau1*H(j,k+1);
            H(j,k  ) -= innerProd;
            H(j,k+1) -= innerProd*Conj(nu1);
        }

        if( ctrl.wantSchurVecs )
        {
            // Accumulate the Schur vectors
            for( Int j=0; j<nZ; ++j )
            {
                Field innerProd = Conj(tau0)*Z(j,k) + tau1*Z(j,k+1);
                Z(j,k  ) -= innerProd;
                Z(j,k+1) -= innerProd*Conj(nu1);
            }
        }

        if( k == shiftStart && shiftStart > winBeg )
        {
            // Make H[shiftStart,shiftStart-1] real by scaling by phase
            // TODO(poulson): Investigate more carefully
            Field subdiagVal = Real(1) - Conj(tau0);
            Field phase = subdiagVal / Abs(subdiagVal);
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
                    if( ctrl.wantSchurVecs )
                    {
                        blas::Scal( nZ, Conj(phase), &Z(0,j), 1 );
                    }
                }
            }
        }
    }
    // Make H(winEnd-1,winEnd-2) real by scaling by phase
    Field subdiagVal = H( winEnd-1, winEnd-2 );
    if( ImagPart(subdiagVal) != Real(0) )
    {
        Real subdiagAbs = Abs(subdiagVal);
        H( winEnd-1, winEnd-2 ) = subdiagAbs;
        Field phase = subdiagVal / subdiagAbs;
        if( winEnd < transformEnd )
        {
            blas::Scal
            ( transformEnd-winEnd, Conj(phase),
              &H(winEnd-1,winEnd), H.LDim() );
        }
        blas::Scal
        ( (winEnd-1)-transformBeg, phase, &H(transformBeg,winEnd-1), 1 );
        if( ctrl.wantSchurVecs )
        {
            blas::Scal( nZ, phase, &Z(0,winEnd-1), 1 );
        }
    }
}

// Unfortunately, it seems to be the case that it is noticeably faster
// for this routine to manually inline the data access than to use the
// (presumably inlined) Matrix::operator()(int,int) calls.
template<typename Real>
void SweepOpt
( Matrix<Complex<Real>>& H,
  Complex<Real> shift,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> Field;
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    Field* HBuf = H.Buffer();
    Field* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    const Int transformBeg = ( ctrl.fullTriangle ? 0 : winBeg );
    const Int transformEnd = ( ctrl.fullTriangle ? n : winEnd );

    auto subInd = IR(winBeg,winEnd);
    auto qrTuple = ChooseStart( H(subInd,subInd), shift );
    const Int shiftStart = winBeg + std::get<0>(qrTuple);
    Complex<Real> nu0 = std::get<1>(qrTuple);
    Complex<Real> nu1 = std::get<2>(qrTuple);

    for( Int k=shiftStart; k<winEnd-1; ++k )
    {
        if( k > shiftStart )
        {
            nu0 = HBuf[k+(k-1)*HLDim];
            nu1 = RealPart( HBuf[(k+1)+(k-1)*HLDim] );
        }
        // TODO(poulson): Assert nu1 is real
        Field tau0 = lapack::Reflector( 2, nu0, &nu1, 1 );
        if( k > shiftStart )
        {
            HBuf[ k   +(k-1)*HLDim] = nu0;
            HBuf[(k+1)+(k-1)*HLDim] = 0;
        }
        // The formulas within lapack::Reflector trivially imply that
        // tau0*Conj(nu1) will be real if nu1 was real on entry to
        // lapack::Reflector (an equivalent claim is made within ZLAHQR)
        Real tau1 = RealPart(tau0*Conj(nu1));

        // Apply the Householder reflector from the left
        for( Int j=k; j<transformEnd; ++j )
        {
            Field innerProd = tau0*HBuf[ k   +j*HLDim] +
                          tau1*HBuf[(k+1)+j*HLDim];
            HBuf[ k   +j*HLDim] -= innerProd;
            HBuf[(k+1)+j*HLDim] -= innerProd*nu1;
        }

        // Apply the Householder reflector from the right
        const Int rightApplyEnd = Min(k+3,winEnd);
        for( Int j=transformBeg; j<rightApplyEnd; ++j )
        {
            Field innerProd = Conj(tau0)*HBuf[j+ k   *HLDim] +
                               tau1 *HBuf[j+(k+1)*HLDim];
            HBuf[j+ k   *HLDim] -= innerProd;
            HBuf[j+(k+1)*HLDim] -= innerProd*Conj(nu1);
        }

        if( ctrl.wantSchurVecs )
        {
            // Accumulate the Schur vectors
            for( Int j=0; j<nZ; ++j )
            {
                Field innerProd =
                  Conj(tau0)*ZBuf[j+ k   *ZLDim] +
                       tau1 *ZBuf[j+(k+1)*ZLDim];
                ZBuf[j+ k   *ZLDim] -= innerProd;
                ZBuf[j+(k+1)*ZLDim] -= innerProd*Conj(nu1);
            }
        }

        if( k == shiftStart && shiftStart > winBeg )
        {
            // Make H[shiftStart,shiftStart-1] real by scaling by phase
            // TODO(poulson): Investigate more carefully
            Field subdiagVal = Real(1) - Conj(tau0);
            Field phase = subdiagVal / Abs(subdiagVal);
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
                    if( ctrl.wantSchurVecs )
                    {
                        blas::Scal( nZ, Conj(phase), &Z(0,j), 1 );
                    }
                }
            }
        }
    }
    // Make H(winEnd-1,winEnd-2) real by scaling by phase
    Field subdiagVal = H( winEnd-1, winEnd-2 );
    if( ImagPart(subdiagVal) != Real(0) )
    {
        Real subdiagAbs = Abs(subdiagVal);
        H( winEnd-1, winEnd-2 ) = subdiagAbs;
        Field phase = subdiagVal / subdiagAbs;
        if( winEnd < transformEnd )
        {
            blas::Scal
            ( transformEnd-winEnd, Conj(phase),
              &H(winEnd-1,winEnd), H.LDim() );
        }
        blas::Scal
        ( (winEnd-1)-transformBeg, phase, &H(transformBeg,winEnd-1), 1 );
        if( ctrl.wantSchurVecs )
        {
            blas::Scal( nZ, phase, &Z(0,winEnd-1), 1 );
        }
    }
}

} // namespace single_shift
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESSQR_SINGLE_SHIFT_SWEEP_HPP
