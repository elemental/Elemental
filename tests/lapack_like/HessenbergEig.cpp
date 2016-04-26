/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Use Ahues and Tissuer's (LAPACK Working Note 122, 1997) refinement of
// Wilkinson's criteria for determining if a subdiagonal entry of a Hessenberg
// matrix is negligible
template<typename Real>
Int DetectSmallSubdiagonal( const Matrix<Complex<Real>>& H )
{
    typedef Complex<Real> F;

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
    //DEBUG_ONLY(
      if( n < 2 )
          LogicError("WilkinsonShift requires n >= 2");
    //)
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
    //DEBUG_ONLY(
      if( ImagPart(eta10) != zero )
          LogicError("Subdiagonal assumed real");
    //)

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
DecideShiftStart( const Matrix<Complex<Real>>& H, Complex<Real> shift )
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

    //DEBUG_ONLY(
      if( windowBeg > shiftStart || shiftStart >= windowEnd )
          LogicError
          ("shiftStart was ",shiftStart,", but window was [",windowBeg,",",
           windowEnd,")");
    //)

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
void HessenbergQR
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

            auto qrTuple = DecideShiftStart( H(subInd,subInd), shift );
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
    const bool fullTriangle=true, wantSchurVecs=true;
    Identity( Z, n, n );
    T = H;
    HessenbergQR
    ( T, 0, T.Height(), w, fullTriangle,
      Z, wantSchurVecs, 0, T.Height() );
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
    const bool fullTriangle=true, wantSchurVecs=true;
    Identity( Z, H.Height(), H.Height() );
    T = H;
    HessenbergQR
    ( T, 0, T.Height(), w, fullTriangle,
      Z, wantSchurVecs, 0, T.Height() );
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

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","random matrix size",200);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        TestAhuesTisseur<float>( print );
        TestAhuesTisseur<double>( print );
#ifdef EL_HAVE_QUAD
        TestAhuesTisseur<Quad>( print );
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
