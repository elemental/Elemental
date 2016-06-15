/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

#include "./HessenbergEig/WindowedSingle.hpp"
#include "./HessenbergEig/SmallBulgeSweep.hpp"
#include "./HessenbergEig/AED.hpp"

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
    H(0,0) = F(1.);
    H(0,1) = F(1.1e5);
    H(0,2) = F(0.);
    H(1,0) = F(1.1e-8);
    H(1,1) = F(1.+1.e-2);
    H(1,2) = F(1.1e5);
    H(2,0) = F(0.);
    H(2,1) = F(1.1e-8);
    H(2,2) = F(1.+2.*1.e-2);
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
    H(0,0) = F(1.);
    H(0,1) = F(1.1e5);
    H(0,2) = F(0.);
    H(1,0) = F(1.1e-8);
    H(1,1) = F(1.+1.e-2);
    H(1,2) = F(1.1e5);
    H(2,0) = F(0.);
    H(2,1) = F(1.1e-8);
    H(2,2) = F(1.+2.*1.e-2);
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
        const Int n = Input("--n","random matrix size",60);
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
#ifdef EL_HAVE_QD
        TestAhuesTisseur<DoubleDouble>( print );
        TestAhuesTisseur<QuadDouble>( print );
#endif
#ifdef EL_HAVE_MPC
        TestAhuesTisseur<BigFloat>( print );
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
#ifdef EL_HAVE_QD
        TestRandom<DoubleDouble>( n, print );
        TestRandom<QuadDouble>( n, print );
#endif
#ifdef EL_HAVE_MPC
        TestRandom<BigFloat>( n, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
