/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Real>
void TestAhuesTisseur( const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
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
    auto info = HessenbergSchur( T, w, Z, ctrl );
    Output("HessenbergSchur: ",timer.Stop()," seconds");
    Output("Convergence achieved after ",info.numIterations," iterations");
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
void TestAhuesTisseurQuasi( const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
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

    Matrix<F> T, Z;
    Matrix<Complex<Real>> w;
    Identity( Z, n, n );
    T = H;
    Timer timer;
    timer.Start();
    auto info = HessenbergSchur( T, w, Z, ctrl );
    Output("HessenbergSchur: ",timer.Stop()," seconds");
    Output("Convergence achieved after ",info.numIterations," iterations");
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

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void TestRandom( Int n, const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
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
    Timer timer;

    T = H;
    w.Resize( n, 1 );
    Identity( Z, n, n ); 
    timer.Start();
    bool multiplyZ = true;
    lapack::HessenbergSchur
    ( n, T.Buffer(), T.LDim(), w.Buffer(), Z.Buffer(), Z.LDim(),
      ctrl.fullTriangle, multiplyZ, ctrl.useAED );
    Output("LAPACK HessenbergSchur: ",timer.Stop()," seconds");

    T = H;
    Identity( Z, n, n );
    timer.Start();
    auto info = HessenbergSchur( T, w, Z, ctrl );
    Output("HessenbergSchur: ",timer.Stop()," seconds");
    Output("Convergence achieved after ",info.numIterations," iterations");
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

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void TestRandom( Int n, const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
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
    Timer timer;

    T = H;
    Identity( Z, n, n );
    timer.Start();
    auto info = HessenbergSchur( T, w, Z, ctrl );
    Output("HessenbergSchur: ",timer.Stop()," seconds");
    Output("Convergence achieved after ",info.numIterations," iterations");
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

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void TestRandomQuasi( Int n, const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
    typedef Real F;
    Output("Testing uniform Hessenberg with ",TypeName<F>());

    Matrix<F> H;
    Uniform( H, n, n );
    MakeTrapezoidal( UPPER, H, -1 );
    const Real HFrob = FrobeniusNorm( H );
    Output("|| H ||_F = ",HFrob);
    if( print )
        Print( H, "H" );

    Matrix<F> T, Z;
    Matrix<Complex<Real>> w;
    Timer timer;

    T = H;
    w.Resize( n, 1 );
    Identity( Z, n, n ); 
    timer.Start();
    bool multiplyZ = true;
    lapack::HessenbergSchur
    ( n, T.Buffer(), T.LDim(), w.Buffer(), Z.Buffer(), Z.LDim(),
      ctrl.fullTriangle, multiplyZ, ctrl.useAED );
    Output("LAPACK HessenbergSchur: ",timer.Stop()," seconds");

    T = H;
    Identity( Z, n, n );
    timer.Start();
    auto info = HessenbergSchur( T, w, Z, ctrl );
    Output("HessenbergSchur: ",timer.Stop()," seconds");
    Output("Convergence achieved after ",info.numIterations," iterations");
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

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void TestRandomQuasi( Int n, const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
    typedef Real F;
    Output("Testing uniform Hessenberg with ",TypeName<F>());

    Matrix<F> H;
    Uniform( H, n, n );
    MakeTrapezoidal( UPPER, H, -1 );
    const Real HFrob = FrobeniusNorm( H );
    Output("|| H ||_F = ",HFrob);
    if( print )
        Print( H, "H" );

    Matrix<F> T, Z;
    Matrix<Complex<Real>> w;
    Timer timer;

    T = H;
    Identity( Z, n, n );
    timer.Start();
    auto info = HessenbergSchur( T, w, Z, ctrl );
    Output("HessenbergSchur: ",timer.Stop()," seconds");
    Output("Convergence achieved after ",info.numIterations," iterations");
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
        const Int n = Input("--n","random matrix size",60);
        const bool useAED = Input("--aed","use Aggressive Early Deflat?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        HessenbergSchurCtrl ctrl;
        ctrl.useAED = useAED;
        ctrl.progress = progress;

        TestAhuesTisseurQuasi<float>( ctrl, print );
        TestAhuesTisseurQuasi<double>( ctrl, print );
#ifdef EL_HAVE_QUAD
        TestAhuesTisseurQuasi<Quad>( ctrl, print );
#endif
#ifdef EL_HAVE_QD
        TestAhuesTisseurQuasi<DoubleDouble>( ctrl, print );
        TestAhuesTisseurQuasi<QuadDouble>( ctrl, print );
#endif
#ifdef EL_HAVE_MPC
        TestAhuesTisseurQuasi<BigFloat>( ctrl, print );
#endif

        TestAhuesTisseur<float>( ctrl, print );
        TestAhuesTisseur<double>( ctrl, print );
#ifdef EL_HAVE_QUAD
        TestAhuesTisseur<Quad>( ctrl, print );
#endif
#ifdef EL_HAVE_QD
        TestAhuesTisseur<DoubleDouble>( ctrl, print );
        TestAhuesTisseur<QuadDouble>( ctrl, print );
#endif
#ifdef EL_HAVE_MPC
        TestAhuesTisseur<BigFloat>( ctrl, print );
#endif

        TestRandomQuasi<float>( n, ctrl, print );
        TestRandomQuasi<double>( n, ctrl, print );
#ifdef EL_HAVE_QUAD
        TestRandomQuasi<Quad>( n, ctrl, print );
#endif
#ifdef EL_HAVE_QD
        TestRandomQuasi<DoubleDouble>( n, ctrl, print );
        TestRandomQuasi<QuadDouble>( n, ctrl, print );
#endif
#ifdef EL_HAVE_MPC
        TestRandomQuasi<BigFloat>( n, ctrl, print );
#endif

        TestRandom<float>( n, ctrl, print );
        TestRandom<double>( n, ctrl, print );
#ifdef EL_HAVE_QUAD
        TestRandom<Quad>( n, ctrl, print );
#endif
#ifdef EL_HAVE_QD
        TestRandom<DoubleDouble>( n, ctrl, print );
        TestRandom<QuadDouble>( n, ctrl, print );
#endif
#ifdef EL_HAVE_MPC
        TestRandom<BigFloat>( n, ctrl, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
