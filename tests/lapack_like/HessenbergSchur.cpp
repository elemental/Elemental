/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.  
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F>
void TestAhuesTisseur( const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = 3;
    const Real eps = limits::Epsilon<Real>();
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
    const Real HFrob = FrobeniusNorm( H );
    const Real relErr = errFrob / (eps*n*HFrob);
    Output("|| H ||_F = ",HFrob);
    Output("|| H Z - Z T ||_F / (eps n || H ||_F) = ",relErr);
    if( print )
        Print( R );
    // TODO(poulson): A more refined failure condition
    if( relErr > Real(100) )
        LogicError("Relative error was unacceptably large");
    else
        Output("Passed test");
    Output("");
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void TestLAPACKHelper
( const Matrix<F>& H, 
  const HessenbergSchurCtrl& ctrl,
  bool testSweep,
  bool print )
{
    DEBUG_CSE
    const Int n = H.Height();
    Matrix<F> T, Z;
    Matrix<Complex<Base<F>>> w;

    T = H;
    w.Resize( n, 1 );
    Identity( Z, n, n ); 
     
    Timer timer;
    timer.Start();
    const bool multiplyZ = true;
    const bool useAED = ( ctrl.alg == HESSENBERG_SCHUR_AED );
    lapack::HessenbergSchur
    ( n, T.Buffer(), T.LDim(), w.Buffer(), Z.Buffer(), Z.LDim(),
      ctrl.fullTriangle, multiplyZ, useAED );
    Output("LAPACK HessenbergSchur: ",timer.Stop()," seconds");
}

template<typename F>
void TestRandomHelper
( const Matrix<F>& H, 
  const HessenbergSchurCtrl& ctrl,
  bool testSweep,
  bool print )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = H.Height();
    const Real eps = limits::Epsilon<Real>();

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

    if( testSweep )
    {
        Matrix<Complex<Real>> wCopy;
        T = H;
        auto ctrlEig( ctrl );
        ctrlEig.fullTriangle = false;
        ctrlEig.wantSchurVecs = false;
        timer.Start();
        auto infoEig = HessenbergSchur( T, wCopy, ctrlEig );
        Output("HessenbergSchur (eigenvalues only): ",timer.Stop()," seconds");
        Output
        ("Convergence achieved after ",infoEig.numIterations," iterations");

        // Test *without* sorting
        Output("Testing perfect-shift eigensolver preconditioner");
        T = H;
        Identity( Z, n, n );
        auto ctrlSweep( ctrl );
        ctrlSweep.wantSchurVecs = true;
        timer.Start();
        w = wCopy;
        hess_schur::Sweep( T, w, Z, ctrlSweep );
        Output("hess_schur::Sweep: ",timer.Stop()," seconds");
        if( print )
        {
            Print( T, "T after sweep" );
            Print( Z, "Z after sweep" );
        }
        timer.Start();
        auto infoFinish = HessenbergSchur( T, w, Z, ctrl );
        Output("HessenbergSchur (after sweep): ",timer.Stop()," seconds");
        Output
        ("Convergence achieved after ",infoFinish.numIterations," iterations");
    }

    Matrix<F> R;
    Gemm( NORMAL, NORMAL, F(1), Z, T, R );
    Gemm( NORMAL, NORMAL, F(1), H, Z, F(-1), R );
    const Real errFrob = FrobeniusNorm( R ); 
    const Real HFrob = FrobeniusNorm( H );
    const Real relErr = errFrob / (eps*n*HFrob);
    Output("|| H ||_F = ",HFrob);
    Output("|| H Z - Z T ||_F / (eps n || H ||_F) = ",relErr);
    if( print )
        Print( R );
    // TODO(poulson): A more refined failure condition
    if( relErr > Real(100) )
        LogicError("Relative error was unacceptably large");
    else
        Output("Passed test");
    Output("");
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void TestRandom
( Int n, const HessenbergSchurCtrl& ctrl, bool testSweep, bool print )
{
    DEBUG_CSE
    Output("Testing uniform Hessenberg with ",TypeName<F>());

    Matrix<F> H;
    Uniform( H, n, n );
    MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "H" );

    TestLAPACKHelper( H, ctrl, testSweep, print );
    TestRandomHelper( H, ctrl, testSweep, print );
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
void TestRandom
( Int n, const HessenbergSchurCtrl& ctrl, bool testSweep, bool print )
{
    DEBUG_CSE
    Output("Testing uniform Hessenberg with ",TypeName<F>());

    Matrix<F> H;
    Uniform( H, n, n );
    MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "H" );

    TestRandomHelper( H, ctrl, testSweep, print );
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","random matrix size",60);
        const Int algInt = Input("--alg","AED: 0, MultiBulge: 1, Simple: 2",0);
        const Int minMultiBulgeSize =
          Input
          ("--minMultiBulgeSize",
           "minimum size for using a multi-bulge algorithm",75);
        const bool accumulate =
          Input("--accumulate","accumulate reflections?",true);
        const bool sortShifts =
          Input("--sortShifts","sort shifts for AED?",true);
        const bool testSweep =
          Input("--testSweep","test pure-shift sweep?",false);
        const bool progress = Input("--progress","print progress?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        HessenbergSchurCtrl ctrl;
        ctrl.alg = static_cast<HessenbergSchurAlg>(algInt);
        ctrl.minMultiBulgeSize = minMultiBulgeSize;
        ctrl.accumulateReflections = accumulate;
        ctrl.sortShifts = sortShifts;
        ctrl.progress = progress;

        TestAhuesTisseur<float>( ctrl, print );
        TestAhuesTisseur<Complex<float>>( ctrl, print );
        TestAhuesTisseur<double>( ctrl, print );
        TestAhuesTisseur<Complex<double>>( ctrl, print );
#ifdef EL_HAVE_QUAD
        TestAhuesTisseur<Quad>( ctrl, print );
        TestAhuesTisseur<Complex<Quad>>( ctrl, print );
#endif
#ifdef EL_HAVE_QD
        TestAhuesTisseur<DoubleDouble>( ctrl, print );
        TestAhuesTisseur<Complex<DoubleDouble>>( ctrl, print );
        TestAhuesTisseur<QuadDouble>( ctrl, print );
        TestAhuesTisseur<Complex<QuadDouble>>( ctrl, print );
#endif
#ifdef EL_HAVE_MPC
        TestAhuesTisseur<BigFloat>( ctrl, print );
        TestAhuesTisseur<Complex<BigFloat>>( ctrl, print );
#endif

        TestRandom<float>( n, ctrl, testSweep, print );
        TestRandom<Complex<float>>( n, ctrl, testSweep, print );
        TestRandom<double>( n, ctrl, testSweep, print );
        TestRandom<Complex<double>>( n, ctrl, testSweep, print );
#ifdef EL_HAVE_QUAD
        TestRandom<Quad>( n, ctrl, testSweep, print );
        TestRandom<Complex<Quad>>( n, ctrl, testSweep, print );
#endif
#ifdef EL_HAVE_QD
        TestRandom<DoubleDouble>( n, ctrl, testSweep, print );
        TestRandom<Complex<DoubleDouble>>( n, ctrl, testSweep, print );
        TestRandom<QuadDouble>( n, ctrl, testSweep, print );
        TestRandom<Complex<QuadDouble>>( n, ctrl, testSweep, print );
#endif
#ifdef EL_HAVE_MPC
        TestRandom<BigFloat>( n, ctrl, testSweep, print );
        TestRandom<Complex<BigFloat>>( n, ctrl, testSweep, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
