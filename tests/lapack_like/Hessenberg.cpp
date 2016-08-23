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
void TestCorrectness
( UpperOrLower uplo, 
  const Matrix<F>& A, 
  const Matrix<F>& phase,
        Matrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Int n = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormAOrig = OneNorm( AOrig );
    Output("Testing error...");
    PushIndent();

    // Set H to the appropriate Hessenberg portion of A
    Matrix<F> H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "Hessenberg" );
    if( display )
        Display( H, "Hessenberg" );

    if( print || display )
    {
        Matrix<F> Q;
        Identity( Q, n, n );
        hessenberg::ApplyQ( LEFT, uplo, NORMAL, A, phase, Q );
        if( print )
            Print( Q, "Q" );
        if( display )
            Display( Q, "Q" );
    }

    // Reverse the accumulated Householder transforms
    hessenberg::ApplyQ( LEFT, uplo, ADJOINT, A, phase, AOrig );
    hessenberg::ApplyQ( RIGHT, uplo, NORMAL, A, phase, AOrig );
    if( print )
        Print( AOrig, "Manual Hessenberg" );
    if( display )
        Display( AOrig, "Manual Hessenberg" );

    // Compare the appropriate portion of AOrig and B
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, AOrig, 1 );
    else
        MakeTrapezoidal( UPPER, AOrig, -1 );
    H -= AOrig;
    if( print )
        Print( H, "Error in rotated Hessenberg" );
    if( display )
        Display( H, "Error in rotated Hessenberg" );
    const Real infError = InfinityNorm( H );
    const Real relError = infError / (n*eps*oneNormAOrig);

    Output("||H - Q^H A Q||_oo / (eps n || A ||_1) = ",relError);

    // TODO: Use a more refined failure condition
    if( relError > Real(1) )
        LogicError("Unacceptably large relative error");

    PopIndent();
}

template<typename F> 
void TestCorrectness
( UpperOrLower uplo, 
  const DistMatrix<F>& A, 
  const DistMatrix<F,STAR,STAR>& phase,
        DistMatrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = AOrig.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormAOrig = OneNorm( AOrig );
    OutputFromRoot(g.Comm(),"Testing error...");
    PushIndent();

    // Set H to the appropriate Hessenberg portion of A
    DistMatrix<F> H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "Hessenberg" );
    if( display )
        Display( H, "Hessenberg" );

    if( print || display )
    {
        DistMatrix<F> Q(g);
        Identity( Q, n, n );
        hessenberg::ApplyQ( LEFT, uplo, NORMAL, A, phase, Q );
        if( print )
            Print( Q, "Q" );
        if( display )
            Display( Q, "Q" );
    }

    // Reverse the accumulated Householder transforms
    hessenberg::ApplyQ( LEFT, uplo, ADJOINT, A, phase, AOrig );
    hessenberg::ApplyQ( RIGHT, uplo, NORMAL, A, phase, AOrig );
    if( print )
        Print( AOrig, "Manual Hessenberg" );
    if( display )
        Display( AOrig, "Manual Hessenberg" );

    // Compare the appropriate portion of AOrig and B
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, AOrig, 1 );
    else
        MakeTrapezoidal( UPPER, AOrig, -1 );
    H -= AOrig;
    if( print )
        Print( H, "Error in rotated Hessenberg" );
    if( display )
        Display( H, "Error in rotated Hessenberg" );
    const Real infError = InfinityNorm( H );
    const Real relError = infError / (n*eps*oneNormAOrig);

    OutputFromRoot
    (g.Comm(),"||H - Q^H A Q||_oo / (eps n || A ||_1) = ",relError);

    // TODO: Use a more refined failure condition
    if( relError > Real(1) )
        LogicError("Unacceptably large relative error");

    PopIndent();
}

template<typename F>
void TestHessenberg
( UpperOrLower uplo,
  Int n,
  bool correctness, 
  bool print,
  bool display )
{
    Matrix<F> A, AOrig;
    Matrix<F> phase;
    Output("Testing with ",TypeName<F>());
    PushIndent();

    Uniform( A, n, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    Output("Starting reduction to Hessenberg form...");
    Timer timer;
    timer.Start();
    Hessenberg( uplo, A, phase );
    const double runTime = timer.Stop();
    // TODO: Flop calculation
    Output(runTime," seconds");
    if( print )
    {
        Print( A, "A after Hessenberg" );
        Print( phase, "phase after Hessenberg" );
    }
    if( display )
    {
        Display( A, "A after Hessenberg" );
        Display( phase, "phase after Hessenberg" );
    }
    if( correctness )
        TestCorrectness( uplo, A, phase, AOrig, print, display );
    PopIndent();
}

template<typename F>
void TestHessenberg
( const Grid& g,
  UpperOrLower uplo,
  Int n,
  bool correctness, 
  bool print,
  bool display )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,STAR,STAR> phase(g);
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    Uniform( A, n, n );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    OutputFromRoot(g.Comm(),"Starting reduction to Hessenberg form...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    Hessenberg( uplo, A, phase );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    // TODO: Flop calculation
    OutputFromRoot(g.Comm(),runTime," seconds");
    if( print )
    {
        Print( A, "A after Hessenberg" );
        Print( phase, "phase after Hessenberg" );
    }
    if( display )
    {
        Display( A, "A after Hessenberg" );
        Display( phase, "phase after Hessenberg" );
    }
    if( correctness )
        TestCorrectness( uplo, A, phase, AOrig, print, display );
    PopIndent();
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool correctness = 
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
#endif

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( sequential && mpi::Rank() == 0 )
        {
            TestHessenberg<float>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<float>>
            ( uplo, n, correctness, print, display );

            TestHessenberg<double>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<double>>
            ( uplo, n, correctness, print, display );

#ifdef EL_HAVE_QD
            TestHessenberg<DoubleDouble>
            ( uplo, n, correctness, print, display );
            TestHessenberg<QuadDouble>
            ( uplo, n, correctness, print, display );

            TestHessenberg<Complex<DoubleDouble>>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<QuadDouble>>
            ( uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_QUAD
            TestHessenberg<Quad>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<Quad>>
            ( uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_MPC
            TestHessenberg<BigFloat>
            ( uplo, n, correctness, print, display );
            TestHessenberg<Complex<BigFloat>>
            ( uplo, n, correctness, print, display );
#endif
        }

        TestHessenberg<float>
        ( g, uplo, n, correctness, print, display );
        TestHessenberg<Complex<float>>
        ( g, uplo, n, correctness, print, display );

        TestHessenberg<double>
        ( g, uplo, n, correctness, print, display );
        TestHessenberg<Complex<double>>
        ( g, uplo, n, correctness, print, display );

#ifdef EL_HAVE_QD
        TestHessenberg<DoubleDouble>
        ( g, uplo, n, correctness, print, display );
        TestHessenberg<QuadDouble>
        ( g, uplo, n, correctness, print, display );

        TestHessenberg<Complex<DoubleDouble>>
        ( g, uplo, n, correctness, print, display );
        TestHessenberg<Complex<QuadDouble>>
        ( g, uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_QUAD
        TestHessenberg<Quad>
        ( g, uplo, n, correctness, print, display );
        TestHessenberg<Complex<Quad>>
        ( g, uplo, n, correctness, print, display );
#endif

#ifdef EL_HAVE_MPC
        TestHessenberg<BigFloat>
        ( g, uplo, n, correctness, print, display );
        TestHessenberg<Complex<BigFloat>>
        ( g, uplo, n, correctness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
