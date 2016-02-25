/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F> 
void TestCorrectness
( UpperOrLower uplo,
  const Matrix<F>& A,
  const Matrix<F>& t,
        Matrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Int n = AOrig.Height();
    const Real infNormAOrig = InfinityNorm( AOrig );
    const Real frobNormAOrig = FrobeniusNorm( AOrig );
    if( mpi::Rank() == 0 )
        Output("Testing error...");

    // Set H to the appropriate Hessenberg portion of A
    Matrix<F> H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );
    if( print && mpi::Rank() == 0 )
        Print( H, "Hessenberg" );
    if( display && mpi::Rank() == 0 )
        Display( H, "Bidiagonal" );

    if( (print || display) && mpi::Rank() == 0 )
    {
        Matrix<F> Q;
        Identity( Q, n, n );
        hessenberg::ApplyQ( LEFT, uplo, NORMAL, A, t, Q );
        if( print )
            Print( Q, "Q" );
        if( display )
            Display( Q, "Q" );
    }

    // Reverse the accumulated Householder transforms
    hessenberg::ApplyQ( LEFT, uplo, ADJOINT, A, t, AOrig );
    hessenberg::ApplyQ( RIGHT, uplo, NORMAL, A, t, AOrig );
    if( print && mpi::Rank() == 0 )
        Print( AOrig, "Manual Hessenberg" );
    if( display && mpi::Rank() == 0 )
        Display( AOrig, "Manual Hessenberg" );

    // Compare the appropriate portion of AOrig and B
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, AOrig, 1 );
    else
        MakeTrapezoidal( UPPER, AOrig, -1 );
    H -= AOrig;
    if( print && mpi::Rank() == 0 )
        Print( H, "Error in rotated Hessenberg" );
    if( display && mpi::Rank() == 0 )
        Display( H, "Error in rotated Hessenberg" );
    const Real infNormError = InfinityNorm( H );
    const Real frobNormError = FrobeniusNorm( H );

    if( mpi::Rank() == 0 )
        Output
        ("    ||A||_oo = ",infNormAOrig,"\n",
         "    ||A||_F  = ",frobNormAOrig,"\n",
         "    ||H - Q^H A Q||_oo = ",infNormError,"\n",
         "    ||H - Q^H A Q||_F  = ",frobNormError);
}

template<typename F>
void TestHessenberg
( UpperOrLower uplo,
  Int n,
  bool testCorrectness,
  bool print,
  bool display )
{
    if( mpi::Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    Matrix<F> A, AOrig;
    Matrix<F> t;

    Uniform( A, n, n );
    if( testCorrectness )
        AOrig = A;
    if( print && mpi::Rank() == 0 )
        Print( A, "A" );
    if( display && mpi::Rank() == 0 )
        Display( A, "A" );

    if( mpi::Rank() == 0 )
        Output("  Starting reduction to Hessenberg form...");
    mpi::Barrier( mpi::COMM_WORLD );
    const double startTime = mpi::Time();
    Hessenberg( uplo, A, t );
    mpi::Barrier( mpi::COMM_WORLD );
    const double runTime = mpi::Time() - startTime;
    // TODO: Flop calculation
    if( mpi::Rank() == 0 )
        Output("  Time = ",runTime," seconds");
    if( print && mpi::Rank() == 0 )
    {
        Print( A, "A after Hessenberg" );
        Print( t, "t after Hessenberg" );
    }
    if( display && mpi::Rank() == 0 )
    {
        Display( A, "A after Hessenberg" );
        Display( t, "t after Hessenberg" );
    }
    if( testCorrectness )
        TestCorrectness( uplo, A, t, AOrig, print, display );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpc::SetPrecision( prec );
#endif

        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        TestHessenberg<float>
        ( uplo, n, testCorrectness, print, display );
        TestHessenberg<Complex<float>>
        ( uplo, n, testCorrectness, print, display );

        TestHessenberg<double>
        ( uplo, n, testCorrectness, print, display );
        TestHessenberg<Complex<double>>
        ( uplo, n, testCorrectness, print, display );

#ifdef EL_HAVE_QD
        TestHessenberg<DoubleDouble>
        ( uplo, n, testCorrectness, print, display );
        TestHessenberg<QuadDouble>
        ( uplo, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_QUAD
        TestHessenberg<Quad>
        ( uplo, n, testCorrectness, print, display );
        TestHessenberg<Complex<Quad>>
        ( uplo, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_MPC
        TestHessenberg<BigFloat>
        ( uplo, n, testCorrectness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
