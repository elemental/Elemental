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
( const Matrix<F>& A,
  const Matrix<F>& tP,
  const Matrix<F>& tQ,
        Matrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Int m = AOrig.Height();
    const Int n = AOrig.Width();
    const Real infNormAOrig = InfinityNorm( AOrig );
    const Real frobNormAOrig = FrobeniusNorm( AOrig );
    if( mpi::Rank() == 0 )
        Output("Testing error...");

    // Grab the diagonal and superdiagonal of the bidiagonal matrix
    auto d = GetDiagonal( A, 0 );
    auto e = GetDiagonal( A, (m>=n ? 1 : -1) );

    // Zero B and then fill its bidiagonal
    Matrix<F> B;
    Zeros( B, m, n );
    SetDiagonal( B, d, 0  );
    SetDiagonal( B, e, (m>=n ? 1 : -1) );
    if( print && mpi::Rank() == 0 )
        Print( B, "Bidiagonal" );
    if( display && mpi::Rank() == 0 )
        Display( B, "Bidiagonal" );

    if( print || display )
    {
        Matrix<F> Q, P;
        Identity( Q, m, m );
        Identity( P, n, n );
        bidiag::ApplyQ( LEFT,  NORMAL, A, tQ, Q );
        bidiag::ApplyP( RIGHT, NORMAL, A, tP, P );
        if( print && mpi::Rank() == 0 )
        {
            Print( Q, "Q" );
            Print( P, "P" );
        }
        if( display && mpi::Rank() == 0 )
        {
            Display( Q, "Q" );
            Display( P, "P" );
        }
    }

    // Reverse the accumulated Householder transforms
    bidiag::ApplyQ( LEFT,  ADJOINT, A, tQ, AOrig );
    bidiag::ApplyP( RIGHT, NORMAL,  A, tP, AOrig );
    if( print && mpi::Rank() == 0 )
        Print( AOrig, "Manual bidiagonal" );
    if( display && mpi::Rank() == 0 )
        Display( AOrig, "Manual bidiagonal" );

    // Compare the appropriate portion of AOrig and B
    if( m >= n )
    {
        MakeTrapezoidal( UPPER, AOrig );
        MakeTrapezoidal( LOWER, AOrig, 1 );
    }
    else
    {
        MakeTrapezoidal( LOWER, AOrig );
        MakeTrapezoidal( UPPER, AOrig, -1 );
    }
    B -= AOrig;
    if( print && mpi::Rank() == 0 )
        Print( B, "Error in rotated bidiagonal" );
    if( display && mpi::Rank() == 0 )
        Display( B, "Error in rotated bidiagonal" );
    const Real infNormError = InfinityNorm( B );
    const Real frobNormError = FrobeniusNorm( B );

    if( mpi::Rank() == 0 )
        Output
        ("    ||A||_oo = ",infNormAOrig,"\n",
         "    ||A||_F  = ",frobNormAOrig,"\n",
         "    ||B - Q^H A P||_oo = ",infNormError,"\n",
         "    ||B - Q^H A P||_F  = ",frobNormError);
}

template<typename F>
void TestBidiag
( Int m, Int n, bool testCorrectness, bool print, bool display )
{
    if( mpi::Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    Matrix<F> A, AOrig;
    Matrix<F> tP, tQ;

    Uniform( A, m, n );
    if( testCorrectness )
        AOrig = A;
    if( print && mpi::Rank() == 0 )
        Print( A, "A" );
    if( display && mpi::Rank() == 0 )
        Display( A, "A" );

    if( mpi::Rank() == 0 )
        Output("  Starting bidiagonalization...");
    mpi::Barrier( mpi::COMM_WORLD );
    const double startTime = mpi::Time();
    Bidiag( A, tP, tQ );
    mpi::Barrier( mpi::COMM_WORLD );
    const double runTime = mpi::Time() - startTime;
    // TODO: Flop calculation
    if( mpi::Rank() == 0 )
        Output("  Time = ",runTime," seconds");
    if( print && mpi::Rank() == 0 )
    {
        Print( A, "A after Bidiag" );
        Print( tP, "tP after Bidiag" );
        Print( tQ, "tQ after Bidiag" );
    }
    if( display && mpi::Rank() == 0 )
    {
        Display( A, "A after Bidiag" );
        Display( tP, "tP after Bidiag" );
        Display( tQ, "tQ after Bidiag" );
    }
    if( testCorrectness )
        TestCorrectness( A, tP, tQ, AOrig, print, display );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
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

        SetBlocksize( nb );
        ComplainIfDebug();

        TestBidiag<float>( m, n, testCorrectness, print, display );
        TestBidiag<Complex<float>>( m, n, testCorrectness, print, display );

        TestBidiag<double>( m, n, testCorrectness, print, display );
        TestBidiag<Complex<double>>( m, n, testCorrectness, print, display );

#ifdef EL_HAVE_QD
        TestBidiag<DoubleDouble>( m, n, testCorrectness, print, display );
        TestBidiag<QuadDouble>( m, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_QUAD
        TestBidiag<Quad>( m, n, testCorrectness, print, display );
        TestBidiag<Complex<Quad>>( m, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_MPC
        TestBidiag<BigFloat>( m, n, testCorrectness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
