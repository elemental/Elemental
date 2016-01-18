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
  const DistMatrix<F>& A, 
  const DistMatrix<F,STAR,STAR>& t,
        DistMatrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = AOrig.Height();
    const Real infNormAOrig = InfinityNorm( AOrig );
    const Real frobNormAOrig = FrobeniusNorm( AOrig );
    if( g.Rank() == 0 )
        Output("Testing error...");

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
        hessenberg::ApplyQ( LEFT, uplo, NORMAL, A, t, Q );
        if( print )
            Print( Q, "Q" );
        if( display )
            Display( Q, "Q" );
    }

    // Reverse the accumulated Householder transforms
    hessenberg::ApplyQ( LEFT, uplo, ADJOINT, A, t, AOrig );
    hessenberg::ApplyQ( RIGHT, uplo, NORMAL, A, t, AOrig );
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
    const Real infNormError = InfinityNorm( H );
    const Real frobNormError = FrobeniusNorm( H );

    if( g.Rank() == 0 )
        Output
        ("    ||A||_oo = ",infNormAOrig,"\n",
         "    ||A||_F  = ",frobNormAOrig,"\n",
         "    ||H - Q^H A Q||_oo = ",infNormError,"\n",
         "    ||H - Q^H A Q||_F  = ",frobNormError);
}

template<typename F>
void TestHessenberg
( const Grid& g,
  UpperOrLower uplo,
  Int n,
  bool testCorrectness, 
  bool print,
  bool display )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,STAR,STAR> t(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    Uniform( A, n, n );
    if( testCorrectness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    if( g.Rank() == 0 )
        Output("  Starting reduction to Hessenberg form...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Hessenberg( uplo, A, t );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    // TODO: Flop calculation
    if( g.Rank() == 0 )
        Output("  ",runTime," seconds");
    if( print )
    {
        Print( A, "A after Hessenberg" );
        Print( t, "t after Hessenberg" );
    }
    if( display )
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
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
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

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        TestHessenberg<float>
        ( g, uplo, n, testCorrectness, print, display );
        TestHessenberg<Complex<float>>
        ( g, uplo, n, testCorrectness, print, display );

        TestHessenberg<double>
        ( g, uplo, n, testCorrectness, print, display );
        TestHessenberg<Complex<double>>
        ( g, uplo, n, testCorrectness, print, display );

#ifdef EL_HAVE_QD
        TestHessenberg<DoubleDouble>
        ( g, uplo, n, testCorrectness, print, display );
        TestHessenberg<QuadDouble>
        ( g, uplo, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_QUAD
        TestHessenberg<Quad>
        ( g, uplo, n, testCorrectness, print, display );
        TestHessenberg<Complex<Quad>>
        ( g, uplo, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_MPC
        TestHessenberg<BigFloat>
        ( g, uplo, n, testCorrectness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
