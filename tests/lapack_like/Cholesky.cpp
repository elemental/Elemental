/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F>
void TestCorrectness
( bool pivot,
  UpperOrLower uplo,
  const DistMatrix<F>& A,
  const DistPermutation& p,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    // Test correctness by multiplying a random set of vectors by A, then
    // using the Cholesky factorization to solve.
    DistMatrix<F> X(g), Y(g);
    Uniform( X, m, 100 );
    Zeros( Y, m, 100 );
    Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
    const Real maxNormL = HermitianMaxNorm( uplo, A );
    const Real maxNormA = HermitianMaxNorm( uplo, AOrig );
    const Real frobNormA = HermitianFrobeniusNorm( uplo, AOrig );
    const Real frobNormY = FrobeniusNorm( Y );

    if( pivot )
        cholesky::SolveAfter( uplo, NORMAL, A, p, Y );
    else
        cholesky::SolveAfter( uplo, NORMAL, A, Y );
    X -= Y;
    const Real frobNormE = FrobeniusNorm( X );

    if( g.Rank() == 0 )
        Output
        ("||L||_max            = ",maxNormL,"\n",
         "||A||_max            = ",maxNormA,"\n",
         "||A||_F              = ",frobNormA,"\n",
         "||Y||_F              = ",frobNormY,"\n",
         "||X - inv(A) X||_F  = ",frobNormE);
}

template<typename F> 
void TestCholesky
( bool testCorrectness,
  bool pivot,
  bool print,
  bool printDiag,
  UpperOrLower uplo,
  Int m,
  const Grid& g,
  bool scalapack )
{
    DistMatrix<F> A(g), AOrig(g);
    DistPermutation p(g);

    HermitianUniformSpectrum( A, m, 1e-9, 10 );
    if( testCorrectness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    if( g.Rank() == 0 )
    {
        if( scalapack && !pivot )
            Output("  ScaLAPACK Cholesky (including round-trip conversion)...");
        else
            Output("  Elemental Cholesky...");
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( pivot )
        Cholesky( uplo, A, p );
    else
        Cholesky( uplo, A, scalapack );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    { 
        Print( A, "A after factorization" );
        if( pivot )
        {
            DistMatrix<Int,VC,STAR> P(g);
            p.ExplicitMatrix( P );
            Print( P, "P" );
        }
    }
    if( printDiag )
        Print( GetRealPartOfDiagonal(A), "diag(A)" );
    if( testCorrectness )
        TestCorrectness( pivot, uplo, A, p, AOrig );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--m","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool pivot = Input("--pivot","use pivoting?",false);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool printDiag = Input("--printDiag","print diag of fact?",false);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
#else
        const bool scalapack = false;
#endif
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalTrrkBlocksize<double>( nbLocal );
        SetLocalTrrkBlocksize<Complex<double>>( nbLocal );
        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Cholesky",uploChar);

        if( commRank == 0 )
            Output("Testing with doubles:");
        if( scalapack )
            TestCholesky<double>
            ( testCorrectness, pivot, print, printDiag, uplo, m, g, true );
        TestCholesky<double>
        ( testCorrectness, pivot, print, printDiag, uplo, m, g, false );

        if( commRank == 0 )
            Output("Testing with double-precision complex:");
        if( scalapack )
            TestCholesky<Complex<double>>
            ( testCorrectness, pivot, print, printDiag, uplo, m, g, true );
        TestCholesky<Complex<double>>
        ( testCorrectness, pivot, print, printDiag, uplo, m, g, false );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
