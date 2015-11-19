/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F,Dist UPerm> 
void TestCorrectness
( bool print, 
  const DistMatrix<F>& A,
  const DistMatrix<Int,UPerm,STAR>& p,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = AOrig.Width();

    if( g.Rank() == 0 )
        Output("Testing error...");

    // Generate random right-hand sides
    DistMatrix<F> X(g);
    Uniform( X, n, 100 );
    auto Y( X );
    PermuteRows( Y, p );
    lu::SolveAfter( NORMAL, A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    const Real infNormX = InfinityNorm( X );
    const Real frobNormX = FrobeniusNorm( X );
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real infNormError = InfinityNorm( X );
    const Real frobNormError = FrobeniusNorm( X );
    const Real infNormA = InfinityNorm( AOrig );
    const Real frobNormA = FrobeniusNorm( AOrig );

    if( g.Rank() == 0 )
        Output
        ("||A||_oo            = ",infNormA,"\n",
         "||A||_F             = ",frobNormA,"\n",
         "||X||_oo            = ",infNormX,"\n",
         "||X||_F             = ",frobNormX,"\n",
         "||A A^-1 X - X||_oo = ",infNormError,"\n",
         "||A A^-1 X - X||_F  = ",frobNormError);
}

template<typename F,Dist UPerm> 
void TestLUMod
( bool conjugate,
  Base<F> tau,
  bool testCorrectness,
  bool print,
  Int m,
  const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<Int,UPerm,STAR> rowPiv(g), p(g);

    Uniform( A, m, m );
    if( testCorrectness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    {
        if( g.Rank() == 0 )
            Output("  Starting LU factorization...");
        mpi::Barrier( g.Comm() );
        const double startTime = mpi::Time();
        LU( A, rowPiv );
        mpi::Barrier( g.Comm() );
        const double runTime = mpi::Time() - startTime;
        const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
        const double gFlops =
          ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
        if( g.Rank() == 0 )
            Output("  ",runTime," seconds (",gFlops," GFlop/s)");
    }

    PivotsToPermutation( rowPiv, p );
    if( print )
    {
        Print( A, "A after original factorization" );
        Print( p, "p after original factorization");
        DistMatrix<Int> P(g);
        ExplicitPermutation( p, P );
        Print( P, "P" );
    }

    // Generate random vectors u and v
    DistMatrix<F> u(g), v(g);
    Uniform( u, m, 1 );
    Uniform( v, m, 1 );
    if( testCorrectness )
    {
        if( conjugate )
            Ger( F(1), u, v, AOrig );
        else
            Geru( F(1), u, v, AOrig );
    }

    { 
        if( g.Rank() == 0 )
            Output("  Starting rank-one LU modification...");
        mpi::Barrier( g.Comm() );
        const double startTime = mpi::Time();
        LUMod( A, p, u, v, conjugate, tau );
        mpi::Barrier( g.Comm() );
        const double runTime = mpi::Time() - startTime;
        if( g.Rank() == 0 )
            Output("  ",runTime," seconds");
    }

    if( print )
    {
        Print( A, "A after modification" );
        Print( p, "p after modification");
        DistMatrix<Int> P(g);
        ExplicitPermutation( p, P );
        Print( P, "P" );
    }

    if( testCorrectness )
        TestCorrectness( print, A, p, AOrig );
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
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const double tau = Input("--tau","pivot threshold",0.1);
        const bool conjugate = Input("--conjugate","conjugate v?",true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            Output("Testing with doubles:");
        TestLUMod<double,VC>( conjugate, tau, testCorrectness, print, m, g );

        if( commRank == 0 )
            Output("Testing with double-precision complex:");
        TestLUMod<Complex<double>,VC>
        ( conjugate, tau, testCorrectness, print, m, g );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
