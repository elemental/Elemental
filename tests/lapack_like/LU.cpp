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
( const DistMatrix<F>& AOrig,
  const DistMatrix<F>& A,
  const DistMatrix<Int,UPerm,STAR>& rowPiv,
  const DistMatrix<Int,UPerm,STAR>& colPiv,
  Int pivoting, bool print )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    if( g.Rank() == 0 )
        Output("Testing error...");

    // Generate random right-hand sides
    DistMatrix<F> X(g);
    Uniform( X, m, 100 );
    auto Y( X );
    if( pivoting == 0 )
        lu::SolveAfter( NORMAL, A, Y );
    else if( pivoting == 1 )
        lu::SolveAfter( NORMAL, A, rowPiv, Y );
    else
        lu::SolveAfter( NORMAL, A, rowPiv, colPiv, Y );

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
void TestLU
( Int m, const Grid& g, Int pivoting, 
  bool testCorrectness, bool forceGrowth, bool print )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<Int,UPerm,STAR> rowPiv(g), colPiv(g);

    if( forceGrowth )
        GEPPGrowth( A, m );
    else
        Uniform( A, m, m );

    if( testCorrectness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    if( g.Rank() == 0 )
        Output("  Starting LU factorization...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( pivoting == 0 )
        LU( A );
    else if( pivoting == 1 )
        LU( A, rowPiv );
    else if( pivoting == 2 )
        LU( A, rowPiv, colPiv );

    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivoting >= 1 )
        {
            DistMatrix<Int> P(g);
            DistMatrix<Int,STAR,STAR> p(g);
            Print( rowPiv, "rowPiv after factorization" );
            PivotsToPermutation( rowPiv, p );
            Print( p, "p after factorization");
            ExplicitPermutation( p, P );
            Print( P, "P" );
        }
        if( pivoting == 2 )
        {
            DistMatrix<Int> Q(g);
            DistMatrix<Int,STAR,STAR> q(g);
            Print( colPiv, "colPiv after factorization" );
            PivotsToPermutation( colPiv, q );
            Print( q, "q after factorization");
            ExplicitPermutation( q, Q );
            Print( Q, "Q" );
        }
    }
    if( testCorrectness )
        TestCorrectness( AOrig, A, rowPiv, colPiv, pivoting, print );
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
        const Int pivot = Input("--pivot","0: none, 1: partial, 2: full",1);
        const bool forceGrowth = Input
            ("--forceGrowth","force element growth?",false);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();
        if( pivot < 0 || pivot > 2 )
            LogicError("Invalid pivot value");

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
        {
            if( pivot == 0 )
                Output("Testing LU with no pivoting");
            else if( pivot == 1 )
                Output("Testing LU with partial pivoting");
            else if( pivot == 2 )
                Output("Testing LU with full pivoting");
        }

        if( commRank == 0 )
            Output("Testing with doubles:");
        TestLU<double,STAR>( m, g, pivot, testCorrectness, forceGrowth, print );

        if( commRank == 0 )
            Output("Testing with double-precision complex:");
        TestLU<Complex<double>,STAR>
        ( m, g, pivot, testCorrectness, forceGrowth, print );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
