/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
using namespace El;

template<typename F,Dist UPerm> 
void TestCorrectness
( const DistMatrix<F>& AOrig,
  const DistMatrix<F>& A,
  const DistMatrix<Int,UPerm,STAR>& p,
  const DistMatrix<Int,UPerm,STAR>& q,
  Int pivoting, bool print )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Generate random right-hand sides
    DistMatrix<F> X(g);
    Uniform( X, m, 100 );
    auto Y( X );
    if( pivoting == 0 )
        lu::SolveAfter( NORMAL, A, Y );
    else if( pivoting == 1 )
        lu::SolveAfter( NORMAL, A, p, Y );
    else
        lu::SolveAfter( NORMAL, A, p, q, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    const Real oneNormOfX = OneNorm( X );
    const Real infNormOfX = InfinityNorm( X );
    const Real frobNormOfX = FrobeniusNorm( X );
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real oneNormOfError = OneNorm( X );
    const Real infNormOfError = InfinityNorm( X );
    const Real frobNormOfError = FrobeniusNorm( X );
    const Real oneNormOfA = OneNorm( AOrig );
    const Real infNormOfA = InfinityNorm( AOrig );
    const Real frobNormOfA = FrobeniusNorm( AOrig );

    if( g.Rank() == 0 )
    {
        cout << "||A||_1             = " << oneNormOfA << "\n"
             << "||A||_oo            = " << infNormOfA << "\n"
             << "||A||_F             = " << frobNormOfA << "\n"
             << "||X||_1             = " << oneNormOfX << "\n"
             << "||X||_oo            = " << infNormOfX << "\n"
             << "||X||_F             = " << frobNormOfX << "\n"
             << "||A A^-1 X - X||_1  = " << oneNormOfError << "\n"
             << "||A A^-1 X - X||_oo = " << infNormOfError << "\n"
             << "||A A^-1 X - X||_F  = " << frobNormOfError << endl;
    }
}

template<typename F,Dist UPerm> 
void TestLU
( Int m, const Grid& g, Int pivoting, 
  bool testCorrectness, bool forceGrowth, bool print )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<Int,UPerm,STAR> p(g), q(g);

    if( forceGrowth )
        GEPPGrowth( A, m );
    else
        Uniform( A, m, m );

    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( print )
        Print( A, "A" );

    if( g.Rank() == 0 )
    {
        cout << "  Starting LU factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( pivoting == 0 )
        LU( A );
    else if( pivoting == 1 )
        LU( A, p );
    else if( pivoting == 2 )
        LU( A, p, q );

    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivoting >= 1 )
        {
            Print( p, "p after factorization");
            DistMatrix<Int> P(g);
            ExplicitPermutation( p, P );
            Print( P, "P" );
        }
        if( pivoting == 2 )
        {
            Print( q, "q after factorization");
            DistMatrix<Int> Q(g);
            ExplicitPermutation( q, Q );
            Print( Q, "Q" );
        }
    }
    if( testCorrectness )
        TestCorrectness( AOrig, A, p, q, pivoting, print );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
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
            cout << "Will test LU with ";
            if( pivot == 0 )
                cout << "no pivoting" << std::endl;
            else if( pivot == 1 )
                cout << "partial pivoting" << std::endl;
            else if( pivot == 2 )
                cout << "full pivoting" << std::endl;
        }

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestLU<double,VC>( m, g, pivot, testCorrectness, forceGrowth, print );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestLU<Complex<double>,VC>
        ( m, g, pivot, testCorrectness, forceGrowth, print );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
