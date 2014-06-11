/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_UNIFORM_INC
using namespace std;
using namespace El;

template<typename F> 
void TestCorrectness
( bool pivoted, bool print, 
  const Matrix<F>& A,
  const Matrix<Int>& pPerm,
  const Matrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Int m = AOrig.Height();

    cout << "Testing error..." << endl;

    // Generate random right-hand sides
    Matrix<F> X, Y;
    Uniform( X, m, 100 );
    Y = X;
    if( pivoted )
        PermuteRows( Y, pPerm );

    // Solve against the (pivoted) right-hand sides
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, Y );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, Y );

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

    cout << "||A||_1                  = " << oneNormOfA << "\n"
         << "||A||_oo                 = " << infNormOfA << "\n"
         << "||A||_F                  = " << frobNormOfA << "\n"
         << "||X||_1                  = " << oneNormOfX << "\n"
         << "||X||_oo                 = " << infNormOfX << "\n"
         << "||X||_F                  = " << frobNormOfX << "\n"
         << "||A U^-1 L^-1 X - X||_1  = " << oneNormOfError << "\n"
         << "||A U^-1 L^-1 X - X||_oo = " << infNormOfError << "\n"
         << "||A U^-1 L^-1 X - X||_F  = " << frobNormOfError << endl;
}

template<typename F> 
void TestLU( bool pivot, bool testCorrectness, bool print, Int m )
{
    Matrix<F> A, ARef;
    Matrix<Int> pPerm;

    Uniform( A, m, m );
    if( testCorrectness )
    {
        cout << "  Making copy of original matrix...";
        cout.flush();
        ARef = A;
        cout << "DONE" << endl;
    }
    if( print )
        Print( A, "A" );

    cout << "  Starting LU factorization...";
    cout.flush();
    const double startTime = mpi::Time();
    if( pivot )
        LU( A, pPerm );
    else
        LU( A );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    cout << "DONE. " << endl
         << "  Time = " << runTime << " seconds. GFlops = " 
         << gFlops << endl;
    if( print )
    {
        Print( A, "A after factorization" );
        if( pivot )
            Print( pPerm, "pPerm after factorization" );
    }
    if( testCorrectness )
        TestCorrectness( pivot, print, A, pPerm, ARef );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

    try
    {
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool pivot = Input("--pivot","pivoted LU?",true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test LU" 
                 << ( pivot ? " with partial pivoting" : " " ) << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestLU<double>( pivot, testCorrectness, print, m );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestLU<Complex<double>>( pivot, testCorrectness, print, m );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
