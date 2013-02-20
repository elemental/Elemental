/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/LU.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( bool pivoted, bool print, 
  const Matrix<F>& A,
  const Matrix<int>& p,
  const Matrix<F>& AOrig )
{
    typedef typename Base<F>::type R;
    const int m = AOrig.Height();

    cout << "Testing error..." << endl;

    // Generate random right-hand sides
    Matrix<F> X, Y;
    Uniform( m, 100, X );
    Y = X;
    if( pivoted )
        ApplyRowPivots( Y, p );

    // Solve against the (pivoted) right-hand sides
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, Y );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    const R oneNormOfX = OneNorm( X );
    const R infNormOfX = InfinityNorm( X );
    const R frobNormOfX = FrobeniusNorm( X );
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const R oneNormOfError = OneNorm( X );
    const R infNormOfError = InfinityNorm( X );
    const R frobNormOfError = FrobeniusNorm( X );
    const R oneNormOfA = OneNorm( AOrig );
    const R infNormOfA = InfinityNorm( AOrig );
    const R frobNormOfA = FrobeniusNorm( AOrig );

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
void TestLU( bool pivot, bool testCorrectness, bool print, int m )
{
    Matrix<F> A, ARef;
    Matrix<int> p;

    Uniform( m, m, A );
    if( testCorrectness )
    {
        cout << "  Making copy of original matrix...";
        cout.flush();
        ARef = A;
        cout << "DONE" << endl;
    }
    if( print )
        A.Print("A");

    cout << "  Starting LU factorization...";
    cout.flush();
    const double startTime = mpi::Time();
    if( pivot )
        LU( A, p );
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
        A.Print("A after factorization");
        if( pivot )
            p.Print("p after factorization");
    }
    if( testCorrectness )
        TestCorrectness( pivot, print, A, p, ARef );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try
    {
        const int m = Input("--height","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool pivot = Input("--pivot","pivoted LU?",true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test LU" 
                 << ( pivot ? " with partial pivoting" : " " ) << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestLU<double>( pivot, testCorrectness, print, m );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestLU<Complex<double> >( pivot, testCorrectness, print, m );
    }
    catch( ArgException& e ) { }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }   
    Finalize();
    return 0;
}
