/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/LU.hpp"
#include "elemental/lapack-like/Norm.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( bool pivoted, bool print, 
  const DistMatrix<F>& A,
  const DistMatrix<int,VC,STAR>& p,
  const DistMatrix<F>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Generate random right-hand sides
    DistMatrix<F> X(g), Y(g);
    Uniform( m, 100, X );
    Y = X;
    if( pivoted )
        ApplyRowPivots( Y, p );

    // Solve against the (pivoted) right-hand sides
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, Y );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    const R oneNormOfX = Norm( X, ONE_NORM );
    const R infNormOfX = Norm( X, INFINITY_NORM );
    const R frobNormOfX = Norm( X, FROBENIUS_NORM );
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const R oneNormOfError = Norm( X, ONE_NORM );
    const R infNormOfError = Norm( X, INFINITY_NORM );
    const R frobNormOfError = Norm( X, FROBENIUS_NORM );
    const R oneNormOfA = Norm( AOrig, ONE_NORM );
    const R infNormOfA = Norm( AOrig, INFINITY_NORM );
    const R frobNormOfA = Norm( AOrig, FROBENIUS_NORM );

    if( g.Rank() == 0 )
    {
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
}

template<typename F> 
void TestLU
( bool pivot, bool testCorrectness, bool print, 
  int m, const Grid& g )
{
    DistMatrix<F> A(g), ARef(g);
    DistMatrix<int,VC,STAR> p(g);

    Uniform( m, m, A );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        ARef = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( print )
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting LU factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( pivot )
        LU( A, p );
    else
        LU( A );
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
        int r = Input("--gridHeight","height of process grid",0);
        const int m = Input("--height","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool pivot = Input("--pivot","pivoted LU?",true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
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
        TestLU<double>( pivot, testCorrectness, print, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestLU<Complex<double> >( pivot, testCorrectness, print, m, g );
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
