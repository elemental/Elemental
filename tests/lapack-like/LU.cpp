/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_GEMM_INC
#include ELEM_LU_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_UNIFORM_INC
#include ELEM_EXPLICITPERMUTATION_INC
using namespace std;
using namespace elem;

template<typename F,Dist UPerm> 
void TestCorrectness
( Int pivoting, bool print, 
  const DistMatrix<F>& A,
  const DistMatrix<Int,UPerm,STAR>& pPerm,
  const DistMatrix<Int,UPerm,STAR>& qPerm,
  const DistMatrix<F>& AOrig )
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
        lu::SolveAfter( NORMAL, A, pPerm, Y );
    else
        lu::SolveAfter( NORMAL, A, pPerm, qPerm, Y );

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
( Int pivoting, bool testCorrectness, bool print, 
  Int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<Int,UPerm,STAR> pPerm(g), qPerm(g);

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
        LU( A, pPerm );
    else if( pivoting == 2 )
        LU( A, pPerm, qPerm );

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
            Print( pPerm, "pPerm after factorization");
            DistMatrix<Int> P(g);
            ExplicitPermutation( pPerm, P );
            Print( P, "P" );
        }
        if( pivoting == 2 )
        {
            Print( qPerm, "qPerm after factorization");
            DistMatrix<Int> Q(g);
            ExplicitPermutation( qPerm, Q );
            Print( Q, "Q" );
        }
    }
    if( testCorrectness )
        TestCorrectness( pivoting, print, A, pPerm, qPerm, AOrig );
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
        TestLU<double,MC>( pivot, testCorrectness, print, m, g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestLU<Complex<double>,MC>( pivot, testCorrectness, print, m, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
