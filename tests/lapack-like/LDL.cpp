/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_HERMITIANUNIFORMSPECTRUM_INC
using namespace std;
using namespace El;

template<typename F,Dist UPerm> 
void TestCorrectness
( bool conjugated, bool print, 
  const DistMatrix<F>& A,
  const DistMatrix<F,MD,STAR>& dSub,
  const DistMatrix<Int,UPerm,STAR>& pPerm,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    DistMatrix<F> X(g), Y(g);
    Uniform( X, m, 100 );
    Y = X;

    // Test correctness by comparing the application of AOrig against a 
    // random set of 100 vectors to the application of tril(A) tril(A)^H
    if( print )
        Print( X, "X" );
    ldl::MultiplyAfter( A, dSub, pPerm, Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X" );
    Symm( LEFT, LOWER, F(-1), AOrig, X, F(1), Y, conjugated );
    if( print )
        Print( Y, "P' L B L' P X - A X" );
    const Real oneNormOfError = OneNorm( Y );
    const Real infNormOfError = InfinityNorm( Y );
    const Real frobNormOfError = FrobeniusNorm( Y );
    const Real infNormOfA = HermitianInfinityNorm( LOWER, AOrig );
    const Real frobNormOfA = HermitianFrobeniusNorm( LOWER, AOrig );
    const Real oneNormOfX = OneNorm( X );
    const Real infNormOfX = InfinityNorm( X );
    const Real frobNormOfX = FrobeniusNorm( X );
    if( g.Rank() == 0 )
    {
        cout << "||A||_1 = ||A||_oo   = " << infNormOfA << "\n"
             << "||A||_F              = " << frobNormOfA << "\n"
             << "||X||_1              = " << oneNormOfX << "\n"
             << "||X||_oo             = " << infNormOfX << "\n"
             << "||X||_F              = " << frobNormOfX << "\n"
             << "||A X - L D L^[T/H] X||_1  = " << oneNormOfError << "\n"
             << "||A X - L D L^[T/H] X||_oo = " << infNormOfError << "\n"
             << "||A X - L D L^[T/H] X||_F  = " << frobNormOfError << endl;
    }
}

template<typename F,Dist UPerm> 
void TestLDL
( bool conjugated, bool testCorrectness, bool print, 
  Int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    if( conjugated )
        HermitianUniformSpectrum( A, m, -100, 100 );
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
        cout << "  Starting LDL^[T/H] factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    DistMatrix<F,MD,STAR> dSub(g);
    DistMatrix<Int,UPerm,STAR> pPerm(g);
    LDL( A, dSub, pPerm, conjugated );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE.\n"
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        Print( A, "A after factorization" );
        Print( pPerm, "pPerm" );
    }
    if( testCorrectness )
        TestCorrectness( conjugated, print, A, dSub, pPerm, AOrig );
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
        Int r = Input("--gridHeight","process grid height",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool conjugated = Input("--conjugate","conjugate LDL?",false);
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
        SetLocalTrrkBlocksize<double>( nbLocal );
        SetLocalTrrkBlocksize<Complex<double>>( nbLocal );
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test LDL" << (conjugated?"^H":"^T") << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestLDL<double,VC>( conjugated, testCorrectness, print, m, g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestLDL<Complex<double>,VC>( conjugated, testCorrectness, print, m, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
