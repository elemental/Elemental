/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_HEMM_INC
#include ELEM_TRMM_INC
#include ELEM_CHOLESKY_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_ONENORM_INC
using namespace std;
using namespace elem;

template<typename F>
void TestCorrectness
( bool pivot, UpperOrLower uplo,
  const DistMatrix<F>& A,
  const DistMatrix<Int,VC,STAR>& p,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    // Test correctness by multiplying a random set of vectors by A, then
    // using the Cholesky factorization to solve.
    auto X = Uniform<F>( g, m, 100 );
    auto Y = Zeros<F>( g, m, 100 );
    Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
    const Real maxNormL = HermitianMaxNorm( uplo, A );
    const Real maxNormA = HermitianMaxNorm( uplo, AOrig );
    const Real infNormA = HermitianInfinityNorm( uplo, AOrig );
    const Real frobNormA = HermitianFrobeniusNorm( uplo, AOrig );
    const Real oneNormY = OneNorm( Y );
    const Real infNormY = InfinityNorm( Y );
    const Real frobNormY = FrobeniusNorm( Y );

    if( pivot )
        cholesky::SolveAfter( uplo, NORMAL, A, p, Y );
    else
        cholesky::SolveAfter( uplo, NORMAL, A, Y );
    Axpy( F(-1), Y, X );
    const Real oneNormE = OneNorm( X );
    const Real infNormE = InfinityNorm( X );
    const Real frobNormE = FrobeniusNorm( X );

    if( g.Rank() == 0 )
    {
        cout << "||L||_max            = " << maxNormL << "\n"
             << "||A||_max            = " << maxNormA << "\n"
             << "||A||_1 = ||A||_oo   = " << infNormA << "\n"
             << "||A||_F              = " << frobNormA << "\n"
             << "||Y||_1              = " << oneNormY << "\n"
             << "||Y||_oo             = " << infNormY << "\n"
             << "||Y||_F              = " << frobNormY << "\n"
             << "||X - inv(A) X||_1  = " << oneNormE << "\n"
             << "||X - inv(A) X||_oo = " << infNormE << "\n"
             << "||X - inv(A) X||_F  = " << frobNormE << endl;
    }
}

template<typename F> 
void TestCholesky
( bool testCorrectness, bool pivot, bool unblocked, bool print, bool printDiag,
  UpperOrLower uplo, Int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<Int,VC,STAR> p(g);

    HermitianUniformSpectrum( A, m, 1e-9, 10 );
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
        cout << "  Starting Cholesky factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( pivot )
    {
        if( unblocked )
        {
            if( uplo == LOWER )
                cholesky::LUnblockedPivoted( A, p );
            else
                cholesky::UUnblockedPivoted( A, p );
        }
        else
            Cholesky( uplo, A, p );
    }
    else
        Cholesky( uplo, A );
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
        if( pivot )
            Print( p, "p" );
    }
    if( printDiag )
        Print( A.GetRealPartOfDiagonal(), "diag(A)" );
    if( testCorrectness )
        TestCorrectness( pivot, uplo, A, p, AOrig );
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
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--m","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool pivot = Input("--pivot","use pivoting?",false);
        const bool unblocked = Input("--unblocked","unblocked pivoting?",false);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool printDiag = Input("--printDiag","print diag of fact?",false);
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
            cout << "Will test Cholesky" << uploChar << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestCholesky<double>
        ( testCorrectness, pivot, unblocked, print, printDiag, uplo, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestCholesky<Complex<double>>
        ( testCorrectness, pivot, unblocked, print, printDiag, uplo, m, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
