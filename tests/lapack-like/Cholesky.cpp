/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <ctime>
#include "elemental.hpp"
using namespace std;
using namespace elem;

template<typename F>
void TestCorrectness
( bool printMatrices, UpperOrLower uplo,
  const DistMatrix<F>& A,
  const DistMatrix<F>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F> X(g), Y(g);
    Uniform( m, 100, X );
    Y = X;

    if( uplo == LOWER )
    {
        // Test correctness by comparing the application of AOrig against a 
        // random set of 100 vectors to the application of tril(A) tril(A)^H
        Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, Y );
        Trmm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, Y );
        Hemm( LEFT, LOWER, F(-1), AOrig, X, F(1), Y );
        const R oneNormOfError = Norm( Y, ONE_NORM );
        const R infNormOfError = Norm( Y, INFINITY_NORM );
        const R frobNormOfError = Norm( Y, FROBENIUS_NORM );
        const R infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        const R frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        const R oneNormOfX = Norm( X, ONE_NORM );
        const R infNormOfX = Norm( X, INFINITY_NORM );
        const R frobNormOfX = Norm( X, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "||A||_1 = ||A||_oo   = " << infNormOfA << "\n"
                 << "||A||_F              = " << frobNormOfA << "\n"
                 << "||X||_1              = " << oneNormOfX << "\n"
                 << "||X||_oo             = " << infNormOfX << "\n"
                 << "||X||_F              = " << frobNormOfX << "\n"
                 << "||A X - L L^H X||_1  = " << oneNormOfError << "\n"
                 << "||A X - L L^H X||_oo = " << infNormOfError << "\n"
                 << "||A X - L L^H X||_F  = " << frobNormOfError << endl;
        }
    }
    else
    {
        // Test correctness by comparing the application of AOrig against a 
        // random set of 100 vectors to the application of triu(A)^H triu(A)
        Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, Y );
        Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, Y );
        Hemm( LEFT, UPPER, F(-1), AOrig, X, F(1), Y );
        const R oneNormOfError = Norm( Y, ONE_NORM );
        const R infNormOfError = Norm( Y, INFINITY_NORM );
        const R frobNormOfError = Norm( Y, FROBENIUS_NORM );
        const R infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        const R frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        const R oneNormOfX = Norm( X, ONE_NORM );
        const R infNormOfX = Norm( X, INFINITY_NORM );
        const R frobNormOfX = Norm( X, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "||A||_1 = ||A||_oo   = " << infNormOfA << "\n"
                 << "||A||_F              = " << frobNormOfA << "\n"
                 << "||X||_1              = " << oneNormOfX << "\n"
                 << "||X||_oo             = " << infNormOfX << "\n"
                 << "||X||_F              = " << frobNormOfX << "\n"
                 << "||A X - U^H U X||_1  = " << oneNormOfError << "\n"
                 << "||A X - U^H U X||_oo = " << infNormOfError << "\n"
                 << "||A X - U^H U X||_F  = " << frobNormOfError << endl;
        }
    }
}

template<typename F> 
void TestCholesky
( bool testCorrectness, bool printMatrices, 
  UpperOrLower uplo, int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);

    HermitianUniformSpectrum( m, A, 1, 10 );
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
    if( printMatrices )
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting Cholesky factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
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
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, uplo, A, AOrig );
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
        int r = Input("--gridHeight","process grid height",0);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const int m = Input("--height","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool printMatrices = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalTrrkBlocksize<double>( nbLocal );
        SetLocalTrrkBlocksize<Complex<double> >( nbLocal );
#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test Cholesky" << uploChar << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestCholesky<double>
        ( testCorrectness, printMatrices, uplo, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestCholesky<Complex<double> >
        ( testCorrectness, printMatrices, uplo, m, g );
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
