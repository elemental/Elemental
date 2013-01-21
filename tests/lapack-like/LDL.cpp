/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Hemm.hpp"
#include "elemental/blas-like/level3/Symm.hpp"
#include "elemental/blas-like/level3/Trmm.hpp"
#include "elemental/lapack-like/HermitianNorm.hpp"
#include "elemental/lapack-like/LDL.hpp"
#include "elemental/lapack-like/Norm.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( bool conjugated, bool print, 
  const DistMatrix<F>& A,
  const DistMatrix<F,MC,STAR>& d,
  const DistMatrix<F>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F> X(g), Y(g);
    Uniform( m, 100, X );
    Y = X;

    // Test correctness by comparing the application of AOrig against a 
    // random set of 100 vectors to the application of tril(A) tril(A)^H
    if( conjugated )
        Trmm( LEFT, LOWER, ADJOINT, UNIT, F(1), A, Y );
    else
        Trmm( LEFT, LOWER, TRANSPOSE, UNIT, F(1), A, Y );
    DiagonalScale( LEFT, NORMAL, d, Y );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, Y );
    if( conjugated )
        Hemm( LEFT, LOWER, F(-1), AOrig, X, F(1), Y );
    else
        Symm( LEFT, LOWER, F(-1), AOrig, X, F(1), Y );
    const R oneNormOfError = Norm( Y, ONE_NORM );
    const R infNormOfError = Norm( Y, INFINITY_NORM );
    const R frobNormOfError = Norm( Y, FROBENIUS_NORM );
    const R infNormOfA = HermitianNorm( LOWER, AOrig, INFINITY_NORM );
    const R frobNormOfA = HermitianNorm( LOWER, AOrig, FROBENIUS_NORM );
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
             << "||A X - L D L^[T/H] X||_1  = " << oneNormOfError << "\n"
             << "||A X - L D L^[T/H] X||_oo = " << infNormOfError << "\n"
             << "||A X - L D L^[T/H] X||_F  = " << frobNormOfError << endl;
    }
}

template<typename F> 
void TestLDL
( bool conjugated, bool testCorrectness, bool print, 
  int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    if( conjugated )
        HermitianUniformSpectrum( m, A, -100, 100 );
    else
        Uniform( m, m, A );
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
        A.Print("A");
    DistMatrix<F,MC,STAR> d(g);

    if( g.Rank() == 0 )
    {
        cout << "  Starting LDL^[T/H] factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( !conjugated )
        LDLT( A, d );
    else
        LDLH( A, d );
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
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( conjugated, print, A, d, AOrig );
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
        const int m = Input("--height","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool conjugated = Input("--conjugate","conjugate LDL?",false);
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
            cout << "Will test LDL" << (conjugated?"^H":"^T") << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestLDL<double>( conjugated, testCorrectness, print, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestLDL<Complex<double> >( conjugated, testCorrectness, print, m, g );
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
