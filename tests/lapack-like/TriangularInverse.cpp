/*
   Copyright (c) 2009-2012, Jack Poulson
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
( bool print,
  UpperOrLower uplo, UnitOrNonUnit diag,
  const DistMatrix<F>& A,
  const DistMatrix<F>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F> X(g), Y(g);
    Uniform( m, 100, X );
    Y = X;

    // Since A o A^-1 = I, test the change introduced by the approximate comp.
    Trmm( LEFT, uplo, NORMAL, diag, F(1), A,     Y );
    Trmm( LEFT, uplo, NORMAL, diag, F(1), AOrig, Y );
    Axpy( F(-1), X, Y );

    const R oneNormOrig = Norm( AOrig, ONE_NORM );
    const R infNormOrig = Norm( AOrig, INFINITY_NORM );
    const R frobNormOrig = Norm( AOrig, FROBENIUS_NORM );
    const R oneNormFinal = Norm( A, ONE_NORM );
    const R infNormFinal = Norm( A, INFINITY_NORM );
    const R frobNormFinal = Norm( A, FROBENIUS_NORM );
    const R oneNormOfError = Norm( Y, ONE_NORM );
    const R infNormOfError = Norm( Y, INFINITY_NORM );
    const R frobNormOfError = Norm( Y, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "||A||_1           = " << oneNormOrig << "\n"
             << "||A||_oo          = " << infNormOrig << "\n"
             << "||A||_F           = " << frobNormOrig << "\n"
             << "||A^-1||_1        = " << oneNormFinal << "\n"
             << "||A^-1||_oo       = " << infNormFinal << "\n"
             << "||A^-1||_F        = " << frobNormFinal << "\n"
             << "||A A^-1 - I||_1  = " << oneNormOfError << "\n"
             << "||A A^-1 - I||_oo = " << infNormOfError << "\n"
             << "||A A^-1 - I||_F  = " << frobNormOfError << endl;
    }
}

template<typename F> 
void TestTriangularInverse
( bool testCorrectness, bool print,
  UpperOrLower uplo, UnitOrNonUnit diag, int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    HermitianUniformSpectrum( m, A, 1, 10 );
    MakeTriangular( uplo, A );
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

    if( g.Rank() == 0 )
    {
        cout << "  Starting triangular inversion...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TriangularInverse( uplo, diag, A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
        A.Print("A after inversion");
    if( testCorrectness )
        TestCorrectness( print, uplo, diag, A, AOrig );
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
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const int m = Input("--height","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
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
            cout << "Will test TriangularInverse" << uploChar << diagChar 
                 << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTriangularInverse<double>
        ( testCorrectness, print, uplo, diag, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestTriangularInverse<Complex<double> >
        ( testCorrectness, print, uplo, diag, m, g );
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
