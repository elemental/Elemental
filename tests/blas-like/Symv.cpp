/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

template<typename T> 
void TestSymv
( const UpperOrLower uplo, const int m, const T alpha, const T beta, 
  const bool print, const Grid& g )
{
    DistMatrix<T> A(g), x(g), y(g);

    Uniform( m, m, A );
    Uniform( m, 1, x );
    Uniform( m, 1, y ); 
    if( print )
    {
        A.Print("A");
        x.Print("x");
        y.Print("y");
    }

    // Test Symm
    if( g.Rank() == 0 )
    {
        cout << "  Starting Parallel Symv...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Symv( uplo, alpha, A, x, beta, y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2.*double(m)*double(m)/(1.e9*runTime);
    const double gFlops = ( IsComplex<T>::val ? 4*realGFlops : gFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        ostringstream msg;
        msg << "y := " << alpha << " Symm(A) x + " << beta << " y";
        y.Print( msg.str() );
    }
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
        int r = Input("--r","height of process grid",0);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const int m = Input("--m","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const int nbLocalDouble = Input
            ("--nbLocalDouble","local blocksize for real doubles",32);
        const int nbLocalComplexDouble = Input
            ("--nbLocalComplexDouble","local blocksize for complex doubles",32);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalSymvBlocksize<double>( nbLocalDouble );
        SetLocalSymvBlocksize<Complex<double> >( nbLocalComplexDouble );

#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test Symv" << uploChar << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestSymv<double>( uplo, m, (double)3, (double)4, print, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestSymv<Complex<double> >
        ( uplo, m, Complex<double>(3), Complex<double>(4), print, g );
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
