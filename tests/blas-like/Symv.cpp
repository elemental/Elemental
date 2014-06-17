/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"

using namespace std;
using namespace El;

template<typename T> 
void TestSymv
( const UpperOrLower uplo, const Int m, const T alpha, const T beta, 
  const bool print, const Grid& g )
{
    DistMatrix<T> A(g), x(g), y(g);

    Uniform( A, m, m );
    Uniform( x, m, 1 );
    Uniform( y, m, 1 ); 
    if( print )
    {
        Print( A, "A" );
        Print( x, "x" );
        Print( y, "y" );
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
    const double gFlops = ( IsComplex<T>::val ? 4*realGFlops : realGFlops );
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
        Print( y, msg.str() );
    }
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
        Int r = Input("--r","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--m","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocalDouble = Input
            ("--nbLocalDouble","local blocksize for real doubles",32);
        const Int nbLocalComplexDouble = Input
            ("--nbLocalComplexDouble","local blocksize for complex doubles",32);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalSymvBlocksize<double>( nbLocalDouble );
        SetLocalSymvBlocksize<Complex<double>>( nbLocalComplexDouble );

        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test Symv" << uploChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestSymv<double>( uplo, m, 3., 4., print, g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestSymv<Complex<double>>
        ( uplo, m, Complex<double>(3), Complex<double>(4), print, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
