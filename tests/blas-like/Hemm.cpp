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
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
using namespace std;
using namespace elem;

template<typename T> 
void TestHemm
( bool print, LeftOrRight side, UpperOrLower uplo,
  Int m, Int n, T alpha, T beta, const Grid& g )
{
    DistMatrix<T> A(g), B(g), C(g);

    if( side == LEFT )
        HermitianUniformSpectrum( A, m, -10, 10 );
    else
        HermitianUniformSpectrum( A, n, -10, 10 );
    Uniform( B, m, n );
    Uniform( C, m, n );
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( C, "C" );
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting Parallel Hemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Hemm( side, uplo, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double realGFlops = 
        ( side==LEFT ? 2.*mD*mD*nD : 2.*mD*nD*nD ) / (1.e9*runTime);
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
        if( side == LEFT )
            msg << "C := " << alpha << " Herm(A) B + " << beta << " C";
        else
            msg << "C := " << alpha << " B Herm(A) + " << beta << " C";
        Print( C, msg.str() );
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
        const char sideChar = Input("--side","side to apply from: L/R",'L');
        const char uploChar = Input("--uplo","lower/upper storage: L/U",'L');
        const Int m = Input("--m","height of result",100);
        const Int n = Input("--n","width of result",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test Hemm" << sideChar << uploChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestHemm<double>( print, side, uplo, m, n, 3., 4., g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestHemm<Complex<double>>
        ( print, side, uplo, m, n, Complex<double>(3), Complex<double>(4), g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
