/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level2/Trsv.hpp"
#include "elemental/blas-like/level3/Trmm.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
using namespace std;
using namespace elem;

template<typename F> 
void TestTrsv
( bool print, UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  Int n, const Grid& g )
{
    typedef BASE(F) Real;
    DistMatrix<F> A(g), x(g), y(g);

    // Generate random A and x
    HermitianUniformSpectrum( A, n, 1, 10 );
    Uniform( x, n, 1 );
    // Either y := op(L) x or y := op(U) x
    y = x;
    Trmm( LEFT, uplo, orientation, diag, F(1), A, y );

    if( print )
    {
        Print( A, "A" );
        Print( x, "x" );
        Print( y, "y" );
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Trsv...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Trsv( uplo, orientation, diag, A, y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = Pow(double(n),2.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
        Print( y, "y after solve" );

    Axpy( F(-1), x, y );
    const Real xNorm = FrobeniusNorm( x );
    const Real yNorm = FrobeniusNorm( y );
    if( g.Rank() == 0 )
    {
        std::cout << "|| x - y ||_2 = " << yNorm << "\n"
                  << "|| x ||_2     = " << xNorm << "\n"
                  << "|| x - y ||_2 / || x ||_2 = " << yNorm/xNorm << "\n"
                  << std::endl;
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::CommRank( comm );
    const Int commSize = mpi::CommSize( comm );

    try
    {
        Int r = Input("--r","height of process grid",0);
        const char uploChar = Input
            ("--uplo","upper or lower triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of triangular matrix: N/T/C",'N');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const Int n = Input("--n","size of triangular matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const Grid g( comm, r );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test Trsv" << uploChar << transChar << diagChar 
                 << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTrsv<double>( print, uplo, orientation, diag, n, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestTrsv<Complex<double>>( print, uplo, orientation, diag, n, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
