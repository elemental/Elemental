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

template<typename F> 
void TestTrsv
( bool print, UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  int n, const Grid& g )
{
    typedef typename Base<F>::type R;
    DistMatrix<F> A(g), x(g), y(g);

    // Generate random A and x
    HermitianUniformSpectrum( n, A, 1, 10 );
    Uniform( n, 1, x );
    // Either y := op(L) x or y := op(U) x
    y = x;
    Trmm( LEFT, uplo, orientation, diag, F(1), A, y );

    if( print )
    {
        A.Print("A");
        x.Print("x");
        y.Print("y");
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
        y.Print("y after solve");

    Axpy( F(-1), x, y );
    const R xNorm = Norm( x );
    const R yNorm = Norm( y );
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
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try
    {
        int r = Input("--r","height of process grid",0);
        const char uploChar = Input
            ("--uplo","upper or lower triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of triangular matrix: N/T/C",'N');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const int n = Input("--n","size of triangular matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
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
        TestTrsv<Complex<double> >( print, uplo, orientation, diag, n, g );
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
