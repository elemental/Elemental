/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

template<typename T> 
void TestGemm
( bool print, Orientation orientA, Orientation orientB,
  int m, int n, int k, T alpha, T beta, const Grid& g )
{
    double startTime, endTime, runTime, realGFlops, gFlops;
    DistMatrix<T> A(g), B(g), C(g);

    if( orientA == NORMAL )
        A.ResizeTo( m, k );
    else
        A.ResizeTo( k, m );
    if( orientB == NORMAL )
        B.ResizeTo( k, n );
    else
        B.ResizeTo( n, k );
    C.ResizeTo( m, n );

    // Test the variant of Gemm that keeps A stationary
    if( g.Rank() == 0 )
        cout << "Stationary A Algorithm:" << endl;
    MakeUniform( A );
    MakeUniform( B );
    MakeUniform( C );
    if( print )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Gemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    internal::GemmA( orientA, orientB, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    runTime = mpi::Time() - startTime;
    realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
    gFlops = ( IsComplex<T>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }

    // Test the variant of Gemm that keeps B stationary
    if( g.Rank() == 0 )
        cout << endl << "Stationary B Algorithm:" << endl;
    MakeUniform( A );
    MakeUniform( B );
    MakeUniform( C );
    if( print )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Gemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    internal::GemmB( orientA, orientB, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    runTime = mpi::Time() - startTime;
    realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
    gFlops = ( IsComplex<T>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl 
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }

    // Test the variant of Gemm that keeps C stationary
    if( g.Rank() == 0 )
        cout << endl << "Stationary C Algorithm:" << endl;
    MakeUniform( A );
    MakeUniform( B );
    MakeUniform( C );
    if( print )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Gemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    internal::GemmC( orientA, orientB, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    runTime = mpi::Time() - startTime;
    realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
    gFlops = ( IsComplex<T>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    
    if( orientA == NORMAL && orientB == NORMAL )
    {
        // Test the variant of Gemm for panel-panel dot products
        if( g.Rank() == 0 )
            cout << endl << "Dot Product Algorithm:" << endl;
        MakeUniform( A );
        MakeUniform( B );
        MakeUniform( C );
        if( print )
        {
            A.Print("A");
            B.Print("B");
            C.Print("C");
        }
        if( g.Rank() == 0 )
        {
            cout << "  Starting Gemm...";
            cout.flush();
        }
        mpi::Barrier( g.Comm() );
        startTime = mpi::Time();
        internal::GemmDot( orientA, orientB, alpha, A, B, beta, C );
        mpi::Barrier( g.Comm() );
        runTime = mpi::Time() - startTime;
        realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
        gFlops = ( IsComplex<T>::val ? 4*realGFlops : realGFlops );
        if( g.Rank() == 0 )
        {
            cout << "DONE. " << endl
                 << "  Time = " << runTime << " seconds. GFlops = " 
                 << gFlops << endl;
        }
        if( print )
        {
            ostringstream msg;
            msg << "C := " << alpha << " A B + " << beta << " C";
            C.Print( msg.str() );
        }
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
        const char transA = Input("--transA","orientation of A: N/T/C",'N');
        const char transB = Input("--transB","orientation of B: N/T/C",'N');
        const int m = Input("--m","height of result",100);
        const int n = Input("--n","width of result",100);
        const int k = Input("--k","inner dimension",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const Orientation orientA = CharToOrientation( transA );
        const Orientation orientB = CharToOrientation( transB );
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
            cout << "Will test Gemm" << transA << transB << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestGemm<double>( print, orientA, orientB, m, n, k, 3., 4., g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestGemm<Complex<double> >
        ( print, orientA, orientB, m, n, k, 
          Complex<double>(3), Complex<double>(4), g );
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
