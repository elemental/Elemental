/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include <ctime>
#include "elemental.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Generates random matrix then solves for its LU factors."
         << endl << endl;
    cout << "  LU <r> <c> <m> <nb> <test correctness?> "
         << "<print matrices?>" << endl << endl;
    cout << "  r: number of process rows      " << endl;
    cout << "  c: number of process cols      " << endl;
    cout << "  m: height of matrix            " << endl;
    cout << "  nb: algorithmic blocksize      " << endl;
    cout << "  test correctness?: false iff 0 " << endl;
    cout << "  print matrices?: false iff 0   " << endl;
    cout << endl;
}

template<typename T>
bool OKRelativeError( T truth, T computed );

template<>
bool OKRelativeError( double truth, double computed )
{ return ( fabs(truth-computed) / max(fabs(truth),(double)1) <= 1e-10 ); }

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{ return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-10 ); }
#endif

template<typename T>
void TestCorrectness
( bool printMatrices,
  const DistMatrix<T,MC,MR>& A,
  const DistMatrix<int,VC,Star>& p,
        DistMatrix<T,Star,Star>& ARef )
{
    const Grid& g = A.GetGrid();
    const int m = ARef.Height();
    DistMatrix<T,Star,Star> A_copy(g);
    DistMatrix<int,Star,Star> p_copy(g);

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    A_copy = A;
    p_copy = p;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    DistMatrix<int,Star,Star> pRef(m,1,g);
    lapack::LU( ARef.LocalMatrix(), pRef.LocalMatrix() );
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
    {
        ARef.Print("True A:");
        pRef.Print("True p:");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    // Test A
    for( int j=0; j<m; ++j )
    {
        for( int i=0; i<m; ++i )
        {
            T truth = ARef.LocalEntry(i,j);
            T computed = A_copy.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                ostringstream msg;
                msg << "FAILED at index (" << i << "," << j << "): truth=" 
                     << truth << ", computed=" << computed;
                throw logic_error( msg.str() );
            }
        }
    }
    for( int i=0; i<m; ++i )
    {
        int truth = pRef.LocalEntry(i,0);
        int computed = p_copy.LocalEntry(i,0);
        if( truth != computed+1 /* 0 vs. 1 indexing */ )
        {
            ostringstream msg;
            msg << "Pivots off at index " << i << ": truth=" << truth 
                 << ", computed=" << computed;
            throw logic_error( msg.str() );
        }
    }
    Barrier( g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestLU
( bool testCorrectness, bool printMatrices,
  int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g);
    DistMatrix<T,Star,Star> ARef(g);
    DistMatrix<int,VC,Star> p(g);

    A.ResizeTo( m, m );
    p.ResizeTo( m, 1 );

    A.SetToRandom();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        ARef = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting LU factorization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::LU( A, p );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::LUGFlops<T>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        A.Print("A after factorization");
        p.Print("p after factorization");
    }
    if( testCorrectness )
    {
        TestCorrectness( printMatrices, A, p, ARef );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 7 )
    {
        if( rank == 0 )
            Usage();
        elemental::Finalize();
        return 0;
    }
    try
    {
        const int   r = atoi( argv[1] );
        const int   c = atoi( argv[2] );
        const int   m = atoi( argv[3] );
        const int   nb = atoi( argv[4] );
        const bool  testCorrectness = atoi( argv[5] );
        const bool  printMatrices = atoi( argv[6] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        Grid g( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
            cout << "Will test LU" << endl;

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestLU<double>
        ( testCorrectness, printMatrices, m, g );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestLU<dcomplex>
        ( testCorrectness, printMatrices, m, g );
        if( rank == 0 )
            cout << endl;
#endif
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:" << endl 
             << e.what() << endl;
    }   
    elemental::Finalize();
    return 0;
}

