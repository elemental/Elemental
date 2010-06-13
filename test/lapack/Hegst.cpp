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
    cout << "Reduced a Hermitian GEneralized EVP to Hermitian STandard EVP" <<
    endl << endl;
    cout << "  Hegst <r> <c> <bothOnLeft> <shape> <m> <nb> <test correctness?>"
         << " <print matrices?>" << endl << endl;
    cout << "  r: number of process rows                   " << endl;
    cout << "  c: number of process cols                   " << endl;
    cout << "  bothOnLeft: we solve A X = B X Lambda iff 0 " << endl;
    cout << "  shape: {L,U}                                " << endl;
    cout << "  m: height of matrix                         " << endl;
    cout << "  nb: algorithmic blocksize                   " << endl;
    cout << "  test correctness?: false iff 0              " << endl;
    cout << "  print matrices?: false iff 0                " << endl;
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
  bool bothOnLeft, 
  Shape shape,  
        DistMatrix<T,Star,Star>& ARef,
  const DistMatrix<T,Star,Star>& BRef )
{
    const Grid& g = A.GetGrid();
    const int m = ARef.Height();
    DistMatrix<T,Star,Star> ACopy(g);

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    ACopy = A;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    lapack::internal::LocalHegst( bothOnLeft, shape, ARef, BRef );
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        ARef.Print("Truth");

    if( g.VCRank() == 0 )
    {
        cout << "  Testing correctness..."; 
        cout.flush();
    }
    if( shape == Lower )
    {
        for( int j=0; j<m; ++j )
        {
            for( int i=j; i<m; ++i )
            {
                T truth = ARef.LocalEntry(i,j);
                T computed = ACopy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j << "): truth="
                         << truth << ", computed=" << computed;
                    throw logic_error( msg.str() );
                }
            }
        }
    }
    else
    {
        for( int j=0; j<m; ++j )
        {
            for( int i=0; i<=j; ++i )
            {
                T truth = ARef.LocalEntry(i,j);
                T computed = ACopy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j << "): truth="
                         << truth << ", computed=" << computed;
                    throw logic_error( msg.str() );
                }
            }
        }
    }
    Barrier( g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestHegst
( bool testCorrectness, bool printMatrices,
  bool bothOnLeft, Shape shape, int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g);
    DistMatrix<T,MC,MR> B(g);
    DistMatrix<T,Star,Star> ARef(g);
    DistMatrix<T,Star,Star> BRef(g);

    A.ResizeTo( m, m );
    B.ResizeTo( m, m );

    A.SetToRandomHPD();
    B.SetToRandomHPD();
    B.MakeTrapezoidal( Left, shape );
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrices...";
            cout.flush();
        }
        ARef = A;
        BRef = B;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting reduction to Hermitian standard EVP...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Hegst( bothOnLeft, shape, A, B );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::HegstGFlops<T>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = "
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after reduction");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, bothOnLeft, shape, ARef, BRef );
}

int main( int argc, char* argv[] )
{
    int rank;
    elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 9 )
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
        const bool  bothOnLeft = ( atoi(argv[3]) != 0 );
        const Shape shape = CharToShape( *argv[4] );
        const int   m = atoi( argv[5] );
        const int   nb = atoi( argv[6] );
        const bool  testCorrectness = ( atoi(argv[7]) != 0 );
        const bool  printMatrices = ( atoi(argv[8]) != 0 );
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
        {
            cout << "Will test Hegst" << bothOnLeft << ShapeToChar(shape)
                 << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestHegst<double>
        ( testCorrectness, printMatrices, bothOnLeft, shape, m, g );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestHegst<dcomplex>
        ( testCorrectness, printMatrices, bothOnLeft, shape, m, g );
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

