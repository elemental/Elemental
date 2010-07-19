/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <ctime>
#include "elemental.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Generates SPD matrix then solves for its Cholesky factor."
         << endl << endl;
    cout << "  Chol <r> <c> <shape> <var2/3> <naive> <m> <nb> "
            "<test correctness?> <print matrices?>" << endl << endl;
    cout << "  r: number of process rows      " << endl;
    cout << "  c: number of process cols      " << endl;
    cout << "  shape: {L,U}                   " << endl;
    cout << "  var2/3: 2 iff 0                " << endl;
    cout << "  naive: smart algorithms iff 0  " << endl;
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
  Shape shape, DistMatrix<T,Star,Star>& ARef )
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
    lapack::internal::LocalChol( shape, ARef );
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
void TestChol
( bool var3, bool naive,
  bool testCorrectness, bool printMatrices, 
  Shape shape, int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g);
    DistMatrix<T,Star,Star> ARef(g);

    A.ResizeTo( m, m );

    A.SetToRandomHPD();
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
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting Cholesky factorization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    if( var3 )
    {
        if( naive )
            lapack::internal::CholVar3Naive( shape, A );
        else
            lapack::internal::CholVar3( shape, A );
    }
    else
    {
        if( naive )
            lapack::internal::CholVar2Naive( shape, A );
        else
            lapack::internal::CholVar2( shape, A );
    }
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::CholGFlops<T>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, shape, ARef );
}

int main( int argc, char* argv[] )
{
    int rank;
    elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 10 )
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
        const Shape shape = CharToShape( *argv[3] );
        const bool  var3 = ( atoi(argv[4]) != 0 );
        const bool  naive = ( atoi(argv[5]) != 0 );
        const int   m = atoi( argv[6] );
        const int   nb = atoi( argv[7] );
        const bool  testCorrectness = ( atoi(argv[8]) != 0 );
        const bool  printMatrices = ( atoi(argv[9]) != 0 );
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
            cout << "Will test Chol" << ShapeToChar(shape) << ", Var"
                 << ( var3 ? "3" : "2" ) << ( naive ? "Naive" : "" ) << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestChol<double>
        ( var3, naive, testCorrectness, printMatrices, shape, m, g );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestChol<dcomplex>
        ( var3, naive, testCorrectness, printMatrices, shape, m, g );
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

