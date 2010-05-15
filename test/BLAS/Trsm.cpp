/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include <cmath>
#include <ctime>
#include <sstream>
#include "Elemental.hpp"
#include "Elemental/BLASInternal.hpp"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

void Usage()
{
    cout << "TRiangular Solve with Multiple right-hand sides." << endl << endl;
    cout << "  Trsm <r> <c> <side> <shape> <orientation> <unit diag?>" << endl
         << "       <m> <n> <nb> <test correctness?> <print matrices?>" << endl;
    cout << endl;
    cout << "  r: number of process rows                  " << endl;
    cout << "  c: number of process cols                  " << endl;
    cout << "  side: {L,R}                                " << endl;
    cout << "  shape: {L,U}                               " << endl;
    cout << "  orientation: {N,T,C}                       " << endl;
    cout << "  diag?: {N,U}                               " << endl;
    cout << "  m: height of right-hand sides              " << endl;
    cout << "  n: number of right-hand sides              " << endl;
    cout << "  nb: algorithmic blocksize                  " << endl;
    cout << "  test correctness?: false iff 0             " << endl;
    cout << "  print matrices?: false iff 0               " << endl;
    cout << endl;
}

template<typename T>
bool OKRelativeError( T truth, T computed );

template<>
bool OKRelativeError( double truth, double computed )
{ return ( fabs(truth-computed) / max(fabs(truth),(double)1) <= 1e-12 ); }

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{ return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-12 ); }
#endif

template<typename T>
void TestCorrectness
( bool printMatrices,
  const DistMatrix<T,MC,MR>& X,
  Side side, Shape shape,
  Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,Star,Star>& ARef,
                 DistMatrix<T,Star,Star>& XRef )
{
    const Grid& grid = X.GetGrid();
    DistMatrix<T,Star,Star> X_copy(grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Copying computed result...";
        cout.flush();
    }
    X_copy = X;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    BLAS::Trsm( side, shape, orientation, diagonal,
                alpha, ARef.LockedLocalMatrix(),
                       XRef.LocalMatrix()         );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        XRef.Print("Truth");

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<X.Width(); ++j )
    {
        for( int i=0; i<X.Height(); ++i )
        {
            T truth = XRef.LocalEntry(i,j);
            T computed = X_copy.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                ostringstream msg;
                msg << "FAILED at index (" << i << "," << j << "): truth="
                     << truth << ", computed=" << computed;
                const string s = msg.str();
                throw s.c_str();
            }
        }
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestTrsm
( bool testCorrectness, bool printMatrices,
  Side side, Shape shape,
  Orientation orientation, Diagonal diagonal,
  int m, int n, T alpha, const Grid& grid )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,MC,MR> X(grid);
    DistMatrix<T,Star,Star> ARef(grid);
    DistMatrix<T,Star,Star> XRef(grid);

    if( side == Left )
        A.ResizeTo( m, m );
    else
        A.ResizeTo( n, n );
    X.ResizeTo( m, n );

    A.SetToRandomDiagDominant();
    X.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        ARef = A;
        XRef = X;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        X.Print("X");
    }
    if( grid.VCRank() == 0 )
    {
        cout << "  Starting Trsm...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    BLAS::Trsm( side, shape, orientation, diagonal, alpha, A, X );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = BLAS::Internal::TrsmGFlops<T>(side,m,n,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops =  " << gFlops << endl;
    if( printMatrices )
        X.Print("X after solve");
    if( testCorrectness )
    {
        TestCorrectness
        ( printMatrices, X,
          side, shape, orientation, diagonal, alpha, ARef, XRef );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 12 )
    {
        if( rank == 0 )
            Usage();
        Elemental::Finalize();
        return 0;
    }
    try
    {
        const int         r = atoi( argv[1] );
        const int         c = atoi( argv[2] );
        const Side        side = CharToSide( *argv[3] );
        const Shape       shape = CharToShape( *argv[4] );
        const Orientation orientation = CharToOrientation( *argv[5] );
        const Diagonal    diagonal = CharToDiagonal( *argv[6] );
        const int         m = atoi( argv[7] );
        const int         n = atoi( argv[8] );
        const int         nb = atoi( argv[9] );
        const bool        testCorrectness = atoi( argv[10] );
        const bool        printMatrices = atoi( argv[11] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        Grid grid( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
        {
            cout << "Will test Trsm" << SideToChar(side) 
                                     << ShapeToChar(shape)
                                     << OrientationToChar(orientation) 
                                     << DiagonalToChar(diagonal) << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestTrsm<double>
        ( testCorrectness, printMatrices,
          side, shape, orientation, diagonal, m, n, (double)3, grid );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestTrsm<dcomplex>
        ( testCorrectness, printMatrices,
          side, shape, orientation, diagonal, m, n, (dcomplex)3, grid );
        if( rank == 0 )
            cout << endl;
#endif
    }
    catch( const char* errorMsg )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:" << endl 
             << errorMsg << endl;
    }   
    Elemental::Finalize();
    return 0;
}

