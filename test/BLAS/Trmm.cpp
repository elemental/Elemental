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
#include "Elemental.h"
#include "ElementalBLASInternal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

void Usage()
{
    cout << "TRiangular Matrix-Matrix multiplication." << endl << endl;
    cout << "  Trmm <r> <c> <side> <shape> <orientation> <unit diag?>" << endl
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
bool OKRelativeError( float truth, float computed )
{
    return ( fabs(truth-computed) / max(fabs(truth),(float)1) <= 1e-3 );
}

template<>
bool OKRelativeError( double truth, double computed )
{
    return ( fabs(truth-computed) / max(fabs(truth),(double)1) <= 1e-13 );
}

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( scomplex truth, scomplex computed )
{
    return ( norm(truth-computed) / max(norm(truth),(float)1) <= 1e-3 );
}

template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{
    return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-13 );
}
#endif

template<typename T>
void TestCorrectness( const Side side,               const Shape shape,
                      const Orientation orientation, const Diagonal diagonal,
                      const T alpha, 
                      const DistMatrix<T,Star,Star>& A_ref,
                            DistMatrix<T,Star,Star>& X_ref,
                      const DistMatrix<T,MC,MR>& X,
                      const bool printMatrices                               )
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
    BLAS::Trmm( side, shape, orientation, diagonal,
                alpha, A_ref.LockedLocalMatrix(),
                       X_ref.LocalMatrix()         );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        X_ref.Print("Truth");

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<X.Width(); ++j )
    {
        for( int i=0; i<X.Height(); ++i )
        {
            T truth = X_ref.LocalEntry(i,j);
            T computed = X_copy.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                cout << "FAILED on process " << grid.VCRank() 
                     << " at index (" << i << "," << j << "): truth="
                     << truth << ", computed=" << computed << endl;
                throw exception();
            }
        }
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestTrmm
( const Side side,               const Shape shape,
  const Orientation orientation, const Diagonal diagonal,
  const int m, const int n, const T alpha,
  const bool testCorrectness, const bool printMatrices, const Grid& grid   )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,MC,MR> X(grid);
    DistMatrix<T,Star,Star> A_ref(grid);
    DistMatrix<T,Star,Star> X_ref(grid);

    if( side == Left )
        A.ResizeTo( m, m );
    else
        A.ResizeTo( n, n );
    X.ResizeTo( m, n );

    A.SetToRandom();
    X.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        A_ref = A;
        X_ref = X;
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
        cout << "  Starting Trmm...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    BLAS::Trmm( side, shape, orientation, diagonal, alpha, A, X );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = BLAS::Internal::TrmmGFlops<T>(side,m,n,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops =  " << gFlops << endl;
    if( printMatrices )
        X.Print("X after solve");
    if( testCorrectness )
    {
        TestCorrectness( side, shape, orientation, diagonal, 
                         alpha, A_ref, X_ref, X, printMatrices );
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
            cout << "Will test Trmm" << SideToChar(side) 
                                     << ShapeToChar(shape)
                                     << OrientationToChar(orientation) 
                                     << DiagonalToChar(diagonal) << endl;
        }

        if( (diagonal == NonUnit && m<=2000) || (diagonal == Unit && m <= 100) )
        {
            if( rank == 0 )
            {
                cout << "--------------------" << endl;
                cout << "Testing with floats:" << endl;
                cout << "--------------------" << endl;
            }
            TestTrmm<float>
            ( side, shape, orientation, diagonal,
              m, n, (float)3, testCorrectness, printMatrices, grid );
            if( rank == 0 )
                cout << endl;
        }
        else
        {
            if( rank == 0 )
            {
                cout << "--------------------------------" << endl;
                cout << "Floats unsuitable for this test." << endl;
                cout << "--------------------------------" << endl;
                cout << endl;
            }
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestTrmm<double>
        ( side, shape, orientation, diagonal,
          m, n, (double)3, testCorrectness, printMatrices, grid );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with single-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestTrmm<scomplex>
        ( side, shape, orientation, diagonal,
          m, n, (scomplex)3, testCorrectness, printMatrices, grid );
        if( rank == 0 )
            cout << endl;

        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestTrmm<dcomplex>
        ( side, shape, orientation, diagonal,
          m, n, (dcomplex)3, testCorrectness, printMatrices, grid );
        if( rank == 0 )
            cout << endl;
#endif
    }
    catch( exception e )
    {
        cerr << "Caught exception on process " << rank << endl;
    }
    Elemental::Finalize();
    return 0;
}

