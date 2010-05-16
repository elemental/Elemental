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
#include "elemental.hpp"
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "SYmmetric Matrix vector multiplication." << endl << endl;;
    cout << "  Symv <r> <c> <Shape> <m> <nb> <correctness?> <print?>"
         << endl << endl;
    cout << "  r: number of process rows    " << endl;
    cout << "  c: number of process cols    " << endl;
    cout << "  Shape: {L,U}                 " << endl;
    cout << "  m: height of C               " << endl;
    cout << "  nb: algorithmic blocksize    " << endl;
    cout << "  correctness?: [0/1]          " << endl;
    cout << "  print?: [0/1]                " << endl;
    cout << endl;
}

template<typename T>
bool OKRelativeError( T truth, T computed );

template<>
bool OKRelativeError( double truth, double computed )
{ return ( fabs(truth-computed) / max(fabs(truth),(double)1)  <= 1e-12 ); }

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{ return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-12 ); }
#endif

template<typename T>
void TestCorrectness
( bool printMatrices,
  const DistMatrix<T,MC,MR>& y,
  Shape shape,
  T alpha, const DistMatrix<T,Star,Star>& ARef,
           const DistMatrix<T,Star,Star>& xRef,
  T beta,        DistMatrix<T,Star,Star>& yRef )
{
    const Grid& grid = y.GetGrid();
    DistMatrix<T,Star,Star> y_copy(grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    y_copy = y;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    blas::Symv( shape,
                alpha, ARef.LockedLocalMatrix(),
                       xRef.LockedLocalMatrix(),
                beta,  yRef.LocalMatrix()       );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        yRef.Print("Truth");

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<y.Width(); ++j )
    {
        for( int i=0; i<y.Height(); ++i )
        {
            T truth = yRef.LocalEntry(i,j);
            T computed = y_copy.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                ostringstream msg;
                msg << "FAILED at index (" << i << "," << j << "): truth="
                     << truth << ", computed=" << computed;
                const string& s = msg.str();
                throw s.c_str();
            }
        }
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestSymv
( const Shape shape,
  const int m, const T alpha, const T beta, 
  const bool testCorrectness, const bool printMatrices, const Grid& grid )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,MC,MR> x(grid);
    DistMatrix<T,MC,MR> y(grid);
    DistMatrix<T,Star,Star> ARef(grid);
    DistMatrix<T,Star,Star> xRef(grid);
    DistMatrix<T,Star,Star> yRef(grid);

    A.ResizeTo( m, m );
    x.ResizeTo( m, 1 );
    y.ResizeTo( m, 1 );

    // Test Symm
    if( grid.VCRank() == 0 )
        cout << "Symm:" << endl;
    A.SetToRandom();
    x.SetToRandom();
    y.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        ARef = A;
        xRef = x;
        yRef = y;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        x.Print("x");
        y.Print("y");
    }
    if( grid.VCRank() == 0 )
    {
        cout << "  Starting Parallel Symv...";
        cout.flush();
    }
    Barrier( grid.VCComm() );
    startTime = Time();
    blas::Symv
    ( shape, alpha, A, x, beta, y );
    Barrier( grid.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::SymvGFlops<T>(m,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        ostringstream msg;
        msg << "y := " << alpha << " Symm(A) x + " << beta << " y";
        y.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestCorrectness
        ( printMatrices, y, shape, alpha, ARef, xRef, beta, yRef );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 8 )
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
        const int   m = atoi( argv[4] );
        const int   nb = atoi( argv[5] );
        const bool  testCorrectness = atoi( argv[6] );
        const bool  printMatrices = atoi( argv[7] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        const Grid grid( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
            cout << "Will test Symv" << ShapeToChar(shape) << endl;

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestSymv<double>
        ( shape, m, (double)3, (double)4, testCorrectness, printMatrices, 
          grid );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestSymv<dcomplex>
        ( shape, m, (dcomplex)3, (dcomplex)4, testCorrectness, printMatrices,
          grid );
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
    elemental::Finalize();
    return 0;
}

