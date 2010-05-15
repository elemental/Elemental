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
    cout << "SYmmetric Rank-2K update." << endl << endl;;
    cout << "  Syr2k <r> <c> <shape> <trans?> <m> <k> <nb> ";
    cout << "<correctness?> <print?>   " << endl << endl;;
    cout << "  r: number of process rows             " << endl;
    cout << "  c: number of process cols             " << endl;
    cout << "  shape?: {L,U}                         " << endl;
    cout << "  trans?: {N,T}                         " << endl;
    cout << "  m: height of C                        " << endl;
    cout << "  k: inner dimension                    " << endl;
    cout << "  nb: algorithmic blocksize             " << endl;
    cout << "  correctness?: false iff 0             " << endl;
    cout << "  print?: false iff 0                   " << endl;
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
  const DistMatrix<T,MC,MR>& C,
  Shape shape, Orientation orientation,
  T alpha, const DistMatrix<T,Star,Star>& ARef,
           const DistMatrix<T,Star,Star>& BRef,
  T beta,        DistMatrix<T,Star,Star>& CRef )
{
    const Grid& grid = C.GetGrid();
    DistMatrix<T,Star,Star> C_copy(grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    C_copy = C;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    BLAS::Syr2k
    ( shape, orientation,
      alpha, ARef.LockedLocalMatrix(),
             BRef.LockedLocalMatrix(),
      beta,  CRef.LocalMatrix()       );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        CRef.Print("Truth:");

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    if( shape == Lower )
    {
        for( int j=0; j<C.Width(); ++j )
        {
            for( int i=j; i<C.Height(); ++i )
            {
                T truth = CRef.LocalEntry(i,j);
                T computed = C_copy.LocalEntry(i,j);

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
    }
    else
    {
        for( int j=0; j<C.Width(); ++j )
        {
            for( int i=0; i<=j; ++i )
            {
                T truth = CRef.LocalEntry(i,j);
                T computed = C_copy.LocalEntry(i,j);

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
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestSyr2k
( bool testCorrectness, bool printMatrices,
  Shape shape, Orientation orientation,
  int m, int k, T alpha, T beta, const Grid& grid )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,MC,MR> B(grid);
    DistMatrix<T,MC,MR> C(grid);
    DistMatrix<T,Star,Star> ARef(grid);
    DistMatrix<T,Star,Star> BRef(grid);
    DistMatrix<T,Star,Star> CRef(grid);

    if( orientation == Normal )
    {
        A.ResizeTo( m, k );
        B.ResizeTo( m, k );
    }
    else
    {
        A.ResizeTo( k, m );
        B.ResizeTo( k, m );
    }

    C.ResizeTo( m, m );

    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
    C.MakeTrapezoidal( Left, shape );
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        ARef = A;
        BRef = B;
        CRef = C;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( grid.VCRank() == 0 )
    {
        cout << "  Starting Syr2k...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    BLAS::Syr2k( shape, orientation, alpha, A, B, beta, C );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = BLAS::Internal::Syr2kGFlops<T>(m,k,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        ostringstream msg;
        if( orientation == Normal )
            msg << "C := " << alpha << " A B' + B A'" << beta << " C";
        else
            msg << "C := " << alpha << " A' B + B' A" << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestCorrectness
        ( printMatrices, C,
          shape, orientation, alpha, ARef, BRef, beta, CRef );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 10 )
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
        const Shape       shape = CharToShape( *argv[3] );
        const Orientation orientation = CharToOrientation( *argv[4] );
        const int         m = atoi( argv[5] );
        const int         k = atoi( argv[6] );
        const int         nb = atoi( argv[7] );
        const bool        testCorrectness = atoi( argv[8] );
        const bool        printMatrices = atoi( argv[9] );
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
            cout << "Will test Syr2k" << ShapeToChar(shape) 
                                      << OrientationToChar(orientation) << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestSyr2k<double>
        ( testCorrectness, printMatrices,
          shape, orientation, m, k, (double)3, (double)4, grid );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestSyr2k<dcomplex>
        ( testCorrectness, printMatrices,
          shape, orientation, m, k, (double)3, (double)4, grid );
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

