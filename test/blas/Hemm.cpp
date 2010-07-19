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
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "HErmitian Matrix Matrix multiplication." << endl << endl;;
    cout << "  Hemm <r> <c> <Side> <Shape> <m> <n> <nb> <correctness?> <print?>"
         << endl << endl;
    cout << "  r: number of process rows    " << endl;
    cout << "  c: number of process cols    " << endl;
    cout << "  Side: {L,R}                  " << endl;
    cout << "  Shape: {L,U}                 " << endl;
    cout << "  m: height of C               " << endl;
    cout << "  n: width  of C               " << endl;
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
  const DistMatrix<T,MC,MR>& C,
  Side side, Shape shape,
  T alpha, const DistMatrix<T,Star,Star>& A_ref,
           const DistMatrix<T,Star,Star>& B_ref,
  T beta,        DistMatrix<T,Star,Star>& C_ref )
{
    const Grid& g = C.GetGrid();
    DistMatrix<T,Star,Star> CCopy(g);

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    CCopy = C;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    blas::Hemm
    ( side, shape,
      alpha, A_ref.LockedLocalMatrix(),
             B_ref.LockedLocalMatrix(),
      beta,  C_ref.LocalMatrix() );
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        C_ref.Print("Truth");

    if( g.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<C.Width(); ++j )
    {
        for( int i=0; i<C.Height(); ++i )
        {
            T truth = C_ref.LocalEntry(i,j);
            T computed = CCopy.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                ostringstream msg;
                msg << "FAILED at index (" << i << "," << j << "): truth=" 
                    << truth << ", computed=" << computed;
                throw logic_error( msg.str() );
            }
        }
    }
    Barrier( g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestHemm
( bool testCorrectness, bool printMatrices,
  Side side, Shape shape,
  int m, int n, T alpha, T beta, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,  MR  > A(g);
    DistMatrix<T,MC,  MR  > B(g);
    DistMatrix<T,MC,  MR  > C(g);
    DistMatrix<T,Star,Star> A_ref(g);
    DistMatrix<T,Star,Star> B_ref(g);
    DistMatrix<T,Star,Star> C_ref(g);

    if( side == Left )
        A.ResizeTo( m, m );
    else
        A.ResizeTo( n, n );
    B.ResizeTo( m, n );
    C.ResizeTo( m, n );

    // Test Hemm
    if( g.VCRank() == 0 )
        cout << "Hemm:" << endl;
    A.SetToRandomHPD();
    B.SetToRandom();
    C.SetToRandom();

    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        A_ref = A;
        B_ref = B;
        C_ref = C;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.VCRank() == 0 )
    {
        cout << "  Starting Parallel Hemm...";
        cout.flush();
    }
    Barrier( g.VCComm() );
    startTime = Time();
    blas::Hemm( side, shape, alpha, A, B, beta, C );
    Barrier( g.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::HemmGFlops<T>(side,m,n,runTime);
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        if( side == Left )
            msg << "C := " << alpha << " Herm(A) B + " << beta << " C";
        else
            msg << "C := " << alpha << " B Herm(A) + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestCorrectness
        ( printMatrices, C,
          side, shape, alpha, A_ref, B_ref, beta, C_ref );
    }
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
        const Side  side = CharToSide( *argv[3] );
        const Shape shape = CharToShape( *argv[4] );
        const int   m = atoi( argv[5] );
        const int   n = atoi( argv[6] );
        const int   nb = atoi( argv[7] );
        const bool  testCorrectness = atoi( argv[8] );
        const bool  printMatrices = atoi( argv[9] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        const Grid g( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
        {
            cout << "Will test Hemm" << SideToChar(side) 
                                     << ShapeToChar(shape) << endl;
        }

        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with doubles:                 " << endl;
            cout << "--------------------------------------" << endl;
        }
        TestHemm<double>
        ( testCorrectness, printMatrices,
          side, shape, m, n, (double)3, (double)4, g );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestHemm<dcomplex>
        ( testCorrectness, printMatrices,
          side, shape, m, n, (dcomplex)3, (dcomplex)4, g );
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

