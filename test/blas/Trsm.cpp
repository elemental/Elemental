/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include <ctime>
#include "elemental.hpp"
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "TRiangular Solve with Multiple right-hand sides.\n\n"
         << "  Trsm <r> <c> <side> <shape> <orientation> <unit diag?> <m> <n> "
            "<nb> <correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  side: {L,R}\n"
         << "  shape: {L,U}\n"
         << "  orientation: {N,T,C}\n"
         << "  diag?: {N,U}\n"
         << "  m: height of right-hand sides\n"
         << "  n: number of right-hand sides\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
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
    const Grid& g = X.Grid();
    DistMatrix<T,Star,Star> XCopy(g);

    if( g.VCRank() == 0 )
    {
        cout << "  Copying computed result...";
        cout.flush();
    }
    XCopy = X;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    blas::Trsm( side, shape, orientation, diagonal,
                alpha, ARef.LockedLocalMatrix(),
                       XRef.LocalMatrix() );
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        XRef.Print("Truth");

    if( g.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<X.Width(); ++j )
    {
        for( int i=0; i<X.Height(); ++i )
        {
            T truth = XRef.GetLocalEntry(i,j);
            T computed = XCopy.GetLocalEntry(i,j);

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
void TestTrsm
( bool testCorrectness, bool printMatrices,
  Side side, Shape shape,
  Orientation orientation, Diagonal diagonal,
  int m, int n, T alpha, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g);
    DistMatrix<T,MC,MR> X(g);
    DistMatrix<T,Star,Star> ARef(g);
    DistMatrix<T,Star,Star> XRef(g);

    if( side == Left )
        A.ResizeTo( m, m );
    else
        A.ResizeTo( n, n );
    X.ResizeTo( m, n );

    A.SetToRandomHPD();
    X.SetToRandom();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        ARef = A;
        XRef = X;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        X.Print("X");
    }
    if( g.VCRank() == 0 )
    {
        cout << "  Starting Trsm...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    blas::Trsm( side, shape, orientation, diagonal, alpha, A, X );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::TrsmGFlops<T>(side,m,n,runTime);
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
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
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 12 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    try
    {
        const int r = atoi(argv[1]);
        const int c = atoi(argv[2]);
        const Side side = CharToSide(*argv[3]);
        const Shape shape = CharToShape(*argv[4]);
        const Orientation orientation = CharToOrientation(*argv[5]);
        const Diagonal diagonal = CharToDiagonal(*argv[6]);
        const int m = atoi(argv[7]);
        const int n = atoi(argv[8]);
        const int nb = atoi(argv[9]);
        const bool testCorrectness = atoi(argv[10]);
        const bool printMatrices = atoi(argv[11]);
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        const Grid g( MPI_COMM_WORLD, r, c );
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
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTrsm<double>
        ( testCorrectness, printMatrices,
          side, shape, orientation, diagonal, m, n, (double)3, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestTrsm<dcomplex>
        ( testCorrectness, printMatrices,
          side, shape, orientation, diagonal, m, n, (dcomplex)3, g );
#endif
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:\n"
             << e.what() << endl;
    }   
    Finalize();
    return 0;
}

