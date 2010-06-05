/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include <ctime>
#include "elemental.hpp"
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "SYmmetric Matrix Matrix multiplication." << endl << endl;;
    cout << "  Symm <r> <c> <Side> <Shape> <m> <n> <nb> <correctness?> <print?>"
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
{ return ( fabs(truth-computed) / max(fabs(truth),(double)1)  <= 1e-13 ); }

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{ return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-13 ); }
#endif

template<typename T>
void TestCorrectness
( bool printMatrices,
  const DistMatrix<T,MC,MR>& C,
  Side side, Shape shape,
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
    blas::Symm( side, shape,
                alpha, ARef.LockedLocalMatrix(),
                       BRef.LockedLocalMatrix(),
                beta,  CRef.LocalMatrix()       );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        CRef.Print("Truth");

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<C.Width(); ++j )
    {
        for( int i=0; i<C.Height(); ++i )
        {
            T truth = CRef.LocalEntry(i,j);
            T computed = C_copy.LocalEntry(i,j);

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
void TestSymm
( const Side side, const Shape shape,
  const int m, const int n, const T alpha, const T beta,
  const bool testCorrectness, const bool printMatrices, const Grid& grid  )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,MC,MR> B(grid);
    DistMatrix<T,MC,MR> C(grid);
    DistMatrix<T,Star,Star> ARef(grid);
    DistMatrix<T,Star,Star> BRef(grid);
    DistMatrix<T,Star,Star> CRef(grid);

    if( side == Left )
        A.ResizeTo( m, m );
    else
        A.ResizeTo( n, n );
    B.ResizeTo( m, n );
    C.ResizeTo( m, n );

    // Test Symm
    if( grid.VCRank() == 0 )
        cout << "Symm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
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
        cout << "  Starting Parallel Symm...";
        cout.flush();
    }
    Barrier( grid.VCComm() );
    startTime = Time();
    blas::Symm
    ( side, shape, alpha, A, B, beta, C );
    Barrier( grid.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::SymmGFlops<T>(side,m,n,runTime);
    if( grid.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        if( side == Left )
            msg << "C := " << alpha << " Symm(A) B + " << beta << " C";
        else
            msg << "C := " << alpha << " B Symm(A) + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestCorrectness
        ( printMatrices, C,
          side, shape, alpha, ARef, BRef, beta, CRef );
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
        const Grid grid( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
        {
            cout << "Will test Symm" << SideToChar(side) 
                                     << ShapeToChar(shape) << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestSymm<double>
        ( side, shape, m, n, (double)3, (double)4, 
          testCorrectness, printMatrices, grid    ); 
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestSymm<dcomplex>
        ( side, shape, m, n, (dcomplex)3, (dcomplex)4, 
          testCorrectness, printMatrices, grid        ); 
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

