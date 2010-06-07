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
    const Grid& g = C.GetGrid();
    DistMatrix<T,Star,Star> C_copy(g);

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    C_copy = C;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    blas::Syr2k
    ( shape, orientation,
      alpha, ARef.LockedLocalMatrix(),
             BRef.LockedLocalMatrix(),
      beta,  CRef.LocalMatrix()       );
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        CRef.Print("Truth:");

    if( g.VCRank() == 0 )
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
                    throw logic_error( msg.str() );
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
void TestSyr2k
( bool testCorrectness, bool printMatrices,
  Shape shape, Orientation orientation,
  int m, int k, T alpha, T beta, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g);
    DistMatrix<T,MC,MR> B(g);
    DistMatrix<T,MC,MR> C(g);
    DistMatrix<T,Star,Star> ARef(g);
    DistMatrix<T,Star,Star> BRef(g);
    DistMatrix<T,Star,Star> CRef(g);

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
        if( g.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        ARef = A;
        BRef = B;
        CRef = C;
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
        cout << "  Starting Syr2k...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    blas::Syr2k( shape, orientation, alpha, A, B, beta, C );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::Syr2kGFlops<T>(m,k,runTime);
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
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
        Grid g( MPI_COMM_WORLD, r, c );
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
          shape, orientation, m, k, (double)3, (double)4, g );
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
          shape, orientation, m, k, (double)3, (double)4, g );
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

