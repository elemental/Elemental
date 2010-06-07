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
    cout << "GEneral Matrix Matrix multiplication." << endl << endl;;
    cout << "  Gemm <r> <c> <orient. of A?> <orient. of B?> <m> <n> <k> <nb> "  
         << endl
         << "       <serial gemm?> <parallel gemm?> <correctness?> <print?> " 
         << endl << endl;
    cout << "  r: number of process rows    " << endl;
    cout << "  c: number of process cols    " << endl;
    cout << "  orient. of A: {N,T,C}        " << endl;
    cout << "  orient. of B: {N,T,C}        " << endl;
    cout << "  m: height of C               " << endl;
    cout << "  n: width  of C               " << endl;
    cout << "  k: inner dimension of AB     " << endl;
    cout << "  nb: algorithmic blocksize    " << endl;
    cout << "  serial gemm? [0/1]           " << endl;
    cout << "  parallel gemm? [0/1]         " << endl;
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
void ManualGemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta, Matrix<T>& C )
{
    const int m = C.Height();
    const int n = C.Width();
    blas::Scal( beta, C );
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        const int k = A.Width();
        for( int j=0; j<n; ++j )
            for( int i=0; i<m; ++i )
                for( int l=0; l<k; ++l )
                    C(i,j) += alpha * A(i,l) * B(l,j);
    }
    else if( orientationOfA == Normal )
    {
        const int k = A.Width();
        if( orientationOfB == Transpose )
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * A(i,l) * B(j,l);
        }
        else
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * A(i,l) * Conj(B(j,l));
        }
    }
    else if( orientationOfB == Normal )
    {
        const int k = A.Height();
        if( orientationOfA == Transpose )
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * A(l,i) * B(l,j);
        }
        else
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * Conj(A(l,i)) * B(l,j);
        }
    }
    else
    {
        const int k = A.Height();
        if( orientationOfA == Transpose && orientationOfB == Transpose )
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * A(l,i) * B(j,l);
        }
        else if( orientationOfA == Transpose )
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * A(l,i) * Conj(B(j,l));
        }
        else if( orientationOfB == Transpose )
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * Conj(A(l,i)) * B(j,l);
        }
        else
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C(i,j) += alpha * Conj(A(l,i)) * Conj(B(j,l));
        }
    }
}

template<typename T>
void TestSerialCorrectness
( bool printMatrices,
  const DistMatrix<T,Star,Star>& C,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,Star,Star>& ARef,
           const DistMatrix<T,Star,Star>& BRef,
  T beta,        DistMatrix<T,Star,Star>& CRef )
{
    const Grid& g = C.GetGrid();

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    ManualGemm( orientationOfA, orientationOfB,
                alpha, ARef.LockedLocalMatrix(),
                       BRef.LockedLocalMatrix(),
                beta,  CRef.LocalMatrix()       );
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        CRef.Print("Truth");

    if( g.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<C.Width(); ++j )
    {
        for( int i=0; i<C.Height(); ++i )
        {
            T truth = CRef.LocalEntry(i,j);
            T computed = C.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                ostringstream msg;
                msg << "FAILED at index (" << i << "," << j << "): truth="
                    << truth << ", computed=" << computed << ", error=" 
                    << Abs(truth-computed) << endl;
                throw logic_error( msg.str() );
            }
        }
    }
    Barrier( g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestParallelCorrectness
( bool printMatrices,
  const DistMatrix<T,MC,MR>& C,
  Orientation orientationOfA, Orientation orientationOfB,
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
    blas::internal::LocalGemm
    ( orientationOfA, orientationOfB, alpha, ARef, BRef, beta, CRef );
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        CRef.Print("Truth");

    if( g.VCRank() == 0 )
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
                    << truth << ", computed=" << computed << ", error="
                    << Abs(truth-computed) << endl;
                throw logic_error( msg.str() );
            }
        }
    }
    Barrier( g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestSerialGemm
( bool testCorrectness, bool printMatrices,
  Orientation orientationOfA, Orientation orientationOfB,
  int m, int n, int k, T alpha, T beta, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,Star,Star> A(g);
    DistMatrix<T,Star,Star> B(g);
    DistMatrix<T,Star,Star> C(g);
    DistMatrix<T,Star,Star> ARef(g);
    DistMatrix<T,Star,Star> BRef(g);
    DistMatrix<T,Star,Star> CRef(g);

    if( orientationOfA == Normal )
        A.ResizeTo( m, k );
    else
        A.ResizeTo( k, m );

    if( orientationOfB == Normal )
        B.ResizeTo( k, n );
    else
        B.ResizeTo( n, k );

    C.ResizeTo( m, n );

    if( g.VCRank() == 0 )
        cout << "Serial Gemm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
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
        cout << "  Starting Gemm...";
        cout.flush();
    }
    Barrier( g.VCComm() );
    startTime = Time();
    blas::Gemm
    ( orientationOfA, orientationOfB, 
      alpha, A.LockedLocalMatrix(), B.LockedLocalMatrix(), 
      beta,  C.LocalMatrix()                              );
    Barrier( g.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::GemmGFlops<T>(m,n,k,runTime);
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl 
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestSerialCorrectness
        ( printMatrices, C,
          orientationOfA, orientationOfB, 
          alpha, ARef, BRef, beta, CRef );
    }
}

template<typename T>
void TestParallelGemm
( bool testCorrectness, bool printMatrices,
  Orientation orientationOfA, Orientation orientationOfB,
  int m, int n, int k, T alpha, T beta, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,  MR  > A(g);
    DistMatrix<T,MC,  MR  > B(g);
    DistMatrix<T,MC,  MR  > C(g);
    DistMatrix<T,Star,Star> ARef(g);
    DistMatrix<T,Star,Star> BRef(g);
    DistMatrix<T,Star,Star> CRef(g);

    if( orientationOfA == Normal )
        A.ResizeTo( m, k );
    else
        A.ResizeTo( k, m );

    if( orientationOfB == Normal )
        B.ResizeTo( k, n );
    else
        B.ResizeTo( n, k );

    C.ResizeTo( m, n );

    // Test the variant of Gemm that keeps A stationary
    if( g.VCRank() == 0 )
        cout << "Stationary A Algorithm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
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
        cout << "  Starting Parallel Gemm...";
        cout.flush();
    }
    Barrier( g.VCComm() );
    startTime = Time();
    blas::internal::GemmA
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    Barrier( g.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::GemmGFlops<T>(m,n,k,runTime);
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestParallelCorrectness
        ( printMatrices, C,
          orientationOfA, orientationOfB, 
          alpha, ARef, BRef, beta, CRef );
    }

    // Test the variant of Gemm that keeps B stationary
    if( g.VCRank() == 0 )
        cout << endl << "Stationary B Algorithm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
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
        cout << "  Starting Parallel Gemm...";
        cout.flush();
    }
    Barrier( g.VCComm() );
    startTime = Time();
    blas::internal::GemmB
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    Barrier( g.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::GemmGFlops<T>(m,n,k,runTime);
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl 
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestParallelCorrectness
        ( printMatrices, C,
          orientationOfA, orientationOfB,
          alpha, ARef, BRef, beta, CRef );
    }

    // Test the variant of Gemm that keeps C stationary
    if( g.VCRank() == 0 )
        cout << endl << "Stationary C Algorithm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
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
        cout << "  Starting Parallel Gemm...";
        cout.flush();
    }
    Barrier( g.VCComm() );
    startTime = Time();
    blas::internal::GemmC
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    Barrier( g.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = blas::internal::GemmGFlops<T>(m,n,k,runTime);
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestParallelCorrectness
        ( printMatrices, C, 
          orientationOfA, orientationOfB,
          alpha, ARef, BRef, beta, CRef );
    }
    
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        // Test the variant of Gemm for panel-panel dot products
        if( g.VCRank() == 0 )
            cout << endl << "Dot Product Algorithm:" << endl;
        A.SetToRandom();
        B.SetToRandom();
        C.SetToRandom();
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
            cout << "  Starting Parallel Gemm...";
            cout.flush();
        }
        Barrier( g.VCComm() );
        startTime = Time();
        blas::internal::GemmDot
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        Barrier( g.VCComm() );
        endTime = Time();
        runTime = endTime - startTime;
        gFlops = blas::internal::GemmGFlops<T>(m,n,k,runTime);
        if( g.VCRank() == 0 )
        {
            cout << "DONE. " << endl
                 << "  Time = " << runTime << " seconds. GFlops = " 
                 << gFlops << endl;
        }
        if( printMatrices )
        {
            ostringstream msg;
            msg << "C := " << alpha << " A B + " << beta << " C";
            C.Print( msg.str() );
        }
        if( testCorrectness )
        {
            TestParallelCorrectness
            ( printMatrices, C,
              orientationOfA, orientationOfB,
              alpha, ARef, BRef, beta, CRef );
        }
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 13 )
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
        const Orientation orientationOfA = CharToOrientation( *argv[3] );
        const Orientation orientationOfB = CharToOrientation( *argv[4] );
        const int         m = atoi( argv[5] );
        const int         n = atoi( argv[6] );
        const int         k = atoi( argv[7] );
        const int         nb = atoi( argv[8] );
        const bool        testSerial = atoi( argv[9] );
        const bool        testParallel = atoi( argv[10] );
        const bool        testCorrectness = atoi( argv[11] );
        const bool        printMatrices = atoi( argv[12] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 )
            cout << "Constructing grid..." << endl;
        const Grid g( MPI_COMM_WORLD, r, c );
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 )
            cout << "Done building grid." << endl;
        SetBlocksize( nb );
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 )
            cout << "Done setting blocksize." << endl;

        if( rank == 0 )
        {
            cout << "Will test Gemm" << OrientationToChar(orientationOfA) 
                                     << OrientationToChar(orientationOfB) 
                                     << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        if( testSerial )
        {
            TestSerialGemm<double>
            ( testCorrectness, printMatrices, 
              orientationOfA, orientationOfB,
              m, n, k, (double)3, (double)4, g );
        }
        if( testParallel )
        {
            TestParallelGemm<double>
            ( testCorrectness, printMatrices,
              orientationOfA, orientationOfB,
              m, n, k, (double)3, (double)4, g );
        }
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        if( testSerial )
        {
            TestSerialGemm<dcomplex>
            ( testCorrectness, printMatrices,
              orientationOfA, orientationOfB,
              m, n, k, (dcomplex)3, (dcomplex)4, g );
        }
        if( testParallel )
        {
            TestParallelGemm<dcomplex>
            ( testCorrectness, printMatrices,
              orientationOfA, orientationOfB,
              m, n, k, (dcomplex)3, (dcomplex)4, g );
        }
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

