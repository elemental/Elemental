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
    cout << "GEneral Matrix Matrix multiplication.\n\n"
         << "  Gemm <r> <c> <orient. of A?> <orient. of B?> <m> <n> <k> <nb>\n" 
         << "       <correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  orient. of A: {N,T,C}\n"
         << "  orient. of B: {N,T,C}\n"
         << "  m: height of C\n" 
         << "  n: width  of C\n"
         << "  k: inner dimension of AB\n"
         << "  nb: algorithmic blocksize\n"
         << "  correctness?: [0/1]\n"
         << "  print?: [0/1]\n" << endl; 
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
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,Star,Star>& ARef,
           const DistMatrix<T,Star,Star>& BRef,
  T beta,        DistMatrix<T,Star,Star>& CRef )
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
            T truth = CRef.GetLocalEntry(i,j);
            T computed = CCopy.GetLocalEntry(i,j);

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
void TestGemm
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
        cout << "  Starting Gemm...";
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
        TestCorrectness
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
        cout << "  Starting Gemm...";
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
        TestCorrectness
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
        cout << "  Starting Gemm...";
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
        TestCorrectness
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
            cout << "  Starting Gemm...";
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
            TestCorrectness
            ( printMatrices, C,
              orientationOfA, orientationOfB,
              alpha, ARef, BRef, beta, CRef );
        }
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 11 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    try
    {
        const int         r = atoi(argv[1]);
        const int         c = atoi(argv[2]);
        const Orientation orientationOfA = CharToOrientation(*argv[3]);
        const Orientation orientationOfB = CharToOrientation(*argv[4]);
        const int         m = atoi(argv[5]);
        const int         n = atoi(argv[6]);
        const int         k = atoi(argv[7]);
        const int         nb = atoi(argv[8]);
        const bool        testCorrectness = atoi(argv[9]);
        const bool        printMatrices = atoi(argv[10]);
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
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
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestGemm<double>
        ( testCorrectness, printMatrices,
          orientationOfA, orientationOfB,
          m, n, k, (double)3, (double)4, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestGemm<dcomplex>
        ( testCorrectness, printMatrices,
          orientationOfA, orientationOfB,
          m, n, k, (dcomplex)3, (dcomplex)4, g );
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

