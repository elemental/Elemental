/*
   Copyright (c) 2009-2012, Jack Poulson
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
#include "elemental.hpp"
using namespace std;
using namespace elem;

void Usage()
{
    cout << "GEneral Matrix Matrix multiplication.\n\n"
         << "  Gemm <r> <c> <orient. of A?> <orient. of B?> <m> <n> <k> <nb> " 
            "<print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  orient. of A: {N,T,C}\n"
         << "  orient. of B: {N,T,C}\n"
         << "  m: height of C\n" 
         << "  n: width  of C\n"
         << "  k: inner dimension of AB\n"
         << "  nb: algorithmic blocksize\n"
         << "  print?: [0/1]\n" << endl; 
}

template<typename T> // represents a real or complex ring
void TestGemm
( bool printMatrices, Orientation orientationOfA, Orientation orientationOfB,
  int m, int n, int k, T alpha, T beta, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g), B(g), C(g);

    if( orientationOfA == NORMAL )
        A.ResizeTo( m, k );
    else
        A.ResizeTo( k, m );
    if( orientationOfB == NORMAL )
        B.ResizeTo( k, n );
    else
        B.ResizeTo( n, k );
    C.ResizeTo( m, n );

    // Test the variant of Gemm that keeps A stationary
    if( g.Rank() == 0 )
        cout << "Stationary A Algorithm:" << endl;
    MakeUniformRandom( A );
    MakeUniformRandom( B );
    MakeUniformRandom( C );
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Gemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    internal::GemmA
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::GemmGFlops<T>(m,n,k,runTime);
    if( g.Rank() == 0 )
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

    // Test the variant of Gemm that keeps B stationary
    if( g.Rank() == 0 )
        cout << endl << "Stationary B Algorithm:" << endl;
    MakeUniformRandom( A );
    MakeUniformRandom( B );
    MakeUniformRandom( C );
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Gemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    internal::GemmB
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::GemmGFlops<T>(m,n,k,runTime);
    if( g.Rank() == 0 )
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

    // Test the variant of Gemm that keeps C stationary
    if( g.Rank() == 0 )
        cout << endl << "Stationary C Algorithm:" << endl;
    MakeUniformRandom( A );
    MakeUniformRandom( B );
    MakeUniformRandom( C );
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Gemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    internal::GemmC
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::GemmGFlops<T>(m,n,k,runTime);
    if( g.Rank() == 0 )
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
    
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        // Test the variant of Gemm for panel-panel dot products
        if( g.Rank() == 0 )
            cout << endl << "Dot Product Algorithm:" << endl;
        MakeUniformRandom( A );
        MakeUniformRandom( B );
        MakeUniformRandom( C );
        if( printMatrices )
        {
            A.Print("A");
            B.Print("B");
            C.Print("C");
        }
        if( g.Rank() == 0 )
        {
            cout << "  Starting Gemm...";
            cout.flush();
        }
        mpi::Barrier( g.Comm() );
        startTime = mpi::Time();
        internal::GemmDot
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        mpi::Barrier( g.Comm() );
        endTime = mpi::Time();
        runTime = endTime - startTime;
        gFlops = internal::GemmGFlops<T>(m,n,k,runTime);
        if( g.Rank() == 0 )
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
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 10 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try
    {
        int argNum = 0;
        const int r = atoi(argv[++argNum]);
        const int c = atoi(argv[++argNum]);
        const Orientation orientationOfA = CharToOrientation(*argv[++argNum]);
        const Orientation orientationOfB = CharToOrientation(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int n = atoi(argv[++argNum]);
        const int k = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        mpi::Barrier( comm );
        const Grid g( comm, r, c );
        mpi::Barrier( comm );
        SetBlocksize( nb );
        mpi::Barrier( comm );

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
        ( printMatrices, orientationOfA, orientationOfB,
          m, n, k, (double)3, (double)4, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestGemm<Complex<double> >
        ( printMatrices, orientationOfA, orientationOfB,
          m, n, k, Complex<double>(3), Complex<double>(4), g );
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

