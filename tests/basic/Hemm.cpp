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
#include <ctime>
#include "elemental.hpp"
using namespace std;
using namespace elem;

void Usage()
{
    cout << "HErmitian Matrix Matrix multiplication.\n\n"
         << "  Hemm <r> <c> <side> <uplo> <m> <n> <nb> <print?>"
         << "\n\n"
         << "  r: number of process rows\n" 
         << "  c: number of process cols\n"
         << "  side: {L,R}\n"
         << "  uplo: {L,U}\n"
         << "  m: height of C\n"
         << "  n: width  of C\n"
         << "  nb: algorithmic blocksize\n"
         << "  print?: [0/1]\n" << endl;
}

template<typename T> // represents a real or complex ring
void TestHemm
( bool printMatrices, Side side, UpperOrLower uplo,
  int m, int n, T alpha, T beta, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g);
    DistMatrix<T,MC,MR> B(g);
    DistMatrix<T,MC,MR> C(g);

    if( side == LEFT )
        A.ResizeTo( m, m );
    else
        A.ResizeTo( n, n );
    B.ResizeTo( m, n );
    C.ResizeTo( m, n );

    // Test Hemm
    if( g.Rank() == 0 )
        cout << "Hemm:" << endl;
    A.SetToRandomHPD();
    B.SetToRandom();
    C.SetToRandom();

    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting Parallel Hemm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    Hemm( side, uplo, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::HemmGFlops<T>(side,m,n,runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        if( side == LEFT )
            msg << "C := " << alpha << " Herm(A) B + " << beta << " C";
        else
            msg << "C := " << alpha << " B Herm(A) + " << beta << " C";
        C.Print( msg.str() );
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 9 )
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
        const Side side = CharToSide(*argv[++argNum]);
        const UpperOrLower uplo = CharToUpperOrLower(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int n = atoi(argv[++argNum]);
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
        const Grid g( comm, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
        {
            cout << "Will test Hemm" << SideToChar(side) 
                                     << UpperOrLowerToChar(uplo) << endl;
        }

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with doubles:                 \n"
                 << "--------------------------------------" << endl;
        }
        TestHemm<double>
        ( printMatrices, side, uplo, m, n, (double)3, (double)4, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestHemm<Complex<double> >
        ( printMatrices, side, uplo, m, n, 
          Complex<double>(3), Complex<double>(4), g );
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

