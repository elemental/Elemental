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
    cout << "SYmmetric Matrix vector multiplication.\n\n"
         << "  Symv <r> <c> <Shape> <m> <nb> <local nb double> "
         << "<local nb complex double> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  Shape: {L,U}\n"
         << "  m: height of C\n"
         << "  nb: algorithmic blocksize\n"
         << "  local nb double: local algorithmic blocksize for double-prec.\n"
         << "  local nb complex double: \" \" complex double-precision\n"
         << "  print?: [0/1]\n" << endl;
}

template<typename T> // represents a real or complex ring
void TestSymv
( const UpperOrLower uplo, const int m, const T alpha, const T beta, 
  const bool printMatrices, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g), x(g), y(g);

    Uniform( m, m, A );
    Uniform( m, 1, x );
    Uniform( m, 1, y ); 
    if( printMatrices )
    {
        A.Print("A");
        x.Print("x");
        y.Print("y");
    }

    // Test Symm
    if( g.Rank() == 0 )
    {
        cout << "  Starting Parallel Symv...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    Symv( uplo, alpha, A, x, beta, y );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::SymvGFlops<T>(m,runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        ostringstream msg;
        msg << "y := " << alpha << " Symm(A) x + " << beta << " y";
        y.Print( msg.str() );
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
        const UpperOrLower uplo = CharToUpperOrLower(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const int nbLocalDouble = atoi(argv[++argNum]);
        const int nbLocalComplexDouble = atoi(argv[++argNum]);
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
        SetLocalSymvBlocksize<double>( nbLocalDouble );
        SetLocalSymvBlocksize<Complex<double> >( nbLocalComplexDouble );

        if( rank == 0 )
            cout << "Will test Symv" << UpperOrLowerToChar(uplo) << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestSymv<double>
        ( uplo, m, (double)3, (double)4, printMatrices, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestSymv<Complex<double> >
        ( uplo, m, Complex<double>(3), Complex<double>(4), printMatrices, g );
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:" << endl 
             << e.what() << endl;
    }
    Finalize();
    return 0;
}

