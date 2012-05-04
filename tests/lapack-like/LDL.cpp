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
    cout << "Generates symmetric matrix then solves for its LDL^T or LDL^H "
            "factorization.\n\n"
         << "  LDL <r> <c> <m> <nb> <rankK local nb> <conjugated?> "
            "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  rankK local nb: local blocksize for triangular rank-k update\n"
         << "  conjugated?: false iff 0\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename F> // represents a real or complex field
void TestCorrectness
( bool conjugated, bool printMatrices, 
  const DistMatrix<F,MC,MR>& A,
  const DistMatrix<F,MC,STAR>& d,
  const DistMatrix<F,MC,MR>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F,MC,MR> X(g), Y(g);
    Uniform( m, 100, X );
    Y = X;

    // Test correctness by comparing the application of AOrig against a 
    // random set of 100 vectors to the application of tril(A) tril(A)^H
    if( conjugated )
        Trmm( LEFT, LOWER, ADJOINT, UNIT, (F)1, A, Y );
    else
        Trmm( LEFT, LOWER, TRANSPOSE, UNIT, (F)1, A, Y );
    DiagonalScale( LEFT, NORMAL, d, Y );
    Trmm( LEFT, LOWER, NORMAL, UNIT, (F)1, A, Y );
    if( conjugated )
        Hemm( LEFT, LOWER, (F)-1, AOrig, X, (F)1, Y );
    else
        Symm( LEFT, LOWER, (F)-1, AOrig, X, (F)1, Y );
    R oneNormOfError = Norm( Y, ONE_NORM );
    R infNormOfError = Norm( Y, INFINITY_NORM );
    R frobNormOfError = Norm( Y, FROBENIUS_NORM );
    R infNormOfA = HermitianNorm( LOWER, AOrig, INFINITY_NORM );
    R frobNormOfA = HermitianNorm( LOWER, AOrig, FROBENIUS_NORM );
    R oneNormOfX = Norm( X, ONE_NORM );
    R infNormOfX = Norm( X, INFINITY_NORM );
    R frobNormOfX = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "||A||_1 = ||A||_oo   = " << infNormOfA << "\n"
             << "||A||_F              = " << frobNormOfA << "\n"
             << "||X||_1              = " << oneNormOfX << "\n"
             << "||X||_oo             = " << infNormOfX << "\n"
             << "||X||_F              = " << frobNormOfX << "\n"
             << "||A X - L D L^[T/H] X||_1  = " << oneNormOfError << "\n"
             << "||A X - L D L^[T/H] X||_oo = " << infNormOfError << "\n"
             << "||A X - L D L^[T/H] X||_F  = " << frobNormOfError << endl;
    }
}

template<typename F> // represents a real or complex field
void TestLDL
( bool conjugated, bool testCorrectness, bool printMatrices, 
  int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;

    DistMatrix<F,MC,MR> A(g), AOrig(g);
    if( conjugated )
        HermitianUniformSpectrum( m, A, -100, 100 );
    else
        Uniform( m, m, A );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");
    DistMatrix<F,MC,STAR> d(g);

    if( g.Rank() == 0 )
    {
        cout << "  Starting LDL^[T/H] factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    if( !conjugated )
        LDLT( A, d );
    else
        LDLH( A, d );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    if( conjugated )
        gFlops = internal::LDLHGFlops<F>( m, runTime );
    else
        gFlops = internal::LDLTGFlops<F>( m, runTime );
    if( g.Rank() == 0 )
    {
        cout << "DONE.\n"
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( conjugated, printMatrices, A, d, AOrig );
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
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const int nbLocal = atoi(argv[++argNum]);
        const bool conjugated = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
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
        SetLocalTrrkBlocksize<double>( nbLocal );
        SetLocalTrrkBlocksize<Complex<double> >( nbLocal );

        if( rank == 0 )
            cout << "Will test LDL" << (conjugated?"^H":"^T") << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestLDL<double>
        ( conjugated, testCorrectness, printMatrices, m, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestLDL<Complex<double> >
        ( conjugated, testCorrectness, printMatrices, m, g );
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

