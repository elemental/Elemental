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
    cout << "Generates SPD matrix then solves for its Cholesky factor.\n\n"
         << "  Cholesky <r> <c> <uplo> <m> <nb> <rankK local nb> "
            "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  uplo: {L,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  rankK local nb: local blocksize for triangular rank-k update\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename F> // represents a real or complex field
void TestCorrectness
( bool printMatrices, UpperOrLower uplo,
  const DistMatrix<F,MC,MR>& A,
  const DistMatrix<F,MC,MR>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F,MC,MR> X(g), Y(g);
    UniformRandom( m, 100, X );
    Y = X;

    if( uplo == LOWER )
    {
        // Test correctness by comparing the application of AOrig against a 
        // random set of 100 vectors to the application of tril(A) tril(A)^H
        Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, A, Y );
        Trmm( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, A, Y );
        Hemm( LEFT, LOWER, (F)-1, AOrig, X, (F)1, Y );
        R oneNormOfError = Norm( Y, ONE_NORM );
        R infNormOfError = Norm( Y, INFINITY_NORM );
        R frobNormOfError = Norm( Y, FROBENIUS_NORM );
        R infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        R frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
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
                 << "||A X - L L^H X||_1  = " << oneNormOfError << "\n"
                 << "||A X - L L^H X||_oo = " << infNormOfError << "\n"
                 << "||A X - L L^H X||_F  = " << frobNormOfError << endl;
        }
    }
    else
    {
        // Test correctness by comparing the application of AOrig against a 
        // random set of 100 vectors to the application of triu(A)^H triu(A)
        Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, Y );
        Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, A, Y );
        Hemm( LEFT, UPPER, (F)-1, AOrig, X, (F)1, Y );
        R oneNormOfError = Norm( Y, ONE_NORM );
        R infNormOfError = Norm( Y, INFINITY_NORM );
        R frobNormOfError = Norm( Y, FROBENIUS_NORM );
        R infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        R frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
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
                 << "||A X - U^H U X||_1  = " << oneNormOfError << "\n"
                 << "||A X - U^H U X||_oo = " << infNormOfError << "\n"
                 << "||A X - U^H U X||_F  = " << frobNormOfError << endl;
        }
    }
}

template<typename F> // represents a real or complex field
void TestCholesky
( bool testCorrectness, bool printMatrices, 
  UpperOrLower uplo, int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<F,MC,MR> A(g), AOrig(g);

    HPDUniformRandom( m, A );
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

    if( g.Rank() == 0 )
    {
        cout << "  Starting Cholesky factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    Cholesky( uplo, A );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::CholeskyGFlops<F>( m, runTime );
    if( g.Rank() == 0 )
    {
        cout << "DONE.\n"
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, uplo, A, AOrig );
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
        const int nbLocal = atoi(argv[++argNum]);
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
            cout << "Will test Cholesky" << UpperOrLowerToChar(uplo) << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestCholesky<double>
        ( testCorrectness, printMatrices, uplo, m, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestCholesky<Complex<double> >
        ( testCorrectness, printMatrices, uplo, m, g );
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

