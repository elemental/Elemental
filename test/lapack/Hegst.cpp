/*
   Copyright (c) 2009-2011, Jack Poulson
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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Reduced a Hermitian GEneralized EVP to Hermitian STandard EVP\n\n"
         << "  Hegst <r> <c> <side> <shape> <m> <nb> <correctness?> "
            "<print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  side: we set up ABX=XW or BAX=XW if side=Left, \n"
         << "                            AX=BXW if side=Right\n"
         << "  shape: {L,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename F> // represents a real or complex field
void TestCorrectness
( bool printMatrices, Side side, Shape shape,
  const DistMatrix<F,MC,MR>& A,
  const DistMatrix<F,MC,MR>& B,
  const DistMatrix<F,MC,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F,MC,MR> X(m,100,g);
    DistMatrix<F,MC,MR> Y(m,100,g);
    DistMatrix<F,MC,MR> Z(m,100,g);
    X.SetToRandom();
    Y = X;

    if( side == Right )
    {
        if( shape == Lower )
        {
            // Test correctness by comparing the application of A against a 
            // random set of 100 vectors to the application of 
            // tril(B)^-1 AOrig tril(B)^-H
            blas::Trsm( Left, Lower, ConjugateTranspose, NonUnit, (F)1, B, Y );
            blas::Hemm( Left, Lower, (F)1, AOrig, Y, (F)0, Z );
            blas::Trsm( Left, Lower, Normal, NonUnit, (F)1, B, Z );
            blas::Hemm( Left, Lower, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = lapack::HermitianInfinityNorm( shape, AOrig );
            F frobNormOfAOrig = lapack::HermitianFrobeniusNorm( shape, AOrig );
            F infNormOfA = lapack::HermitianInfinityNorm( shape, A );
            F frobNormOfA = lapack::HermitianFrobeniusNorm( shape, A );
            F oneNormOfError = lapack::OneNorm( Z );
            F infNormOfError = lapack::InfinityNorm( Z );
            F frobNormOfError = lapack::FrobeniusNorm( Z );
            if( g.VCRank() == 0 )
            {
                cout << "||AOrig||_1 = ||AOrig||_oo     = "
                     << Abs(infNormOfAOrig) << "\n"
                     << "||AOrig||_F                    = "
                     << Abs(frobNormOfAOrig) << "\n"
                     << "||A||_1 = ||A||_oo             = "
                     << Abs(infNormOfA) << "\n"
                     << "||A||_F                        = "
                     << Abs(frobNormOfA) << "\n"
                     << "||A X - L^-1 AOrig L^-H X||_1  = "
                     << Abs(oneNormOfError) << "\n"
                     << "||A X - L^-1 AOrig L^-H X||_oo = " 
                     << Abs(infNormOfError) << "\n"
                     << "||A X - L^-1 AOrig L^-H X||_F  = "
                     << Abs(frobNormOfError) << endl;
            }
        }
        else
        {
            // Test correctness by comparing the application of A against a 
            // random set of 100 vectors to the application of 
            // triu(B)^-H AOrig triu(B)^-1
            blas::Trsm( Left, Upper, Normal, NonUnit, (F)1, B, Y );
            blas::Hemm( Left, Upper, (F)1, AOrig, Y, (F)0, Z );
            blas::Trsm( Left, Upper, ConjugateTranspose, NonUnit, (F)1, B, Z );
            blas::Hemm( Left, Upper, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = lapack::HermitianInfinityNorm( shape, AOrig );
            F frobNormOfAOrig = lapack::HermitianFrobeniusNorm( shape, AOrig );
            F infNormOfA = lapack::HermitianInfinityNorm( shape, A );
            F frobNormOfA = lapack::HermitianFrobeniusNorm( shape, A );
            F oneNormOfError = lapack::OneNorm( Z );
            F infNormOfError = lapack::InfinityNorm( Z );
            F frobNormOfError = lapack::FrobeniusNorm( Z );
            if( g.VCRank() == 0 )
            {
                cout << "||AOrig||_1 = ||AOrig||_oo     = "
                     << Abs(infNormOfAOrig) << "\n"
                     << "||AOrig||_F                    = "
                     << Abs(frobNormOfAOrig) << "\n"
                     << "||A||_1 = ||A||_oo             = "
                     << Abs(infNormOfA) << "\n"
                     << "||A||_F                        = "
                     << Abs(frobNormOfA) << "\n"
                     << "||A X - U^-H AOrig U^-1 X||_1  = "
                     << Abs(oneNormOfError) << "\n"
                     << "||A X - U^-H AOrig U^-1 X||_oo = " 
                     << Abs(infNormOfError) << "\n"
                     << "||A X - U^-H AOrig U^-1 X||_F  = "
                     << Abs(frobNormOfError) << endl;
            }
        }
    }
    else
    {
        if( shape == Lower )
        {
            // Test correctness by comparing the application of A against a 
            // random set of 100 vectors to the application of 
            // tril(B)^H AOrig tril(B)
            blas::Trmm( Left, Lower, Normal, NonUnit, (F)1, B, Y );
            blas::Hemm( Left, Lower, (F)1, AOrig, Y, (F)0, Z );
            blas::Trmm( Left, Lower, ConjugateTranspose, NonUnit, (F)1, B, Z );
            blas::Hemm( Left, Lower, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = lapack::HermitianInfinityNorm( shape, AOrig );
            F frobNormOfAOrig = lapack::HermitianFrobeniusNorm( shape, AOrig );
            F infNormOfA = lapack::HermitianInfinityNorm( shape, AOrig );
            F frobNormOfA = lapack::HermitianFrobeniusNorm( shape, AOrig );
            F oneNormOfError = lapack::OneNorm( Z );
            F infNormOfError = lapack::InfinityNorm( Z );
            F frobNormOfError = lapack::FrobeniusNorm( Z );
            if( g.VCRank() == 0 )
            {
                cout << "||AOrig||_1 = ||AOrig||_oo = "
                     << Abs(infNormOfAOrig) << "\n"
                     << "||AOrig||_F                = "
                     << Abs(frobNormOfAOrig) << "\n"
                     << "||A||_1 = ||A||_F          = "
                     << Abs(infNormOfA) << "\n"
                     << "||A||_F                    = "
                     << Abs(frobNormOfA) << "\n"
                     << "||A X - L^H AOrig L X||_1  = "
                     << Abs(oneNormOfError) << "\n"
                     << "||A X - L^H AOrig L X||_oo = " 
                     << Abs(infNormOfError) << "\n"
                     << "||A X - L^H AOrig L X||_F  = "
                     << Abs(frobNormOfError) << endl;
            }
        }
        else
        {
            // Test correctness by comparing the application of A against a 
            // random set of 100 vectors to the application of 
            // triu(B) AOrig triu(B)^H
            blas::Trmm( Left, Upper, ConjugateTranspose, NonUnit, (F)1, B, Y );
            blas::Hemm( Left, Upper, (F)1, AOrig, Y, (F)0, Z );
            blas::Trmm( Left, Upper, Normal, NonUnit, (F)1, B, Z );
            blas::Hemm( Left, Upper, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = lapack::HermitianInfinityNorm( shape, AOrig );
            F frobNormOfAOrig = lapack::HermitianFrobeniusNorm( shape, AOrig );
            F infNormOfA = lapack::HermitianInfinityNorm( shape, AOrig );
            F frobNormOfA = lapack::HermitianFrobeniusNorm( shape, AOrig );
            F oneNormOfError = lapack::OneNorm( Z );
            F infNormOfError = lapack::InfinityNorm( Z );
            F frobNormOfError = lapack::FrobeniusNorm( Z );
            if( g.VCRank() == 0 )
            {
                cout << "||AOrig||_1 = ||AOrig||_oo = "
                     << Abs(infNormOfAOrig) << "\n"
                     << "||AOrig||_F                = "
                     << Abs(frobNormOfAOrig) << "\n"
                     << "||A||_1 = ||A||_F          = "
                     << Abs(infNormOfA) << "\n"
                     << "||A||_F                    = "
                     << Abs(frobNormOfA) << "\n"
                     << "||A X - U AOrig U^H X||_1  = "
                     << Abs(oneNormOfError) << "\n"
                     << "||A X - U AOrig U^H X||_oo = " 
                     << Abs(infNormOfError) << "\n"
                     << "||A X - U AOrig U^H X||_F  = "
                     << Abs(frobNormOfError) << endl;
            }
        }
    }
}

template<typename F> // represents a real or complex field
void TestHegst
( bool testCorrectness, bool printMatrices,
  Side side, Shape shape, int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<F,MC,MR> A(g);
    DistMatrix<F,MC,MR> B(g);
    DistMatrix<F,MC,MR> AOrig(g);

    A.ResizeTo( m, m );
    B.ResizeTo( m, m );

    A.SetToRandomHPD();
    B.SetToRandomHPD();
    B.MakeTrapezoidal( Left, shape );
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting reduction to Hermitian standard EVP...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Hegst( side, shape, A, B );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::HegstGFlops<F>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = "
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after reduction");
    if( testCorrectness )
        TestCorrectness( printMatrices, side, shape, A, B, AOrig );
}

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 9 )
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
        const Shape shape = CharToShape(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
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
        const Grid g( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
        {
            cout << "Will test Hegst" << SideToChar(side) << ShapeToChar(shape)
                 << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestHegst<double>
        ( testCorrectness, printMatrices, side, shape, m, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestHegst<dcomplex>
        ( testCorrectness, printMatrices, side, shape, m, g );
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

