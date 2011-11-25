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
using namespace std;
using namespace elemental;

void Usage()
{
    cout << "Reduced a Hermitian GEneralized EVP to Hermitian STandard EVP\n\n"
         << "  Hegst <r> <c> <side> <uplo> <m> <nb> <correctness?> "
            "<print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  side: we set up ABX=XW or BAX=XW if side=LEFT, \n"
         << "                            AX=BXW if side=RIGHT\n"
         << "  uplo: {L,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename F> // represents a real or complex field
void TestCorrectness
( bool printMatrices, Side side, UpperOrLower uplo,
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

    if( side == RIGHT )
    {
        if( uplo == LOWER )
        {
            // Test correctness by comparing the application of A against a 
            // random set of 100 vectors to the application of 
            // tril(B)^-1 AOrig tril(B)^-H
            basic::Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, B, Y );
            basic::Hemm( LEFT, LOWER, (F)1, AOrig, Y, (F)0, Z );
            basic::Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, B, Z );
            basic::Hemm( LEFT, LOWER, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
            F frobNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
            F infNormOfA = advanced::HermitianNorm( uplo, A, INFINITY_NORM );
            F frobNormOfA = advanced::HermitianNorm( uplo, A, FROBENIUS_NORM );
            F oneNormOfError = advanced::Norm( Z, ONE_NORM );
            F infNormOfError = advanced::Norm( Z, INFINITY_NORM );
            F frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
            if( g.Rank() == 0 )
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
            basic::Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, B, Y );
            basic::Hemm( LEFT, UPPER, (F)1, AOrig, Y, (F)0, Z );
            basic::Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, B, Z );
            basic::Hemm( LEFT, UPPER, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
            F frobNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
            F infNormOfA = advanced::HermitianNorm( uplo, A, INFINITY_NORM );
            F frobNormOfA = advanced::HermitianNorm( uplo, A, FROBENIUS_NORM );
            F oneNormOfError = advanced::Norm( Z, ONE_NORM );
            F infNormOfError = advanced::Norm( Z, INFINITY_NORM );
            F frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
            if( g.Rank() == 0 )
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
        if( uplo == LOWER )
        {
            // Test correctness by comparing the application of A against a 
            // random set of 100 vectors to the application of 
            // tril(B)^H AOrig tril(B)
            basic::Trmm( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, B, Y );
            basic::Hemm( LEFT, LOWER, (F)1, AOrig, Y, (F)0, Z );
            basic::Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, B, Z );
            basic::Hemm( LEFT, LOWER, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
            F frobNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
            F infNormOfA = 
                advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
            F frobNormOfA = 
                advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
            F oneNormOfError = advanced::Norm( Z, ONE_NORM );
            F infNormOfError = advanced::Norm( Z, INFINITY_NORM );
            F frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
            if( g.Rank() == 0 )
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
            basic::Trmm( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, B, Y );
            basic::Hemm( LEFT, UPPER, (F)1, AOrig, Y, (F)0, Z );
            basic::Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, B, Z );
            basic::Hemm( LEFT, UPPER, (F)-1, A, X, (F)1, Z );
            F infNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
            F frobNormOfAOrig = 
                advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
            F infNormOfA = 
                advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
            F frobNormOfA = 
                advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
            F oneNormOfError = advanced::Norm( Z, ONE_NORM );
            F infNormOfError = advanced::Norm( Z, INFINITY_NORM );
            F frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
            if( g.Rank() == 0 )
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
  Side side, UpperOrLower uplo, int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<F,MC,MR> A(g);
    DistMatrix<F,MC,MR> B(g);
    DistMatrix<F,MC,MR> AOrig(g);

    A.ResizeTo( m, m );
    B.ResizeTo( m, m );

    if( testCorrectness )
    {
        DistMatrix<F,MC,MR> C( m, m, g );
        C.SetToRandom();
        basic::Herk( uplo, NORMAL, (F)1, C, (F)0, A );
        C.SetToRandom();
        basic::Herk( uplo, NORMAL, (F)1, C, (F)0, B );
        advanced::Cholesky( uplo, B );
    }
    else
    {
        A.SetToRandomHPD();
        B.SetToRandomHPD();
    }
    B.MakeTrapezoidal( LEFT, uplo );
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
    {
        A.Print("A");
        B.Print("B");
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting reduction to Hermitian standard EVP...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    advanced::Hegst( side, uplo, A, B );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = advanced::internal::HegstGFlops<F>( m, runTime );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = "
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after reduction");
    if( testCorrectness )
        TestCorrectness( printMatrices, side, uplo, A, B, AOrig );
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
        const Grid g( comm, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
        {
            cout << "Will test Hegst" << SideToChar(side) 
                 << UpperOrLowerToChar(uplo) << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestHegst<double>
        ( testCorrectness, printMatrices, side, uplo, m, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestHegst<dcomplex>
        ( testCorrectness, printMatrices, side, uplo, m, g );
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

