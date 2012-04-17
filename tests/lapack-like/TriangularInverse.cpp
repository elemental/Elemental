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
    cout << "Inverts a triangular matrix.\n\n"
         << "  TriangularInverse <r> <c> <uplo> <diag> <m> <nb> "
            "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  uplo: {L,U}\n"
         << "  diag: {N,U}\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename F> // represents a real or complex number
void TestCorrectness
( bool printMatrices,
  UpperOrLower uplo, UnitOrNonUnit diag,
  const DistMatrix<F,MC,MR>& A,
  const DistMatrix<F,MC,MR>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    DistMatrix<F,MC,MR> X(g), Y(g);
    UniformRandom( m, 100, X );
    Y = X;

    // Since A o A^-1 = I, test the change introduced by the approximate comp.
    Trmm( LEFT, uplo, NORMAL, diag, (F)1, A,     Y );
    Trmm( LEFT, uplo, NORMAL, diag, (F)1, AOrig, Y );
    Axpy( (F)-1, X, Y );

    R oneNormOrig = Norm( AOrig, ONE_NORM );
    R infNormOrig = Norm( AOrig, INFINITY_NORM );
    R frobNormOrig = Norm( AOrig, FROBENIUS_NORM );
    R oneNormFinal = Norm( A, ONE_NORM );
    R infNormFinal = Norm( A, INFINITY_NORM );
    R frobNormFinal = Norm( A, FROBENIUS_NORM );
    R oneNormOfError = Norm( Y, ONE_NORM );
    R infNormOfError = Norm( Y, INFINITY_NORM );
    R frobNormOfError = Norm( Y, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "||A||_1           = " << oneNormOrig << "\n"
             << "||A||_oo          = " << infNormOrig << "\n"
             << "||A||_F           = " << frobNormOrig << "\n"
             << "||A^-1||_1        = " << oneNormFinal << "\n"
             << "||A^-1||_oo       = " << infNormFinal << "\n"
             << "||A^-1||_F        = " << frobNormFinal << "\n"
             << "||A A^-1 - I||_1  = " << oneNormOfError << "\n"
             << "||A A^-1 - I||_oo = " << infNormOfError << "\n"
             << "||A A^-1 - I||_F  = " << frobNormOfError << endl;
    }
}

template<typename F> // represents a real or complex number
void TestTriangularInverse
( bool testCorrectness, bool printMatrices,
  UpperOrLower uplo, UnitOrNonUnit diag, int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<F,MC,MR> A(g), AOrig(g);
    HPDUniformRandom( m, A );
    MakeTrapezoidal( LEFT, uplo, 0, A );
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
        cout << "  Starting triangular inversion...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    TriangularInverse( uplo, diag, A );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::TriangularInverseGFlops<F>( m, runTime );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after inversion");
    if( testCorrectness )
        TestCorrectness( printMatrices, uplo, diag, A, AOrig );
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
        const UnitOrNonUnit diag = CharToUnitOrNonUnit(*argv[++argNum]);
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
            cout << "Will test TriangularInverse" << UpperOrLowerToChar(uplo) 
                 << UnitOrNonUnitToChar(diag) << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTriangularInverse<double>
        ( testCorrectness, printMatrices, uplo, diag, m, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestTriangularInverse<Complex<double> >
        ( testCorrectness, printMatrices, uplo, diag, m, g );
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

