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
    cout << "Generates random matrix then solves for its LQ factorization.\n\n"
         << "  LQ <r> <c> <m> <n> <nb> <correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  m: height of matrix\n"
         << "  n: width of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename R> // represents a real number
void TestCorrectness
( bool printMatrices,
  const DistMatrix<R,MC,MR>& A,
        DistMatrix<R,MC,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int minDim = std::min(m,n);

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<R,MC,MR> Z(m,n,g);
    MakeIdentity( Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, FORWARD, 0, A, Z );

    DistMatrix<R,MC,MR> ZUpper(g);
    ZUpper.View( Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<R,MC,MR> X(minDim,minDim,g);
    MakeIdentity( X );

    // Form X := I - Q^H Q
    Axpy( (R)-1, ZUpper, X );

    R oneNormOfError = Norm( X, ONE_NORM );
    R infNormOfError = Norm( X, INFINITY_NORM );
    R frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
             << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
    }

    if( g.Rank() == 0 )
    {
        cout << "  Testing if A = LQ..." << endl;
    }

    // Form L Q
    DistMatrix<R,MC,MR> L( A );
    MakeTrapezoidal( LEFT, LOWER, 0, L );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, L );

    // Form L Q - A
    Axpy( (R)-1, AOrig, L );
    
    R oneNormOfA = Norm( AOrig, ONE_NORM );
    R infNormOfA = Norm( AOrig, INFINITY_NORM );
    R frobNormOfA = Norm( AOrig, FROBENIUS_NORM );
    oneNormOfError = Norm( L, ONE_NORM );
    infNormOfError = Norm( L, INFINITY_NORM );
    frobNormOfError = Norm( L, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - LQ||_1  = " << oneNormOfError << "\n"
             << "    ||A - LQ||_oo = " << infNormOfError << "\n"
             << "    ||A - LQ||_F  = " << frobNormOfError << endl;
    }
}

template<typename R> // represents a real number
void TestCorrectness
( bool printMatrices,
  const DistMatrix<Complex<R>,MC,MR  >& A,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R>,MC,MR  >& AOrig )
{
    typedef Complex<R> C;

    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int minDim = std::min(m,n);

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<C,MC,MR> Z(m,n,g);
    MakeIdentity( Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, CONJUGATED, 0, A, t, Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, FORWARD, CONJUGATED, 0, A, t, Z );
    
    DistMatrix<C,MC,MR> ZUpper(g);
    ZUpper.View( Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<C,MC,MR> X(minDim,minDim,g);
    MakeIdentity( X );

    // Form X := I - Q^H Q
    Axpy( (C)-1, ZUpper, X );

    R oneNormOfError = Norm( X, ONE_NORM );
    R infNormOfError = Norm( X, INFINITY_NORM );
    R frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
             << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
    }

    if( g.Rank() == 0 )
        cout << "  Testing if A = LQ..." << endl;

    // Form L Q
    DistMatrix<C,MC,MR> L( A );
    MakeTrapezoidal( LEFT, LOWER, 0, L );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, CONJUGATED, 0, A, t, L );

    // Form L Q - A
    Axpy( (C)-1, AOrig, L );
    
    R oneNormOfA = Norm( AOrig, ONE_NORM );
    R infNormOfA = Norm( AOrig, INFINITY_NORM );
    R frobNormOfA = Norm( AOrig, FROBENIUS_NORM );
    oneNormOfError = Norm( L, ONE_NORM );
    infNormOfError = Norm( L, INFINITY_NORM );
    frobNormOfError = Norm( L, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - LQ||_1  = " << oneNormOfError << "\n"
             << "    ||A - LQ||_oo = " << infNormOfError << "\n"
             << "    ||A - LQ||_F  = " << frobNormOfError << endl;
    }
}

template<typename F> // represents a real or complex field
void TestLQ
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g );

template<>
void TestLQ<double>
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> A(g), AOrig(g);
    UniformRandom( m, n, A );

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
        cout << "  Starting LQ factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    LQ( A );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::LQGFlops<R>( m, n, runTime );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, AOrig );
}

template<>
void TestLQ<Complex<double> >
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g )
{
    typedef Complex<double> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,MR  > A(g), AOrig(g);
    UniformRandom( m, n, A );

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
    DistMatrix<C,MD,STAR> t(g);

    if( g.Rank() == 0 )
    {
        cout << "  Starting LQ factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    LQ( A, t );
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    gFlops = internal::LQGFlops<C>( m, n, runTime );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, t, AOrig );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 8 )
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
        const int n = atoi(argv[++argNum]);
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
            cout << "Will test LQ" << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestLQ<double>( testCorrectness, printMatrices, m, n, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestLQ<Complex<double> >( testCorrectness, printMatrices, m, n, g );
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

