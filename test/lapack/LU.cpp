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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Generates random matrix then solves for its LU factors.\n\n"
         << "  LU <r> <c> <m> <nb> <correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename T>
void TestCorrectness
( bool printMatrices,
  const DistMatrix<T,MC,MR>& A,
  const DistMatrix<int,VC,Star>& p,
  const DistMatrix<T,MC,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = AOrig.Height();
    DistMatrix<int,Star,Star> p_Star_Star(g);
    vector<int> image;
    vector<int> preimage;

    if( g.VCRank() == 0 )
        cout << "Testing error..." << endl;

    // Compose the pivots
    p_Star_Star = p;
    lapack::internal::ComposePivots( p_Star_Star, image, preimage, 0 );

    // Apply the pivots to our random right-hand sides
    DistMatrix<T,MC,MR> X(m,100,g);
    DistMatrix<T,MC,MR> Y(g);
    X.SetToRandom();
    T oneNormOfX = lapack::OneNorm( X );
    T infNormOfX = lapack::InfinityNorm( X );
    T frobNormOfX = lapack::FrobeniusNorm( X );
    Y = X;
    lapack::internal::ApplyRowPivots( Y, image, preimage, 0 );

    // Solve against the pivoted right-hand sides
    blas::Trsm( Left, Lower, Normal, Unit, (T)1, A, Y );
    blas::Trsm( Left, Upper, Normal, NonUnit, (T)1, A, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    blas::Gemm( Normal, Normal, (T)-1, AOrig, Y, (T)1, X );
    T oneNormOfError = lapack::OneNorm( X );
    T infNormOfError = lapack::InfinityNorm( X );
    T frobNormOfError = lapack::FrobeniusNorm( X );
    T oneNormOfA = lapack::OneNorm( AOrig );
    T infNormOfA = lapack::InfinityNorm( AOrig );
    T frobNormOfA = lapack::FrobeniusNorm( AOrig );

    if( g.VCRank() == 0 )
    {
        cout << "||A||_1                  = " << Abs(oneNormOfA) << "\n"
             << "||A||_oo                 = " << Abs(infNormOfA) << "\n"
             << "||A||_F                  = " << Abs(frobNormOfA) << "\n"
             << "||X||_1                  = " << Abs(oneNormOfX) << "\n"
             << "||X||_oo                 = " << Abs(infNormOfX) << "\n"
             << "||X||_F                  = " << Abs(frobNormOfX) << "\n"
             << "||A U^-1 L^-1 X - X||_1  = " << Abs(oneNormOfError) << "\n"
             << "||A U^-1 L^-1 X - X||_oo = " << Abs(infNormOfError) << "\n"
             << "||A U^-1 L^-1 X - X||_F  = " << Abs(frobNormOfError) << endl;
    }
}

template<typename T>
void TestLU
( bool testCorrectness, bool printMatrices,
  int m, const Grid& g )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(g);
    DistMatrix<T,MC,MR> ARef(g);
    DistMatrix<int,VC,Star> p(g);

    A.ResizeTo( m, m );
    p.ResizeTo( m, 1 );

    A.SetToRandom();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        ARef = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting LU factorization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::LU( A, p );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::LUGFlops<T>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
#ifdef TIMING
    if( g.VCRank() == 0 )
        elemental::lapack::lu::PrintTimings();
#endif
    if( printMatrices )
    {
        A.Print("A after factorization");
        p.Print("p after factorization");
    }
    if( testCorrectness )
    {
        TestCorrectness( printMatrices, A, p, ARef );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 7 )
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
            cout << "Will test LU" << endl;

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestLU<double>
        ( testCorrectness, printMatrices, m, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestLU<dcomplex>
        ( testCorrectness, printMatrices, m, g );
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

