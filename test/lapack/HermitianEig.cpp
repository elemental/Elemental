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
    cout << "Generates random HPD matrix then solves for its eigenpairs.\n\n"
         << "  HermitianEig <r> <c> <shape> <m> <nb> <correctness?> "
         << "<print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  shape: Upper iff 1\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

template<typename R>
void TestCorrectness
( bool printMatrices,
  Shape shape,
  const DistMatrix<R,MC,  MR>& A,
  const DistMatrix<R,Star,VR>& w,
  const DistMatrix<R,MC,  MR>& Z,
  const DistMatrix<R,MC  ,MR>& ARef )
{
    const Grid& g = A.GetGrid();
    const int k = Z.Width();

    if( printMatrices )
    {
        w.Print("Computed eigenvalues:");
        Z.Print("Computed eigenvectors:");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<R,Star,Star> w_Star_Star(g); 
    w_Star_Star = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Testing orthogonality of eigenvectors...";
        cout.flush();
    }
    DistMatrix<R,MC,MR> X(k,k,g);
    X.SetToIdentity();
    blas::Herk( shape, ConjugateTranspose, (R)-1, Z, (R)1, X );
    R myResidual = 0; 
    for( int j=0; j<X.LocalWidth(); ++j )
        for( int i=0; i<X.LocalHeight(); ++i )
            myResidual = max( Abs(X.GetLocalEntry(i,j)),myResidual );
    R residual;
    Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "max deviation from I is " << residual << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Testing for deviation of AX from WX...";
        cout.flush();
    }
    // Set X := WZ, where W is the diagonal eigenvalue matrix
    X = Z;
    for( int j=0; j<Z.LocalWidth(); ++j )
    {
        R omega = w.GetLocalEntry(0,j);
        elemental::wrappers::blas::Scal
        ( X.LocalHeight(), omega, X.LocalBuffer(0,j), 1 );
    }
    blas::Gemm( Normal, Normal, (R)-1, ARef, Z, (R)1, X );
    myResidual = 0;
    for( int j=0; j<X.LocalWidth(); ++j )
        for( int i=0; i<X.LocalHeight(); ++i )
            myResidual = max( Abs(X.GetLocalEntry(i,j)),myResidual );
    Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "max deviation of AZ from WZ is " << residual << endl;
}

template<typename R>
void TestHermitianEig
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<R,MC,MR> A(m,m,g);
    DistMatrix<R,MC,MR> ARef(m,m,g);
    DistMatrix<R,Star,VR> w(1,m,g);
    DistMatrix<R,MC,MR> Z(m,m,g);

    A.SetToRandomHPD();
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
    {
        A.Print("A");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting Hermitian eigensolver.\n";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::HermitianEig( shape, A, w, Z );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << endl;
    }
    if( printMatrices )
    {
        A.Print("A after eigensolve");
    }
    if( testCorrectness )
    {
        TestCorrectness( printMatrices, shape, A, w, Z, ARef );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 8 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    try
    {
        const int r = atoi(argv[1]);
        const int c = atoi(argv[2]);
        const Shape shape = CharToShape(*argv[3]);
        const int m = atoi(argv[4]);
        const int nb = atoi(argv[5]);
        const bool testCorrectness = atoi(argv[6]);
        const bool printMatrices = atoi(argv[7]);
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
            cout << "Will test " << ( shape==Lower ? "lower" : "upper" )
                 << " HermitianEig." << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestHermitianEig<double>
        ( testCorrectness, printMatrices, shape, m, g );
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

