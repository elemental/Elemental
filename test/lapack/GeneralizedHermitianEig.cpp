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
    cout << "Generates random Hermitian A and random HPD B then solves for "
         << "their eigenpairs.\n\n"
         << "  GeneralizedHermitianEig <r> <c> <side> <shape> <m> <nb> "
         << "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  side: L/R\n"
         << "  shape: L/U\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

void TestCorrectnessDouble
( bool printMatrices,
  Side side,
  Shape shape,
  const DistMatrix<double,MC,  MR>& A,
  const DistMatrix<double,MC,  MR>& B,
  const DistMatrix<double,Star,VR>& w,
  const DistMatrix<double,MC,  MR>& Z,
  const DistMatrix<double,MC  ,MR>& ARef,
  const DistMatrix<double,MC,  MR>& BRef )
{
    const Grid& g = A.GetGrid();
    const int n = Z.Height();
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
    DistMatrix<double,Star,MR> w_Star_MR(g); 
    w_Star_MR.AlignWith( Z );
    w_Star_MR = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( side == Right )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of AZ from BZW...";
            cout.flush();
        }
        // Set X := BZW, where W is the diagonal eigenvalue matrix
        DistMatrix<double,MC,MR> X( g );
        X.AlignWith( Z );
        X.ResizeTo( n, k );
        blas::Hemm( Left, shape, (double)1, BRef, Z, (double)0, X );
        for( int j=0; j<X.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j);
            elemental::wrappers::blas::Scal
            ( X.LocalHeight(), omega, X.LocalBuffer(0,j), 1 );
        }
        // X := X - AZ = BZW - AZ
        blas::Hemm( Left, shape, (double)-1, ARef, Z, (double)1, X );
        // Compute the residual, ||BZW-AZ||_oo
        double myResidual = 0;
        for( int j=0; j<X.LocalWidth(); ++j )
            for( int i=0; i<X.LocalHeight(); ++i )
                myResidual = max( Abs(X.GetLocalEntry(i,j)),myResidual );
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||A Z - B Z W||_oo = " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        DistMatrix<double,MC,MR> Y(g);
        Y = Z;
        if( shape == Lower )
        {
            blas::Trmm
            ( Left, Lower, ConjugateTranspose, NonUnit, (double)1, B, Y );
        }
        else
        {
            blas::Trmm
            ( Left, Upper, Normal, NonUnit, (double)1, B, Y );
        }
        X.ResizeTo( k, k );
        X.SetToIdentity();
        blas::Herk( shape, ConjugateTranspose, (double)-1, Y, (double)1, X );
        myResidual = 0; 
        for( int j=0; j<X.LocalWidth(); ++j )
            for( int i=0; i<X.LocalHeight(); ++i )
                myResidual = max( Abs(X.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||Z^H B Z - I||_oo = " << residual << endl;
    }
    else
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of BAZ from ZW...";
            cout.flush();
        }
        // Set X := AZ
        DistMatrix<double,MC,MR> X( g );
        X.AlignWith( Z );
        X.ResizeTo( n, k );
        blas::Hemm( Left, shape, (double)1, ARef, Z, (double)0, X );
        // Set Y := BX = BAZ
        DistMatrix<double,MC,MR> Y( n, k, g );
        blas::Hemm( Left, shape, (double)1, BRef, X, (double)0, Y );
        // Compute the residual, ||Y - ZW||_oo = ||BAZ - ZW||_oo
        double myResidual = 0;
        for( int j=0; j<Y.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Y.LocalHeight(); ++i )
            {
                double thisResidual = 
                    Abs(omega*Z.GetLocalEntry(i,j)-Y.GetLocalEntry(i,j));
                myResidual = max( thisResidual, myResidual );
            }
        }
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||B A Z - Z W||_oo =  " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        Y = Z;
        if( shape == Lower )
        {
            blas::Trsm
            ( Left, Lower, Normal, NonUnit, (double)1, B, Y );
        }
        else
        {
            blas::Trsm
            ( Left, Upper, ConjugateTranspose, NonUnit, (double)1, B, Y );
        }
        X.ResizeTo( k, k );
        X.SetToIdentity();
        blas::Herk( shape, ConjugateTranspose, (double)-1, Y, (double)1, X );
        myResidual = 0; 
        for( int j=0; j<X.LocalWidth(); ++j )
            for( int i=0; i<X.LocalHeight(); ++i )
                myResidual = max( Abs(X.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||Z^H B^-1 Z - I||_oo = " << residual << endl;
    }
}

#ifndef WITHOUT_COMPLEX
void TestCorrectnessDoubleComplex
( bool printMatrices,
  Side side,
  Shape shape,
  const DistMatrix<std::complex<double>,MC,  MR>& A,
  const DistMatrix<std::complex<double>,MC,  MR>& B,
  const DistMatrix<             double, Star,VR>& w,
  const DistMatrix<std::complex<double>,MC,  MR>& Z,
  const DistMatrix<std::complex<double>,MC  ,MR>& ARef,
  const DistMatrix<std::complex<double>,MC  ,MR>& BRef )
{
    const Grid& g = A.GetGrid();
    const int n = Z.Height();
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
    DistMatrix<double,Star,MR> w_Star_MR(true,Z.RowAlignment(),g); 
    w_Star_MR = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( side == Right )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of AZ from BZW...";
            cout.flush();
        }
        // Set X := BZW, where W is the diagonal eigenvalue matrix
        DistMatrix<std::complex<double>,MC,MR> X( g );
        X.AlignWith( Z );
        X.ResizeTo( n, k );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), BRef, Z, 
          std::complex<double>(0), X );
        for( int j=0; j<X.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j);
            elemental::wrappers::blas::Scal
            ( 2*X.LocalHeight(), omega, (double*)X.LocalBuffer(0,j), 1 );
        }
        // X := X - AZ = BZW - AZ
        blas::Hemm
        ( Left, shape, std::complex<double>(-1), ARef, Z, 
        std::complex<double>(1), X );
        // Compute the residual, ||BZW-AZ||_oo
        double myResidual = 0;
        for( int j=0; j<X.LocalWidth(); ++j )
            for( int i=0; i<X.LocalHeight(); ++i )
                myResidual = max( Abs(X.GetLocalEntry(i,j)),myResidual );
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||A Z - B Z W||_oo = " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        DistMatrix<std::complex<double>,MC,MR> Y(g);
        Y = Z;
        if( shape == Lower )
        {
            blas::Trmm
            ( Left, Lower, ConjugateTranspose, NonUnit, 
              std::complex<double>(1), B, Y );
        }
        else
        {
            blas::Trmm
            ( Left, Upper, Normal, NonUnit,
              std::complex<double>(1), B, Y );
        }
        X.ResizeTo( k, k );
        X.SetToIdentity();
        blas::Herk
        ( shape, ConjugateTranspose, std::complex<double>(-1), Y, 
          std::complex<double>(1), X );
        myResidual = 0; 
        for( int j=0; j<X.LocalWidth(); ++j )
            for( int i=0; i<X.LocalHeight(); ++i )
                myResidual = max( Abs(X.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||Z^H B Z - I||_oo = " << residual << endl;
    }
    else
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of BAZ from ZW...";
            cout.flush();
        }
        // Set X := AZ
        DistMatrix<std::complex<double>,MC,MR> X( g );
        X.AlignWith( Z );
        X.ResizeTo( n, k );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), ARef, Z, 
          std::complex<double>(0), X );
        // Set Y := BX = BAZ
        DistMatrix<std::complex<double>,MC,MR> Y( n, k, g );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), BRef, X, 
          std::complex<double>(0), Y );
        // Compute the residual, ||Y - ZW||_oo = ||BAZ - ZW||_oo
        double myResidual = 0;
        for( int j=0; j<Y.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Y.LocalHeight(); ++i )
            {
                double thisResidual = 
                    Abs(omega*Z.GetLocalEntry(i,j)-Y.GetLocalEntry(i,j));
                myResidual = max( thisResidual, myResidual );
            }
        }
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||B A Z - Z W||_oo =  " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        Y = Z;
        if( shape == Lower )
        {
            blas::Trsm
            ( Left, Lower, Normal, NonUnit, std::complex<double>(1), B, Y );
        }
        else
        {
            blas::Trsm
            ( Left, Upper, ConjugateTranspose, NonUnit, 
              std::complex<double>(1), B, Y );
        }
        X.ResizeTo( k, k );
        X.SetToIdentity();
        blas::Herk
        ( shape, ConjugateTranspose, std::complex<double>(-1), Y, 
          std::complex<double>(1), X );
        myResidual = 0; 
        for( int j=0; j<X.LocalWidth(); ++j )
            for( int i=0; i<X.LocalHeight(); ++i )
                myResidual = max( Abs(X.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||Z^H B^-1 Z - I||_oo = " << residual << endl;
    }
}
#endif // WITHOUT_COMPLEX

void TestGeneralizedHermitianEigDouble
( bool testCorrectness, bool printMatrices,
  Side side, Shape shape, int m, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<double,MC,MR> A(m,m,g);
    DistMatrix<double,MC,MR> B(m,m,g);
    DistMatrix<double,MC,MR> ARef(g);
    DistMatrix<double,MC,MR> BRef(g);
    DistMatrix<double,Star,VR> w(g);
    DistMatrix<double,MC,MR> Z(g);

    A.SetToRandomHPD(); 
    B.SetToRandomHPD();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        ARef = A;
        BRef = B;
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
        cout << "  Starting Generalized Hermitian Eigensolver.\n";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::GeneralizedHermitianEig( side, shape, A, B, w, Z );
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
        B.Print("B after eigensolve");
    }
    if( testCorrectness )
    {
        TestCorrectnessDouble
        ( printMatrices, side, shape, A, B, w, Z, ARef, BRef );
    }
}
    
#ifndef WITHOUT_COMPLEX
void TestGeneralizedHermitianEigDoubleComplex
( bool testCorrectness, bool printMatrices,
  Side side, Shape shape, int m, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<std::complex<double>,MC,  MR> A(m,m,g);
    DistMatrix<std::complex<double>,MC,  MR> B(m,m,g);
    DistMatrix<std::complex<double>,MC,  MR> ARef(g);
    DistMatrix<std::complex<double>,MC,  MR> BRef(g);
    DistMatrix<             double, Star,VR> w(g);
    DistMatrix<std::complex<double>,MC,  MR> Z(g);

    A.SetToRandomHPD();
    B.SetToRandomHPD();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        ARef = A;
        BRef = B;
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
        cout << "  Starting Generalized Hermitian Eigensolver.\n";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::GeneralizedHermitianEig( side, shape, A, B, w, Z );
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
        B.Print("B after eigensolve");
    }
    if( testCorrectness )
    {
        TestCorrectnessDoubleComplex
        ( printMatrices, side, shape, A, B, w, Z, ARef, BRef );
    }
}
#endif // WITHOUT_COMPLEX

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
        const int r = atoi(argv[1]);
        const int c = atoi(argv[2]);
        const Side side = CharToSide(*argv[3]);
        const Shape shape = CharToShape(*argv[4]);
        const int m = atoi(argv[5]);
        const int nb = atoi(argv[6]);
        const bool testCorrectness = atoi(argv[7]);
        const bool printMatrices = atoi(argv[8]);
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
            cout << "Will test " << ( side==Left ? "left" : "right" ) << " "
                 << ( shape==Lower ? "lower" : "upper" )
                 << " GeneralizedHermitianEig." << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestGeneralizedHermitianEigDouble
        ( testCorrectness, printMatrices, side, shape, m, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestGeneralizedHermitianEigDoubleComplex
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

