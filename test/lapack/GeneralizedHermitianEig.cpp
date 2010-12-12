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
         << "  GeneralizedHermitianEig <r> <c> <genEigType> <only eigenvalues?>"
            " <range> <a> <b> <highAccuracy?> <shape> <m> <nb> <symv local nb> "
            "<hemv local nb> <correctness?> "
            "<print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  genEigType: 1 -> AX=BXW, 2 -> ABX=XW, 3-> BAX=XW\n"
         << "  only eigenvalues?: 0/1\n"
         << "  range: 'A' for all, 'I' for index range, "
            "'V' for floating-point range\n"
         << "  a: if range=='I', 0-indexed first eigenpair to compute\n"
            "     if range=='V', lower-bound on eigenvalues\n"
         << "  b: if range=='I', 0-indexed last eigenpair to compute\n"
            "     if range=='V', upper-bound on eigenvalues\n"
         << "  highAccuracy? try for high acc. iff != 0\n"
         << "  shape: L/U\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  Symv local nb: local blocksize for Symv, double-prec.\n"
         << "  Hemv local nb: \" \", complex double-precision\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

void TestCorrectnessDouble
( bool printMatrices,
  lapack::GenEigType genEigType,
  Shape shape,
  const DistMatrix<double,MC,  MR>& A,
  const DistMatrix<double,MC,  MR>& B,
  const DistMatrix<double,Star,VR>& w,
  const DistMatrix<double,MC,  MR>& X,
  const DistMatrix<double,MC  ,MR>& AOrig,
  const DistMatrix<double,MC,  MR>& BOrig )
{
    const Grid& g = A.Grid();
    const int n = X.Height();
    const int k = X.Width();

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<double,Star,MR> w_Star_MR(g); 
    w_Star_MR.AlignWith( X );
    w_Star_MR = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( genEigType == lapack::AXBX )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of AX from BXW...";
            cout.flush();
        }
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        blas::Hemm( Left, shape, (double)1, BOrig, X, (double)0, Y );
        for( int j=0; j<X.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j);
            elemental::wrappers::blas::Scal
            ( Y.LocalHeight(), omega, Y.LocalBuffer(0,j), 1 );
        }
        // Y := Y - AX = BXW - AX
        blas::Hemm( Left, shape, (double)-1, AOrig, X, (double)1, Y );
        // Compute the residual, ||BXW-AX||_oo
        double myResidual = 0;
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)),myResidual );
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||A X - B X W||_oo = " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        DistMatrix<double,MC,MR> Z(g);
        Z = X;
        if( shape == Lower )
        {
            blas::Trmm
            ( Left, Lower, ConjugateTranspose, NonUnit, (double)1, B, Z );
        }
        else
        {
            blas::Trmm
            ( Left, Upper, Normal, NonUnit, (double)1, B, Z );
        }
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        blas::Herk( shape, ConjugateTranspose, (double)-1, Z, (double)1, Y );
        myResidual = 0; 
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||X^H B X - I||_oo = " << residual << endl;
    }
    else if( genEigType == lapack::ABX )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of ABX from XW...";
            cout.flush();
        }
        // Set Y := BX
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        blas::Hemm( Left, shape, (double)1, BOrig, X, (double)0, Y );
        // Set Z := AY = ABX
        DistMatrix<double,MC,MR> Z( n, k, g );
        blas::Hemm( Left, shape, (double)1, AOrig, Y, (double)0, Z );
        // Compute the residual, ||Z - XW||_oo = ||ABX - XW||_oo
        double myResidual = 0;
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
            {
                double thisResidual = 
                    Abs(omega*X.GetLocalEntry(i,j)-Z.GetLocalEntry(i,j));
                myResidual = max( thisResidual, myResidual );
            }
        }
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||A B X - X W||_oo =  " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        Z = X;
        if( shape == Lower )
        {
            blas::Trmm
            ( Left, Lower, ConjugateTranspose, NonUnit, (double)1, B, Z );
        }
        else
        {
            blas::Trmm
            ( Left, Upper, Normal, NonUnit, (double)1, B, Z );
        }
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        blas::Herk( shape, ConjugateTranspose, (double)-1, Z, (double)1, Y );
        myResidual = 0; 
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||X^H B X - I||_oo = " << residual << endl;
    }
    else /* genEigType == lapack::BAX */
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of BAX from XW...";
            cout.flush();
        }
        // Set Y := AX
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        blas::Hemm( Left, shape, (double)1, AOrig, X, (double)0, Y );
        // Set Z := BY = BAX
        DistMatrix<double,MC,MR> Z( n, k, g );
        blas::Hemm( Left, shape, (double)1, BOrig, Y, (double)0, Z );
        // Compute the residual, ||Z - XW||_oo = ||BAX - XW||_oo
        double myResidual = 0;
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
            {
                double thisResidual = 
                    Abs(omega*X.GetLocalEntry(i,j)-Z.GetLocalEntry(i,j));
                myResidual = max( thisResidual, myResidual );
            }
        }
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||B A X - X W||_oo =  " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B^-1...";
            cout.flush();
        }
        Z = X;
        if( shape == Lower )
        {
            blas::Trsm
            ( Left, Lower, Normal, NonUnit, (double)1, B, Z );
        }
        else
        {
            blas::Trsm
            ( Left, Upper, ConjugateTranspose, NonUnit, (double)1, B, Z );
        }
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        blas::Herk( shape, ConjugateTranspose, (double)-1, Z, (double)1, Y );
        myResidual = 0; 
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||X^H B^-1 X - I||_oo = " << residual << endl;
    }
}

#ifndef WITHOUT_COMPLEX
void TestCorrectnessDoubleComplex
( bool printMatrices,
  lapack::GenEigType genEigType,
  Shape shape,
  const DistMatrix<std::complex<double>,MC,  MR>& A,
  const DistMatrix<std::complex<double>,MC,  MR>& B,
  const DistMatrix<             double, Star,VR>& w,
  const DistMatrix<std::complex<double>,MC,  MR>& X,
  const DistMatrix<std::complex<double>,MC  ,MR>& AOrig,
  const DistMatrix<std::complex<double>,MC  ,MR>& BOrig )
{
    const Grid& g = A.Grid();
    const int n = X.Height();
    const int k = X.Width();

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<double,Star,MR> w_Star_MR(true,X.RowAlignment(),g); 
    w_Star_MR = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( genEigType == lapack::AXBX )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of AX from BXW...";
            cout.flush();
        }
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<std::complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), BOrig, X, 
          std::complex<double>(0), Y );
        for( int j=0; j<Y.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j);
            elemental::wrappers::blas::Scal
            ( 2*Y.LocalHeight(), omega, (double*)Y.LocalBuffer(0,j), 1 );
        }
        // Y := Y - AX = BXW - AX
        blas::Hemm
        ( Left, shape, std::complex<double>(-1), AOrig, X, 
        std::complex<double>(1), Y );
        // Compute the residual, ||BXW-AX||_oo
        double myResidual = 0;
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)),myResidual );
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||A X - B X W||_oo = " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        DistMatrix<std::complex<double>,MC,MR> Z(g);
        Z = X;
        if( shape == Lower )
        {
            blas::Trmm
            ( Left, Lower, ConjugateTranspose, NonUnit, 
              std::complex<double>(1), B, Z );
        }
        else
        {
            blas::Trmm
            ( Left, Upper, Normal, NonUnit,
              std::complex<double>(1), B, Z );
        }
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        blas::Herk
        ( shape, ConjugateTranspose, std::complex<double>(-1), Z, 
          std::complex<double>(1), Y );
        myResidual = 0; 
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||X^H B X - I||_oo = " << residual << endl;
    }
    else if( genEigType == lapack::ABX )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of ABX from XW...";
            cout.flush();
        }
        // Set Y := BX
        DistMatrix<std::complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), BOrig, X, 
          std::complex<double>(0), Y );
        // Set Z := AY = ABX
        DistMatrix<std::complex<double>,MC,MR> Z( n, k, g );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), AOrig, Y, 
          std::complex<double>(0), Z );
        // Compute the residual, ||Z - XW||_oo = ||ABX - XW||_oo
        double myResidual = 0;
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
            {
                double thisResidual = 
                    Abs(omega*X.GetLocalEntry(i,j)-Z.GetLocalEntry(i,j));
                myResidual = max( thisResidual, myResidual );
            }
        }
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||A B X - X W||_oo =  " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B...";
            cout.flush();
        }
        Z = X;
        if( shape == Lower )
        {
            blas::Trmm
            ( Left, Lower, ConjugateTranspose, NonUnit, 
              std::complex<double>(1), B, Z );
        }
        else
        {
            blas::Trmm
            ( Left, Upper, Normal, NonUnit, 
              std::complex<double>(1), B, Z );
        }
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        blas::Herk
        ( shape, ConjugateTranspose, std::complex<double>(-1), Z, 
          std::complex<double>(1), Y );
        myResidual = 0; 
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||X^H B X - I||_oo = " << residual << endl;
    }
    else /* genEigType == lapack::BAX */
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Testing for deviation of BAX from XW...";
            cout.flush();
        }
        // Set Y := AX
        DistMatrix<std::complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), AOrig, X, 
          std::complex<double>(0), Y );
        // Set Z := BY = BAX
        DistMatrix<std::complex<double>,MC,MR> Z( n, k, g );
        blas::Hemm
        ( Left, shape, std::complex<double>(1), BOrig, Y, 
          std::complex<double>(0), Z );
        // Compute the residual, ||Z - XW||_oo = ||BAX - XW||_oo
        double myResidual = 0;
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_Star_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
            {
                double thisResidual = 
                    Abs(omega*X.GetLocalEntry(i,j)-Z.GetLocalEntry(i,j));
                myResidual = max( thisResidual, myResidual );
            }
        }
        double residual;
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||B A X - X W||_oo =  " << residual << endl;
        
        if( g.VCRank() == 0 )
        {
            cout << "  Testing orthonormality of eigenvectors w.r.t. B^-1...";
            cout.flush();
        }
        Z = X;
        if( shape == Lower )
        {
            blas::Trsm
            ( Left, Lower, Normal, NonUnit, std::complex<double>(1), B, Z );
        }
        else
        {
            blas::Trsm
            ( Left, Upper, ConjugateTranspose, NonUnit, 
              std::complex<double>(1), B, Z );
        }
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        blas::Herk
        ( shape, ConjugateTranspose, std::complex<double>(-1), Z, 
          std::complex<double>(1), Y );
        myResidual = 0; 
        for( int j=0; j<Y.LocalWidth(); ++j )
            for( int i=0; i<Y.LocalHeight(); ++i )
                myResidual = max( Abs(Y.GetLocalEntry(i,j)), myResidual );
        Reduce( &myResidual, &residual, 1, MPI_MAX, 0, g.VCComm() );
        if( g.VCRank() == 0 )
            cout << "||X^H B^-1 X - I||_oo = " << residual << endl;
    }
}
#endif // WITHOUT_COMPLEX

void TestGeneralizedHermitianEigDouble
( bool testCorrectness, bool printMatrices,
  lapack::GenEigType genEigType, bool onlyEigenvalues, Shape shape, 
  int m, char range, double vl, double vu, int il, int iu,
  bool tryForHighAccuracy, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<double,MC,MR> A(m,m,g);
    DistMatrix<double,MC,MR> B(m,m,g);
    DistMatrix<double,MC,MR> AOrig(g);
    DistMatrix<double,MC,MR> BOrig(g);
    DistMatrix<double,Star,VR> w(g);
    DistMatrix<double,MC,MR> X(g);

    A.SetToRandomHPD();
    if( genEigType == lapack::BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        DistMatrix<double,MC,MR> C(m,m,g);
        C.SetToRandom();
        blas::Herk( shape, ConjugateTranspose, (double)1, C, (double)0, B );
    }
    else
    {
        B.SetToRandomHPD();
    }

    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        AOrig = A;
        BOrig = B;
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
        cout << "  Starting Generalized Hermitian Eigensolver...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, tryForHighAccuracy );
        }
        else if( range == 'I' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, il, iu, tryForHighAccuracy );
        }
        else
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, vl, vu, tryForHighAccuracy );
        }
    }
    else
    {
        if( range == 'A' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, X, tryForHighAccuracy );
        }
        else if( range == 'I' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, X, il, iu, tryForHighAccuracy );
        }
        else
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, X, vl, vu, tryForHighAccuracy );
        }
    }
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
        w.Print("eigenvalues:");
        if( !onlyEigenvalues )
            X.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigenvalues )
    {
        TestCorrectnessDouble
        ( printMatrices, genEigType, shape, A, B, w, X, AOrig, BOrig );
    }
}
    
#ifndef WITHOUT_COMPLEX
void TestGeneralizedHermitianEigDoubleComplex
( bool testCorrectness, bool printMatrices,
  lapack::GenEigType genEigType, bool onlyEigenvalues, Shape shape, 
  int m, char range, double vl, double vu, int il, int iu, 
  bool tryForHighAccuracy, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<std::complex<double>,MC,  MR> A(m,m,g);
    DistMatrix<std::complex<double>,MC,  MR> B(m,m,g);
    DistMatrix<std::complex<double>,MC,  MR> AOrig(g);
    DistMatrix<std::complex<double>,MC,  MR> BOrig(g);
    DistMatrix<             double, Star,VR> w(g);
    DistMatrix<std::complex<double>,MC,  MR> X(g);

    A.SetToRandomHPD();
    if( genEigType == lapack::BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        DistMatrix<std::complex<double>,MC,MR> C(m,m,g);
        C.SetToRandom();
        blas::Herk
        ( shape, ConjugateTranspose, 
          std::complex<double>(1), C, std::complex<double>(0), B );
    }
    else
    {
        B.SetToRandomHPD();
    }

    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        AOrig = A;
        BOrig = B;
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
        cout << "  Starting Generalized Hermitian Eigensolver...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, tryForHighAccuracy );
        }
        else if( range == 'I' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, il, iu, tryForHighAccuracy );
        }
        else
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, vl, vu, tryForHighAccuracy );
        }
    }
    else
    {
        if( range == 'A' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, X, tryForHighAccuracy );
        }
        else if( range == 'I' )
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, X, il, iu, tryForHighAccuracy );
        }
        else
        {
            lapack::GeneralizedHermitianEig
            ( genEigType, shape, A, B, w, X, vl, vu, tryForHighAccuracy );
        }
    }
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
        w.Print("eigenvalues:");
        if( !onlyEigenvalues )
            X.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigenvalues )
    {
        TestCorrectnessDoubleComplex
        ( printMatrices, genEigType, shape, A, B, w, X, AOrig, BOrig );
    }
}
#endif // WITHOUT_COMPLEX

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 16 )
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
        const int genEigInt = atoi(argv[++argNum]);
        const bool onlyEigenvalues = atoi(argv[++argNum]);
        const char range = *argv[++argNum];
        if( range != 'A' && range != 'I' && range != 'V' )
            throw std::runtime_error("'range' must be 'A', 'I', or 'V'");
        double vl = 0, vu = 0;
        int il = 0, iu = 0;
        if( range == 'I' )
        {
            il = atoi(argv[++argNum]);
            iu = atoi(argv[++argNum]);
        }
        else if( range == 'V' )
        {
            vl = atof(argv[++argNum]);
            vu = atof(argv[++argNum]);
        }
        else
        {
            argNum += 2;
        }
        const bool tryForHighAccuracy = atoi(argv[++argNum]);
        const Shape shape = CharToShape(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const int nbLocalSymvDouble = atoi(argv[++argNum]);
#ifndef WITHOUT_COMPLEX
        const int nbLocalHemvComplexDouble = atoi(argv[++argNum]);
#else
        ++argNum;
#endif
        const bool testCorrectness = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);

        if( testCorrectness && onlyEigenvalues && rank==0 )
            cout << "Cannot test correctness with only eigenvalues." << endl;

        lapack::GenEigType genEigType;
        std::string genEigString;
        if( genEigInt == 1 )
        {
            genEigType = lapack::AXBX;
            genEigString = "AXBX";
        }
        else if( genEigInt == 2 )
        {
            genEigType = lapack::ABX;
            genEigString = "ABX";
        }
        else if( genEigInt == 3 )
        {
            genEigType = lapack::BAX;
            genEigString = "BAX";
        }
        else
            throw std::runtime_error
                  ( "Invalid GenEigType, choose from {1,2,3}" );
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
        blas::SetLocalSymvDoubleBlocksize( nbLocalSymvDouble );
#ifndef WITHOUT_COMPLEX
        blas::SetLocalHemvComplexDoubleBlocksize( nbLocalHemvComplexDouble );
#endif

        if( rank == 0 )
        {
            cout << "Will test " 
                 << ( shape==Lower ? "lower" : "upper" )
                 << " " << genEigString << " GeneralizedHermitianEig." << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestGeneralizedHermitianEigDouble
        ( testCorrectness, printMatrices, 
          genEigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu,
          tryForHighAccuracy, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestGeneralizedHermitianEigDoubleComplex
        ( testCorrectness, printMatrices, 
          genEigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu, 
          tryForHighAccuracy, g );
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

