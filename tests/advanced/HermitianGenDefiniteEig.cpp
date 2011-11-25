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
    cout << "Generates random Hermitian A and random HPD B then solves for "
         << "their eigenpairs.\n\n"
         << "  HermitianGenDefiniteEig <r> <c> <eigType> <only eigenvalues?>"
            " <range> <a> <b> <uplo> <m> <nb> <local nb symv/hemv> "
            "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  eigType: 1 -> AX=BXW, 2 -> ABX=XW, 3-> BAX=XW\n"
         << "  only eigenvalues?: 0/1\n"
         << "  range: 'A' for all, 'I' for index range, "
            "'V' for floating-point range\n"
         << "  a: if range=='I', 0-indexed first eigenpair to compute\n"
            "     if range=='V', lower-bound on eigenvalues\n"
         << "  b: if range=='I', 0-indexed last eigenpair to compute\n"
            "     if range=='V', upper-bound on eigenvalues\n"
         << "  uplo: L/U\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  local nb symv/hemv: local blocksize for symv/hemv\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

void TestCorrectnessDouble
( bool printMatrices,
  HermitianGenDefiniteEigType eigType,
  UpperOrLower uplo,
  const DistMatrix<double,MC,MR  >& A,
  const DistMatrix<double,MC,MR  >& B,
  const DistMatrix<double,VR,STAR>& w,
  const DistMatrix<double,MC,MR  >& X,
  const DistMatrix<double,MC,MR  >& AOrig,
  const DistMatrix<double,MC,MR  >& BOrig )
{
    const Grid& g = A.Grid();
    const int n = X.Height();
    const int k = X.Width();

    if( g.Rank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<double,MR,STAR> w_MR_STAR(g); 
    w_MR_STAR.AlignWith( X );
    w_MR_STAR = w;
    if( g.Rank() == 0 )
        cout << "DONE" << endl;

    if( eigType == AXBX )
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of AX from BXW..." << endl;
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm( LEFT, uplo, (double)1, BOrig, X, (double)0, Y );
        for( int jLocal=0; jLocal<X.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocalEntry(jLocal,0);
            elemental::blas::Scal
            ( Y.LocalHeight(), omega, Y.LocalBuffer(0,jLocal), 1 );
        }
        // Y := Y - AX = BXW - AX
        basic::Hemm( LEFT, uplo, (double)-1, AOrig, X, (double)1, Y );
        // Find the infinity norms of A, B, and X, and ||BXW-AX||
        double infNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Y, ONE_NORM );
        double infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"       
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||A B X - X W||_1  = " << oneNormOfError << "\n"
                 << "    ||A B X - X W||_oo = " << infNormOfError << "\n"
                 << "    ||A B X - X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B..." 
                 << endl; 
        }
        DistMatrix<double,MC,MR> Z(g);
        Z = X;
        if( uplo == LOWER )
            basic::Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (double)1, B, Z );
        else
            basic::Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (double)1, B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk( uplo, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
        }
    }
    else if( eigType == ABX )
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of ABX from XW..." << endl;
        // Set Y := BX
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm( LEFT, uplo, (double)1, BOrig, X, (double)0, Y );
        // Set Z := AY = ABX
        DistMatrix<double,MC,MR> Z( n, k, g );
        basic::Hemm( LEFT, uplo, (double)1, AOrig, Y, (double)0, Z );
        // Set Z := Z - XW = ABX - XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocalEntry(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const double chi = X.GetLocalEntry(iLocal,jLocal);
                const double zeta = Z.GetLocalEntry(iLocal,jLocal);
                Z.SetLocalEntry(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        double infNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||A B X - X W||_1  = " << oneNormOfError << "\n"
                 << "    ||A B X - X W||_oo = " << infNormOfError << "\n"
                 << "    ||A B X - X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B..." 
                 << endl;
        }
        Z = X;
        if( uplo == LOWER )
            basic::Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (double)1, B, Z );
        else
            basic::Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (double)1, B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk( uplo, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else /* eigType == BAX */
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of BAX from XW..." << endl;
        // Set Y := AX
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm( LEFT, uplo, (double)1, AOrig, X, (double)0, Y );
        // Set Z := BY = BAX
        DistMatrix<double,MC,MR> Z( n, k, g );
        basic::Hemm( LEFT, uplo, (double)1, BOrig, Y, (double)0, Z );
        // Set Z := Z - XW = BAX - XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocalEntry(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const double chi = X.GetLocalEntry(iLocal,jLocal);
                const double zeta = Z.GetLocalEntry(iLocal,jLocal);
                Z.SetLocalEntry(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        double infNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||B A X - X W||_1  = " << oneNormOfError << "\n"
                 << "    ||B A X - X W||_oo = " << infNormOfError << "\n"
                 << "    ||B A X - X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B^-1..."
                 << endl;
        }
        Z = X;
        if( uplo == LOWER )
            basic::Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (double)1, B, Z );
        else
            basic::Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (double)1, B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk( uplo, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
            cout << "    ||X^H B^-1 X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_F  = " << frobNormOfError << endl;
    }
}

#ifndef WITHOUT_COMPLEX
void TestCorrectnessDoubleComplex
( bool printMatrices,
  HermitianGenDefiniteEigType eigType,
  UpperOrLower uplo,
  const DistMatrix<complex<double>,MC,MR  >& A,
  const DistMatrix<complex<double>,MC,MR  >& B,
  const DistMatrix<        double, VR,STAR>& w,
  const DistMatrix<complex<double>,MC,MR  >& X,
  const DistMatrix<complex<double>,MC,MR  >& AOrig,
  const DistMatrix<complex<double>,MC,MR  >& BOrig )
{
    const Grid& g = A.Grid();
    const int n = X.Height();
    const int k = X.Width();

    if( g.Rank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<double,MR,STAR> w_MR_STAR(true,X.RowAlignment(),g); 
    w_MR_STAR = w;
    if( g.Rank() == 0 )
        cout << "DONE" << endl;

    if( eigType == AXBX )
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of AX from BXW..." << endl;
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm
        ( LEFT, uplo, 
          complex<double>(1), BOrig, X, 
          complex<double>(0), Y );
        for( int jLocal=0; jLocal<Y.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocalEntry(jLocal,0);
            elemental::blas::Scal
            ( 2*Y.LocalHeight(), omega, (double*)Y.LocalBuffer(0,jLocal), 1 );
        }
        // Y := Y - AX = BXW - AX
        basic::Hemm
        ( LEFT, uplo, 
          complex<double>(-1), AOrig, X, 
          complex<double>(1), Y );
        // Find the infinity norms of A, B, X, and AX-BXW
        double infNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Y, ONE_NORM );
        double infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||A X - B X W||_1  = " << oneNormOfError << "\n"
                 << "    ||A X - B X W||_oo = " << infNormOfError << "\n"
                 << "    ||A X - B X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B..."
                 << endl;
        }
        DistMatrix<complex<double>,MC,MR> Z(g);
        Z = X;
        if( uplo == LOWER )
            basic::Trmm
            ( LEFT, LOWER, ADJOINT, NON_UNIT, complex<double>(1), B, Z );
        else
            basic::Trmm
            ( LEFT, UPPER, NORMAL, NON_UNIT, complex<double>(1), B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk
        ( uplo, ADJOINT, 
          complex<double>(-1), Z, 
          complex<double>(1), Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else if( eigType == ABX )
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of ABX from XW..." << endl;
        // Set Y := BX
        DistMatrix<complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm
        ( LEFT, uplo, 
          complex<double>(1), BOrig, X, 
          complex<double>(0), Y );
        // Set Z := AY = ABX
        DistMatrix<complex<double>,MC,MR> Z( n, k, g );
        basic::Hemm
        ( LEFT, uplo, 
          complex<double>(1), AOrig, Y, 
          complex<double>(0), Z );
        // Set Z := Z - XW = ABX - XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocalEntry(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const complex<double> chi = X.GetLocalEntry(iLocal,jLocal);
                const complex<double> zeta = 
                    Z.GetLocalEntry(iLocal,jLocal);
                Z.SetLocalEntry(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        double infNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||A B X - X W||_1  = " << oneNormOfError << "\n"
                 << "    ||A B X - X W||_oo = " << infNormOfError << "\n"
                 << "    ||A B X - X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B..."
                 << endl;
        }
        Z = X;
        if( uplo == LOWER )
            basic::Trmm
            ( LEFT, LOWER, ADJOINT, NON_UNIT, complex<double>(1), B, Z );
        else
            basic::Trmm
            ( LEFT, UPPER, NORMAL, NON_UNIT, complex<double>(1), B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk
        ( uplo, ADJOINT, 
          complex<double>(-1), Z, 
          complex<double>(1), Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else /* eigType == BAX */
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of BAX from XW..." << endl;
        // Set Y := AX
        DistMatrix<complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm
        ( LEFT, uplo, 
          complex<double>(1), AOrig, X, 
          complex<double>(0), Y );
        // Set Z := BY = BAX
        DistMatrix<complex<double>,MC,MR> Z( n, k, g );
        basic::Hemm
        ( LEFT, uplo, 
          complex<double>(1), BOrig, Y, 
          complex<double>(0), Z );
        // Set Z := Z - XW = BAX-XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocalEntry(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const complex<double> chi = X.GetLocalEntry(iLocal,jLocal);
                const complex<double> zeta = Z.GetLocalEntry(iLocal,jLocal);
                Z.SetLocalEntry(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        double infNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||B A X - X W||_1  = " << oneNormOfError << "\n"
                 << "    ||B A X - X W||_oo = " << infNormOfError << "\n"
                 << "    ||B A X - X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B^-1..."
                 << endl;
        }
        Z = X;
        if( uplo == LOWER )
            basic::Trsm
            ( LEFT, LOWER, NORMAL, NON_UNIT, complex<double>(1), B, Z );
        else
            basic::Trsm
            ( LEFT, UPPER, ADJOINT, NON_UNIT, complex<double>(1), B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk
        ( uplo, ADJOINT, 
          complex<double>(-1), Z, 
          complex<double>(1), Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
            cout << "    ||X^H B^-1 X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_F  = " << frobNormOfError << endl;
    }
}
#endif // WITHOUT_COMPLEX

void TestHermitianGenDefiniteEigDouble
( bool testCorrectness, bool printMatrices,
  HermitianGenDefiniteEigType eigType, 
  bool onlyEigenvalues, UpperOrLower uplo, 
  int m, char range, double vl, double vu, int il, int iu, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<double,MC,MR  > A(m,m,g);
    DistMatrix<double,MC,MR  > B(m,m,g);
    DistMatrix<double,MC,MR  > AOrig(g);
    DistMatrix<double,MC,MR  > BOrig(g);
    DistMatrix<double,VR,STAR> w(g);
    DistMatrix<double,MC,MR  > X(g);

    A.SetToRandomHPD();
    if( eigType == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        DistMatrix<double,MC,MR> C(m,m,g);
        C.SetToRandom();
        basic::Herk( uplo, ADJOINT, (double)1, C, (double)0, B );
    }
    else
        B.SetToRandomHPD();

    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        AOrig = A;
        BOrig = B;
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
        cout << "  Starting Hermitian Generalized-Definite Eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig( eigType, uplo, A, B, w );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, il, iu );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, vl, vu );
    }
    else
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, X );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, X, il, iu );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, X, vl, vu );
    }
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    if( g.Rank() == 0 )
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
        ( printMatrices, eigType, uplo, A, B, w, X, AOrig, BOrig );
    }
}
    
#ifndef WITHOUT_COMPLEX
void TestHermitianGenDefiniteEigDoubleComplex
( bool testCorrectness, bool printMatrices,
  HermitianGenDefiniteEigType eigType, 
  bool onlyEigenvalues, UpperOrLower uplo, 
  int m, char range, double vl, double vu, int il, int iu, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<complex<double>,MC,MR  > A(m,m,g);
    DistMatrix<complex<double>,MC,MR  > B(m,m,g);
    DistMatrix<complex<double>,MC,MR  > AOrig(g);
    DistMatrix<complex<double>,MC,MR  > BOrig(g);
    DistMatrix<        double, VR,STAR> w(g);
    DistMatrix<complex<double>,MC,MR  > X(g);

    A.SetToRandomHPD();
    if( eigType == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        DistMatrix<complex<double>,MC,MR> C(m,m,g);
        C.SetToRandom();
        basic::Herk
        ( uplo, ADJOINT, 
          complex<double>(1), C, 
          complex<double>(0), B );
    }
    else
    {
        B.SetToRandomHPD();
    }

    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        AOrig = A;
        BOrig = B;
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
        cout << "  Starting Hermitian Generalized-Definite Eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, il, iu );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, vl, vu );
    }
    else
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, X );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, X, il, iu );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, uplo, A, B, w, X, vl, vu );
    }
    mpi::Barrier( g.Comm() );
    endTime = mpi::Time();
    runTime = endTime - startTime;
    if( g.Rank() == 0 )
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
        ( printMatrices, eigType, uplo, A, B, w, X, AOrig, BOrig );
    }
}
#endif // WITHOUT_COMPLEX

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 14 )
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
        const int eigInt = atoi(argv[++argNum]);
        const bool onlyEigenvalues = atoi(argv[++argNum]);
        const char range = *argv[++argNum];
        if( range != 'A' && range != 'I' && range != 'V' )
            throw runtime_error("'range' must be 'A', 'I', or 'V'");
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
        const UpperOrLower uplo = CharToUpperOrLower(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const int nbLocalSymv = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);

        if( testCorrectness && onlyEigenvalues && rank==0 )
            cout << "Cannot test correctness with only eigenvalues." << endl;

        HermitianGenDefiniteEigType eigType;
        string eigTypeString;
        if( eigInt == 1 )
        {
            eigType = AXBX;
            eigTypeString = "AXBX";
        }
        else if( eigInt == 2 )
        {
            eigType = ABX;
            eigTypeString = "ABX";
        }
        else if( eigInt == 3 )
        {
            eigType = BAX;
            eigTypeString = "BAX";
        }
        else
            throw runtime_error
            ("Invalid HermitianGenDefiniteEigType, choose from {1,2,3}");
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
        basic::SetLocalSymvBlocksize<double>( nbLocalSymv );
#ifndef WITHOUT_COMPLEX
        basic::SetLocalHemvBlocksize<complex<double> >( nbLocalSymv );
#endif

        if( rank == 0 )
        {
            cout << "Will test " 
                 << ( uplo==LOWER ? "lower" : "upper" )
                 << " " << eigTypeString << " HermitianGenDefiniteEig." << endl;
        }

        if( rank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision normal tridiag algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        advanced::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestHermitianGenDefiniteEigDouble
        ( testCorrectness, printMatrices, 
          eigType, onlyEigenvalues, uplo, m, range, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision square tridiag algorithm, \n"
                 << "row-major grid:\n"
                 << "-------------------------------------------"
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        advanced::SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianGenDefiniteEigDouble
        ( testCorrectness, printMatrices, 
          eigType, onlyEigenvalues, uplo, m, range, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision square tridiag algorithm, \n"
                 << "col-major grid:\n"
                 << "-------------------------------------------"
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        advanced::SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianGenDefiniteEigDouble
        ( testCorrectness, printMatrices, 
          eigType, onlyEigenvalues, uplo, m, range, vl, vu, il, iu, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "-------------------------------------------------------\n"
                 << "Testing with double-precision complex normal algorithm:\n"
                 << "-------------------------------------------------------" 
                 << endl;
        }
        TestHermitianGenDefiniteEigDoubleComplex
        ( testCorrectness, printMatrices, 
          eigType, onlyEigenvalues, uplo, m, range, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "row-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        advanced::SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianGenDefiniteEigDoubleComplex
        ( testCorrectness, printMatrices, 
          eigType, onlyEigenvalues, uplo, m, range, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "col-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        advanced::SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        advanced::SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianGenDefiniteEigDoubleComplex
        ( testCorrectness, printMatrices, 
          eigType, onlyEigenvalues, uplo, m, range, vl, vu, il, iu, g );
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

