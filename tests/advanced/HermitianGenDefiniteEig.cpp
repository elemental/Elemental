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
            " <range> <a> <b> <highAccuracy?> <shape> <m> <nb> "
            "<local nb symv/hemv> <correctness?> "
            "<print?>\n\n"
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
         << "  highAccuracy? try for high acc. iff != 0\n"
         << "  shape: L/U\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  local nb symv/hemv: local blocksize for symv/hemv\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

void TestCorrectnessDouble
( bool printMatrices,
  HermitianGenDefiniteEigType eigType,
  Shape shape,
  const DistMatrix<double,MC,  MR>& A,
  const DistMatrix<double,MC,  MR>& B,
  const DistMatrix<double,STAR,VR>& w,
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
    DistMatrix<double,STAR,MR> w_STAR_MR(g); 
    w_STAR_MR.AlignWith( X );
    w_STAR_MR = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( eigType == AXBX )
    {
        if( g.VCRank() == 0 )
            cout << "  Testing for deviation of AX from BXW..." << endl;
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm( LEFT, shape, (double)1, BOrig, X, (double)0, Y );
        for( int j=0; j<X.LocalWidth(); ++j )
        {
            double omega = w_STAR_MR.GetLocalEntry(0,j);
            elemental::blas::Scal
            ( Y.LocalHeight(), omega, Y.LocalBuffer(0,j), 1 );
        }
        // Y := Y - AX = BXW - AX
        basic::Hemm( LEFT, shape, (double)-1, AOrig, X, (double)1, Y );
        // Find the infinity norms of A, B, and X, and ||BXW-AX||
        double infNormOfA = 
            advanced::HermitianNorm( shape, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( shape, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( shape, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( shape, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Y, ONE_NORM );
        double infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
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
        if( shape == LOWER )
            basic::Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (double)1, B, Z );
        else
            basic::Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (double)1, B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk( shape, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
        {
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
        }
    }
    else if( eigType == ABX )
    {
        if( g.VCRank() == 0 )
            cout << "  Testing for deviation of ABX from XW..." << endl;
        // Set Y := BX
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm( LEFT, shape, (double)1, BOrig, X, (double)0, Y );
        // Set Z := AY = ABX
        DistMatrix<double,MC,MR> Z( n, k, g );
        basic::Hemm( LEFT, shape, (double)1, AOrig, Y, (double)0, Z );
        // Set Z := Z - XW = ABX - XW
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_STAR_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
                Z.SetLocalEntry(i,j,
                    Z.GetLocalEntry(i,j)-omega*X.GetLocalEntry(i,j));
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        double infNormOfA = 
            advanced::HermitianNorm( shape, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( shape, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( shape, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( shape, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
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
        if( shape == LOWER )
            basic::Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (double)1, B, Z );
        else
            basic::Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (double)1, B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk( shape, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else /* eigType == BAX */
    {
        if( g.VCRank() == 0 )
            cout << "  Testing for deviation of BAX from XW..." << endl;
        // Set Y := AX
        DistMatrix<double,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm( LEFT, shape, (double)1, AOrig, X, (double)0, Y );
        // Set Z := BY = BAX
        DistMatrix<double,MC,MR> Z( n, k, g );
        basic::Hemm( LEFT, shape, (double)1, BOrig, Y, (double)0, Z );
        // Set Z := Z - XW = BAX - XW
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_STAR_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
                Z.SetLocalEntry(i,j,
                    Z.GetLocalEntry(i,j)-omega*X.GetLocalEntry(i,j));
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        double infNormOfA = 
            advanced::HermitianNorm( shape, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( shape, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( shape, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( shape, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
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
        if( shape == LOWER )
            basic::Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (double)1, B, Z );
        else
            basic::Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (double)1, B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk( shape, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
            cout << "    ||X^H B^-1 X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_F  = " << frobNormOfError << endl;
    }
}

#ifndef WITHOUT_COMPLEX
void TestCorrectnessDoubleComplex
( bool printMatrices,
  HermitianGenDefiniteEigType eigType,
  Shape shape,
  const DistMatrix<std::complex<double>,MC,  MR>& A,
  const DistMatrix<std::complex<double>,MC,  MR>& B,
  const DistMatrix<             double, STAR,VR>& w,
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
    DistMatrix<double,STAR,MR> w_STAR_MR(true,X.RowAlignment(),g); 
    w_STAR_MR = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( eigType == AXBX )
    {
        if( g.VCRank() == 0 )
            cout << "  Testing for deviation of AX from BXW..." << endl;
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<std::complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm
        ( LEFT, shape, std::complex<double>(1), BOrig, X, 
          std::complex<double>(0), Y );
        for( int j=0; j<Y.LocalWidth(); ++j )
        {
            double omega = w_STAR_MR.GetLocalEntry(0,j);
            elemental::blas::Scal
            ( 2*Y.LocalHeight(), omega, (double*)Y.LocalBuffer(0,j), 1 );
        }
        // Y := Y - AX = BXW - AX
        basic::Hemm
        ( LEFT, shape, std::complex<double>(-1), AOrig, X, 
        std::complex<double>(1), Y );
        // Find the infinity norms of A, B, X, and AX-BXW
        double infNormOfA = 
            advanced::HermitianNorm( shape, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( shape, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( shape, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( shape, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Y, ONE_NORM );
        double infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
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
        DistMatrix<std::complex<double>,MC,MR> Z(g);
        Z = X;
        if( shape == LOWER )
            basic::Trmm
            ( LEFT, LOWER, ADJOINT, NON_UNIT, std::complex<double>(1), B, Z );
        else
            basic::Trmm
            ( LEFT, UPPER, NORMAL, NON_UNIT, std::complex<double>(1), B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk
        ( shape, ADJOINT, std::complex<double>(-1), Z, 
          std::complex<double>(1), Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else if( eigType == ABX )
    {
        if( g.VCRank() == 0 )
            cout << "  Testing for deviation of ABX from XW..." << endl;
        // Set Y := BX
        DistMatrix<std::complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm
        ( LEFT, shape, std::complex<double>(1), BOrig, X, 
          std::complex<double>(0), Y );
        // Set Z := AY = ABX
        DistMatrix<std::complex<double>,MC,MR> Z( n, k, g );
        basic::Hemm
        ( LEFT, shape, std::complex<double>(1), AOrig, Y, 
          std::complex<double>(0), Z );
        // Set Z := Z - XW = ABX - XW
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_STAR_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
                Z.SetLocalEntry(i,j,
                    Z.GetLocalEntry(i,j)-omega*X.GetLocalEntry(i,j));
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        double infNormOfA = 
            advanced::HermitianNorm( shape, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( shape, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( shape, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( shape, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
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
        if( shape == LOWER )
            basic::Trmm
            ( LEFT, LOWER, ADJOINT, NON_UNIT, std::complex<double>(1), B, Z );
        else
            basic::Trmm
            ( LEFT, UPPER, NORMAL, NON_UNIT, std::complex<double>(1), B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk
        ( shape, ADJOINT, std::complex<double>(-1), Z, 
          std::complex<double>(1), Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else /* eigType == BAX */
    {
        if( g.VCRank() == 0 )
            cout << "  Testing for deviation of BAX from XW..." << endl;
        // Set Y := AX
        DistMatrix<std::complex<double>,MC,MR> Y( g );
        Y.AlignWith( X );
        Y.ResizeTo( n, k );
        basic::Hemm
        ( LEFT, shape, std::complex<double>(1), AOrig, X, 
          std::complex<double>(0), Y );
        // Set Z := BY = BAX
        DistMatrix<std::complex<double>,MC,MR> Z( n, k, g );
        basic::Hemm
        ( LEFT, shape, std::complex<double>(1), BOrig, Y, 
          std::complex<double>(0), Z );
        // Set Z := Z - XW = BAX-XW
        for( int j=0; j<Z.LocalWidth(); ++j )
        {
            double omega = w_STAR_MR.GetLocalEntry(0,j); 
            for( int i=0; i<Z.LocalHeight(); ++i )
                Z.SetLocalEntry(i,j,
                    Z.GetLocalEntry(i,j)-omega*X.GetLocalEntry(i,j));
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        double infNormOfA = 
            advanced::HermitianNorm( shape, AOrig, INFINITY_NORM );
        double frobNormOfA = 
            advanced::HermitianNorm( shape, AOrig, FROBENIUS_NORM );
        double infNormOfB = 
            advanced::HermitianNorm( shape, BOrig, INFINITY_NORM );
        double frobNormOfB = 
            advanced::HermitianNorm( shape, BOrig, FROBENIUS_NORM );
        double oneNormOfX = advanced::Norm( X, ONE_NORM );
        double infNormOfX = advanced::Norm( X, INFINITY_NORM );
        double frobNormOfX = advanced::Norm( X, FROBENIUS_NORM );
        double oneNormOfError = advanced::Norm( Z, ONE_NORM );
        double infNormOfError = advanced::Norm( Z, INFINITY_NORM );
        double frobNormOfError = advanced::Norm( Z, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
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
        if( shape == LOWER )
            basic::Trsm
            ( LEFT, LOWER, NORMAL, NON_UNIT, std::complex<double>(1), B, Z );
        else
            basic::Trsm
            ( LEFT, UPPER, ADJOINT, NON_UNIT, std::complex<double>(1), B, Z );
        Y.ResizeTo( k, k );
        Y.SetToIdentity();
        basic::Herk
        ( shape, ADJOINT, std::complex<double>(-1), Z, 
          std::complex<double>(1), Y );
        oneNormOfError = advanced::Norm( Y, ONE_NORM );
        infNormOfError = advanced::Norm( Y, INFINITY_NORM );
        frobNormOfError = advanced::Norm( Y, FROBENIUS_NORM );
        if( g.VCRank() == 0 )
            cout << "    ||X^H B^-1 X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_F  = " << frobNormOfError << endl;
    }
}
#endif // WITHOUT_COMPLEX

void TestHermitianGenDefiniteEigDouble
( bool testCorrectness, bool printMatrices,
  HermitianGenDefiniteEigType eigType, 
  bool onlyEigenvalues, Shape shape, 
  int m, char range, double vl, double vu, int il, int iu,
  bool tryForHighAccuracy, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<double,MC,MR> A(m,m,g);
    DistMatrix<double,MC,MR> B(m,m,g);
    DistMatrix<double,MC,MR> AOrig(g);
    DistMatrix<double,MC,MR> BOrig(g);
    DistMatrix<double,STAR,VR> w(g);
    DistMatrix<double,MC,MR> X(g);

    A.SetToRandomHPD();
    if( eigType == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        DistMatrix<double,MC,MR> C(m,m,g);
        C.SetToRandom();
        basic::Herk( shape, ADJOINT, (double)1, C, (double)0, B );
    }
    else
        B.SetToRandomHPD();

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
        cout << "  Starting Hermitian Generalized-Definite Eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, tryForHighAccuracy );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, il, iu, tryForHighAccuracy );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, vl, vu, tryForHighAccuracy );
    }
    else
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, X, tryForHighAccuracy );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, X, il, iu, tryForHighAccuracy );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, X, vl, vu, tryForHighAccuracy );
    }
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
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
        ( printMatrices, eigType, shape, A, B, w, X, AOrig, BOrig );
    }
}
    
#ifndef WITHOUT_COMPLEX
void TestHermitianGenDefiniteEigDoubleComplex
( bool testCorrectness, bool printMatrices,
  HermitianGenDefiniteEigType eigType, 
  bool onlyEigenvalues, Shape shape, 
  int m, char range, double vl, double vu, int il, int iu, 
  bool tryForHighAccuracy, const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<std::complex<double>,MC,  MR> A(m,m,g);
    DistMatrix<std::complex<double>,MC,  MR> B(m,m,g);
    DistMatrix<std::complex<double>,MC,  MR> AOrig(g);
    DistMatrix<std::complex<double>,MC,  MR> BOrig(g);
    DistMatrix<             double, STAR,VR> w(g);
    DistMatrix<std::complex<double>,MC,  MR> X(g);

    A.SetToRandomHPD();
    if( eigType == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        DistMatrix<std::complex<double>,MC,MR> C(m,m,g);
        C.SetToRandom();
        basic::Herk
        ( shape, ADJOINT, 
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
        cout << "  Starting Hermitian Generalized-Definite Eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.VCComm() );
    startTime = mpi::Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, tryForHighAccuracy );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, il, iu, tryForHighAccuracy );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, vl, vu, tryForHighAccuracy );
    }
    else
    {
        if( range == 'A' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, X, tryForHighAccuracy );
        else if( range == 'I' )
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, X, il, iu, tryForHighAccuracy );
        else
            advanced::HermitianGenDefiniteEig
            ( eigType, shape, A, B, w, X, vl, vu, tryForHighAccuracy );
    }
    mpi::Barrier( g.VCComm() );
    endTime = mpi::Time();
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
        ( printMatrices, eigType, shape, A, B, w, X, AOrig, BOrig );
    }
}
#endif // WITHOUT_COMPLEX

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::CommRank( comm );

    if( argc < 15 )
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
        const int nbLocalSymv = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);

        if( testCorrectness && onlyEigenvalues && rank==0 )
            cout << "Cannot test correctness with only eigenvalues." << endl;

        HermitianGenDefiniteEigType eigType;
        std::string eigTypeString;
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
            throw std::runtime_error
                  ( "Invalid HermitianGenDefiniteEigType, choose from {1,2,3}" );
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
        basic::SetLocalHemvBlocksize< std::complex<double> >( nbLocalSymv );
#endif

        if( rank == 0 )
        {
            cout << "Will test " 
                 << ( shape==LOWER ? "lower" : "upper" )
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
          eigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu,
          tryForHighAccuracy, g );

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
          eigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu,
          tryForHighAccuracy, g );

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
          eigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu,
          tryForHighAccuracy, g );

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
          eigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu, 
          tryForHighAccuracy, g );

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
          eigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu, 
          tryForHighAccuracy, g );

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
          eigType, onlyEigenvalues, shape, m, range, vl, vu, il, iu, 
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

