/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <ctime>
#include "elemental.hpp"
using namespace std;
using namespace elem;

void TestCorrectness
( bool print,
  HermitianGenDefiniteEigType eigType,
  UpperOrLower uplo,
  const DistMatrix<double>& A,
  const DistMatrix<double>& B,
  const DistMatrix<double,VR,STAR>& w,
  const DistMatrix<double>& X,
  const DistMatrix<double>& AOrig,
  const DistMatrix<double>& BOrig )
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
        DistMatrix<double> Y( g );
        Y.AlignWith( X );
        Zeros( n, k, Y );
        Hemm( LEFT, uplo, (double)1, BOrig, X, (double)0, Y );
        for( int jLocal=0; jLocal<X.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocal(jLocal,0);
            blas::Scal( Y.LocalHeight(), omega, Y.LocalBuffer(0,jLocal), 1 );
        }
        // Y := Y - AX = BXW - AX
        Hemm( LEFT, uplo, (double)-1, AOrig, X, (double)1, Y );
        // Find the infinity norms of A, B, and X, and ||BXW-AX||
        double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = Norm( X, ONE_NORM );
        double infNormOfX = Norm( X, INFINITY_NORM );
        double frobNormOfX = Norm( X, FROBENIUS_NORM );
        double oneNormOfError = Norm( Y, ONE_NORM );
        double infNormOfError = Norm( Y, INFINITY_NORM );
        double frobNormOfError = Norm( Y, FROBENIUS_NORM );
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
        DistMatrix<double> Z(g);
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (double)1, B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (double)1, B, Z );
        Identity( k, k, Y );
        Herk( uplo, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = Norm( Y, ONE_NORM );
        infNormOfError = Norm( Y, INFINITY_NORM );
        frobNormOfError = Norm( Y, FROBENIUS_NORM );
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
        DistMatrix<double> Y( g );
        Y.AlignWith( X );
        Zeros( n, k, Y );
        Hemm( LEFT, uplo, (double)1, BOrig, X, (double)0, Y );
        // Set Z := AY = ABX
        DistMatrix<double> Z( n, k, g );
        Hemm( LEFT, uplo, (double)1, AOrig, Y, (double)0, Z );
        // Set Z := Z - XW = ABX - XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocal(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const double chi = X.GetLocal(iLocal,jLocal);
                const double zeta = Z.GetLocal(iLocal,jLocal);
                Z.SetLocal(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = Norm( X, ONE_NORM );
        double infNormOfX = Norm( X, INFINITY_NORM );
        double frobNormOfX = Norm( X, FROBENIUS_NORM );
        double oneNormOfError = Norm( Z, ONE_NORM );
        double infNormOfError = Norm( Z, INFINITY_NORM );
        double frobNormOfError = Norm( Z, FROBENIUS_NORM );
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
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, (double)1, B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, (double)1, B, Z );
        Identity( k, k, Y );
        Herk( uplo, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = Norm( Y, ONE_NORM );
        infNormOfError = Norm( Y, INFINITY_NORM );
        frobNormOfError = Norm( Y, FROBENIUS_NORM );
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
        DistMatrix<double> Y( g );
        Y.AlignWith( X );
        Zeros( n, k, Y );
        Hemm( LEFT, uplo, (double)1, AOrig, X, (double)0, Y );
        // Set Z := BY = BAX
        DistMatrix<double> Z( n, k, g );
        Hemm( LEFT, uplo, (double)1, BOrig, Y, (double)0, Z );
        // Set Z := Z - XW = BAX - XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocal(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const double chi = X.GetLocal(iLocal,jLocal);
                const double zeta = Z.GetLocal(iLocal,jLocal);
                Z.SetLocal(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = Norm( X, ONE_NORM );
        double infNormOfX = Norm( X, INFINITY_NORM );
        double frobNormOfX = Norm( X, FROBENIUS_NORM );
        double oneNormOfError = Norm( Z, ONE_NORM );
        double infNormOfError = Norm( Z, INFINITY_NORM );
        double frobNormOfError = Norm( Z, FROBENIUS_NORM );
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
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (double)1, B, Z );
        else
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (double)1, B, Z );
        Identity( k, k, Y );
        Herk( uplo, ADJOINT, (double)-1, Z, (double)1, Y );
        oneNormOfError = Norm( Y, ONE_NORM );
        infNormOfError = Norm( Y, INFINITY_NORM );
        frobNormOfError = Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
            cout << "    ||X^H B^-1 X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_F  = " << frobNormOfError << endl;
    }
}

void TestCorrectness
( bool print,
  HermitianGenDefiniteEigType eigType,
  UpperOrLower uplo,
  const DistMatrix<Complex<double> >& A,
  const DistMatrix<Complex<double> >& B,
  const DistMatrix<double,VR,STAR>& w,
  const DistMatrix<Complex<double> >& X,
  const DistMatrix<Complex<double> >& AOrig,
  const DistMatrix<Complex<double> >& BOrig )
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
        DistMatrix<Complex<double> > Y( g );
        Y.AlignWith( X );
        Zeros( n, k, Y );
        Hemm
        ( LEFT, uplo, 
          Complex<double>(1), BOrig, X, 
          Complex<double>(0), Y );
        for( int jLocal=0; jLocal<Y.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocal(jLocal,0);
            blas::Scal
            ( 2*Y.LocalHeight(), omega, (double*)Y.LocalBuffer(0,jLocal), 1 );
        }
        // Y := Y - AX = BXW - AX
        Hemm
        ( LEFT, uplo, 
          Complex<double>(-1), AOrig, X, 
          Complex<double>(1), Y );
        // Find the infinity norms of A, B, X, and AX-BXW
        double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = Norm( X, ONE_NORM );
        double infNormOfX = Norm( X, INFINITY_NORM );
        double frobNormOfX = Norm( X, FROBENIUS_NORM );
        double oneNormOfError = Norm( Y, ONE_NORM );
        double infNormOfError = Norm( Y, INFINITY_NORM );
        double frobNormOfError = Norm( Y, FROBENIUS_NORM );
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
        DistMatrix<Complex<double> > Z(g);
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, Complex<double>(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, Complex<double>(1), B, Z );
        Identity( k, k, Y );
        Herk
        ( uplo, ADJOINT, 
          Complex<double>(-1), Z, 
          Complex<double>(1), Y );
        oneNormOfError = Norm( Y, ONE_NORM );
        infNormOfError = Norm( Y, INFINITY_NORM );
        frobNormOfError = Norm( Y, FROBENIUS_NORM );
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
        DistMatrix<Complex<double> > Y( g );
        Y.AlignWith( X );
        Zeros( n, k, Y );
        Hemm
        ( LEFT, uplo, 
          Complex<double>(1), BOrig, X, 
          Complex<double>(0), Y );
        // Set Z := AY = ABX
        DistMatrix<Complex<double> > Z( n, k, g );
        Hemm
        ( LEFT, uplo, 
          Complex<double>(1), AOrig, Y, 
          Complex<double>(0), Z );
        // Set Z := Z - XW = ABX - XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocal(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const Complex<double> chi = X.GetLocal(iLocal,jLocal);
                const Complex<double> zeta = Z.GetLocal(iLocal,jLocal);
                Z.SetLocal(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = Norm( X, ONE_NORM );
        double infNormOfX = Norm( X, INFINITY_NORM );
        double frobNormOfX = Norm( X, FROBENIUS_NORM );
        double oneNormOfError = Norm( Z, ONE_NORM );
        double infNormOfError = Norm( Z, INFINITY_NORM );
        double frobNormOfError = Norm( Z, FROBENIUS_NORM );
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
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, Complex<double>(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, Complex<double>(1), B, Z );
        Identity( k, k, Y );
        Herk
        ( uplo, ADJOINT, 
          Complex<double>(-1), Z, 
          Complex<double>(1), Y );
        oneNormOfError = Norm( Y, ONE_NORM );
        infNormOfError = Norm( Y, INFINITY_NORM );
        frobNormOfError = Norm( Y, FROBENIUS_NORM );
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
        DistMatrix<Complex<double> > Y( g );
        Y.AlignWith( X );
        Zeros( n, k, Y );
        Hemm
        ( LEFT, uplo, 
          Complex<double>(1), AOrig, X, 
          Complex<double>(0), Y );
        // Set Z := BY = BAX
        DistMatrix<Complex<double> > Z( n, k, g );
        Hemm
        ( LEFT, uplo, 
          Complex<double>(1), BOrig, Y, 
          Complex<double>(0), Z );
        // Set Z := Z - XW = BAX-XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const double omega = w_MR_STAR.GetLocal(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const Complex<double> chi = X.GetLocal(iLocal,jLocal);
                const Complex<double> zeta = Z.GetLocal(iLocal,jLocal);
                Z.SetLocal(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        double infNormOfB = HermitianNorm( uplo, BOrig, INFINITY_NORM );
        double frobNormOfB = HermitianNorm( uplo, BOrig, FROBENIUS_NORM );
        double oneNormOfX = Norm( X, ONE_NORM );
        double infNormOfX = Norm( X, INFINITY_NORM );
        double frobNormOfX = Norm( X, FROBENIUS_NORM );
        double oneNormOfError = Norm( Z, ONE_NORM );
        double infNormOfError = Norm( Z, INFINITY_NORM );
        double frobNormOfError = Norm( Z, FROBENIUS_NORM );
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
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, Complex<double>(1), B, Z );
        else
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, Complex<double>(1), B, Z );
        Identity( k, k, Y );
        Herk
        ( uplo, ADJOINT, 
          Complex<double>(-1), Z, 
          Complex<double>(1), Y );
        oneNormOfError = Norm( Y, ONE_NORM );
        infNormOfError = Norm( Y, INFINITY_NORM );
        frobNormOfError = Norm( Y, FROBENIUS_NORM );
        if( g.Rank() == 0 )
            cout << "    ||X^H B^-1 X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_F  = " << frobNormOfError << endl;
    }
}

void TestHermitianGenDefiniteEigDouble
( bool testCorrectness, bool print,
  HermitianGenDefiniteEigType eigType, 
  bool onlyEigvals, UpperOrLower uplo, 
  int m, char range, double vl, double vu, int il, int iu, const Grid& g )
{
    DistMatrix<double> A(g), AOrig(g);
    DistMatrix<double> B(g), BOrig(g);
    DistMatrix<double,VR,STAR> w(g);
    DistMatrix<double> X(g);

    HermitianUniformSpectrum( m, A, 1, 10 );
    if( eigType == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        Zeros( m, m, B );
        DistMatrix<double> C(g);
        Uniform( m, m, C );
        Herk( uplo, ADJOINT, (double)1, C, (double)0, B );
    }
    else
        HermitianUniformSpectrum( m, B, 1, 10 );

    if( testCorrectness && !onlyEigvals )
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
    if( print )
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
    const double startTime = mpi::Time();
    if( onlyEigvals )
    {
        if( range == 'A' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w );
        else if( range == 'I' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, il, iu );
        else
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, vl, vu );
    }
    else
    {
        if( range == 'A' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X );
        else if( range == 'I' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X, il, iu );
        else
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X, vl, vu );
    }
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << endl;
    }
    if( print )
    {
        w.Print("eigenvalues:");
        if( !onlyEigvals )
            X.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigvals )
    {
        TestCorrectness
        ( print, eigType, uplo, A, B, w, X, AOrig, BOrig );
    }
}
    
void TestHermitianGenDefiniteEigDoubleComplex
( bool testCorrectness, bool print,
  HermitianGenDefiniteEigType eigType, 
  bool onlyEigvals, UpperOrLower uplo, 
  int m, char range, double vl, double vu, int il, int iu, const Grid& g )
{
    DistMatrix<Complex<double> > A(g), AOrig(g);
    DistMatrix<Complex<double> > B(g), BOrig(g);
    DistMatrix<double,VR,STAR> w(g);
    DistMatrix<Complex<double> > X(g);

    HermitianUniformSpectrum( m, A, 1, 10 );
    if( eigType == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        Zeros( m, m, B );
        DistMatrix<Complex<double> > C(g);
        Uniform( m, m, C );
        Herk
        ( uplo, ADJOINT, 
          Complex<double>(1), C, 
          Complex<double>(0), B );
    }
    else
        HermitianUniformSpectrum( m, B, 1, 10 );

    if( testCorrectness && !onlyEigvals )
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
    if( print )
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
    const double startTime = mpi::Time();
    if( onlyEigvals )
    {
        if( range == 'A' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w );
        else if( range == 'I' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, il, iu );
        else
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, vl, vu );
    }
    else
    {
        if( range == 'A' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X );
        else if( range == 'I' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X, il, iu );
        else
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X, vl, vu );
    }
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << endl;
    }
    if( print )
    {
        w.Print("eigenvalues:");
        if( !onlyEigvals )
            X.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigvals )
    {
        TestCorrectness
        ( print, eigType, uplo, A, B, w, X, AOrig, BOrig );
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try
    {
        int r = Input("--gridHeight","height of process grid",0);
        const int eigInt = Input("--eigType",
             "1 is A x = lambda B x, "
             "2 is A B x = lambda x, "
             "3 is B A x = lambda x",1);
        const bool onlyEigvals = Input 
            ("--onlyEigvals","only compute eigenvalues?",false);
        const char range = Input
            ("--range",
             "range of eigenpairs: 'A' for all, 'I' for index range, "
             "'V' for value range",'A');
        const int il = Input("--il","lower bound of index range",0);
        const int iu = Input("--iu","upper bound of index range",100);
        const double vl = Input("--vl","lower bound of value range",0.);
        const double vu = Input("--vu","upper bound of value range",100.);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const int m = Input("--height","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalSymvBlocksize<double>( nbLocal );
        SetLocalHemvBlocksize<Complex<double> >( nbLocal );
        if( range != 'A' && range != 'I' && range != 'V' )
            throw runtime_error("'range' must be 'A', 'I', or 'V'");
        if( testCorrectness && onlyEigvals && commRank==0 )
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
            throw logic_error("Invalid eigenvalue problem integer");
#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test " 
                 << ( uplo==LOWER ? "lower" : "upper" )
                 << " " << eigTypeString << " HermitianGenDefiniteEig." << endl;

        if( commRank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision normal tridiag algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestHermitianGenDefiniteEigDouble
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision square tridiag algorithm, \n"
                 << "row-major grid:\n"
                 << "-------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianGenDefiniteEigDouble
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision square tridiag algorithm, \n"
                 << "col-major grid:\n"
                 << "-------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianGenDefiniteEigDouble
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------------------\n"
                 << "Testing with double-precision complex normal algorithm:\n"
                 << "-------------------------------------------------------" 
                 << endl;
        }
        TestHermitianGenDefiniteEigDoubleComplex
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "row-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianGenDefiniteEigDoubleComplex
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "col-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianGenDefiniteEigDoubleComplex
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );
    }
    catch( ArgException& e ) { }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }   
    Finalize();
    return 0;
}
