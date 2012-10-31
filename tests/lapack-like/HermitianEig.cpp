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
    cout << "Generates random Hermitian matrix then solves for its eigenpairs."
         << "\n\n"
         << "  HermitianEig <r> <c> <only eigenvalues?> <range> <a> <b> "
            "<clustered?> <uplo> <m> <nb> <local nb symv/hemv> "
            "<correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  only eigenvalues?: 0/1\n"
         << "  range: 'A' for all, 'I' for index range, "
            "'V' for floating-point range\n"
         << "  a: if range=='I', 0-indexed first eigenpair to compute\n"
            "     if range=='V', lower-bound on eigenvalues\n"
         << "  b: if range=='I', 0-indexed last eigenpair to compute\n"
            "     if range=='V', upper-bound on eigenvalues\n"
         << "  clustered: uniform spectrum if 0, clustered spectrum otherwise\n"
         << "  uplo: 'L' for lower-triangular storage, 'U' for upper\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  local nb symv/hemv: local blocksize for symv/hemv\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

void TestCorrectness
( bool printMatrices,
  UpperOrLower uplo,
  const DistMatrix<double>& A,
  const DistMatrix<double,VR,STAR>& w,
  const DistMatrix<double>& Z,
  const DistMatrix<double>& AOrig )
{
    const Grid& g = A.Grid();
    const int n = Z.Height();
    const int k = Z.Width();

    if( g.Rank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<double,MR,STAR> w_MR_STAR(g); 
    w_MR_STAR.AlignWith( Z );
    w_MR_STAR = w;
    if( g.Rank() == 0 )
        cout << "DONE" << endl;

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of eigenvectors..." << endl;
    DistMatrix<double> X(g);
    Identity( k, k, X );
    Herk( uplo, ADJOINT, (double)-1, Z, (double)1, X );
    double oneNormOfError = Norm( X, ONE_NORM );
    double infNormOfError = Norm( X, INFINITY_NORM );
    double frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||Z^H Z - I||_1  = " << oneNormOfError << "\n"
             << "    ||Z^H Z - I||_oo = " << infNormOfError << "\n"
             << "    ||Z^H Z - I||_F  = " << frobNormOfError << "\n\n"
             << "  Testing for deviation of AZ from ZW..." << endl;
    }
    // Set X := AZ
    X.AlignWith( Z );
    Zeros( n, k, X );
    Hemm( LEFT, uplo, (double)1, AOrig, Z, (double)0, X );
    // Set X := X - ZW = AZ - ZW
    for( int jLocal=0; jLocal<X.LocalWidth(); ++jLocal )
    {
        const double omega = w_MR_STAR.GetLocal(jLocal,0);
        for( int iLocal=0; iLocal<X.LocalHeight(); ++iLocal )
        {
            const double chi = X.GetLocal(iLocal,jLocal);
            const double zeta = Z.GetLocal(iLocal,jLocal);
            X.SetLocal(iLocal,jLocal,chi-omega*zeta);
        }
    }
    // Find the infinity norms of A, Z, and AZ-ZW
    double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
    double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
    double oneNormOfZ = Norm( Z, ONE_NORM );
    double infNormOfZ = Norm( Z, INFINITY_NORM );
    double frobNormOfZ = Norm( Z, FROBENIUS_NORM );
    oneNormOfError = Norm( X, ONE_NORM );
    infNormOfError = Norm( X, INFINITY_NORM );
    frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
             << "    ||A||_F            = " << frobNormOfA << "\n"
             << "    ||Z||_1            = " << oneNormOfZ << "\n"
             << "    ||Z||_oo           = " << infNormOfZ << "\n"
             << "    ||Z||_F            = " << frobNormOfZ << "\n"
             << "    ||A Z - Z W||_1    = " << oneNormOfError << "\n"
             << "    ||A Z - Z W||_oo   = " << infNormOfError << "\n"
             << "    ||A Z - Z W||_F    = " << frobNormOfError << endl;
    }
}

void TestCorrectness
( bool printMatrices,
  UpperOrLower uplo,
  const DistMatrix<Complex<double> >& A,
  const DistMatrix<double,VR,STAR>& w,
  const DistMatrix<Complex<double> >& Z,
  const DistMatrix<Complex<double> >& AOrig )
{
    const Grid& g = A.Grid();
    const int n = Z.Height();
    const int k = Z.Width();

    if( g.Rank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<double,MR,STAR> w_MR_STAR(true,Z.RowAlignment(),g); 
    w_MR_STAR = w;
    if( g.Rank() == 0 )
        cout << "DONE" << endl;

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of eigenvectors..." << endl;
    DistMatrix<Complex<double> > X( g );
    Identity( k, k, X );
    Herk
    ( uplo, ADJOINT, 
      Complex<double>(-1), Z, 
      Complex<double>(1), X );
    double oneNormOfError = Norm( X, ONE_NORM );
    double infNormOfError = Norm( X, INFINITY_NORM );
    double frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||Z^H Z - I||_1  = " << oneNormOfError << "\n"
             << "    ||Z^H Z - I||_oo = " << infNormOfError << "\n"
             << "    ||Z^H Z - I||_F  = " << frobNormOfError << "\n\n"
             << "  Testing for deviation of AZ from ZW..." << endl;
    }
    // X := AZ
    X.AlignWith( Z );
    Zeros( n, k, X );
    Hemm
    ( LEFT, uplo, 
      Complex<double>(1), AOrig, Z, 
      Complex<double>(0), X );
    // Find the residual ||X-ZW||_oo = ||AZ-ZW||_oo
    for( int jLocal=0; jLocal<X.LocalWidth(); ++jLocal )
    {
        const double omega = w_MR_STAR.GetLocal(jLocal,0);
        for( int iLocal=0; iLocal<X.LocalHeight(); ++iLocal )
        {
            const Complex<double> chi = X.GetLocal(iLocal,jLocal);
            const Complex<double> zeta = Z.GetLocal(iLocal,jLocal);
            X.SetLocal(iLocal,jLocal,chi-omega*zeta);
        }
    }
    // Find the infinity norms of A, Z, and AZ-ZW
    double infNormOfA = HermitianNorm( uplo, AOrig, INFINITY_NORM );
    double frobNormOfA = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
    double oneNormOfZ = Norm( Z, ONE_NORM );
    double infNormOfZ = Norm( Z, INFINITY_NORM );
    double frobNormOfZ = Norm( Z, FROBENIUS_NORM );
    oneNormOfError = Norm( X, ONE_NORM );
    infNormOfError = Norm( X, INFINITY_NORM );
    frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
             << "    ||A||_F            = " << frobNormOfA << "\n"
             << "    ||Z||_1            = " << oneNormOfZ << "\n"
             << "    ||Z||_oo           = " << infNormOfZ << "\n"
             << "    ||Z||_F            = " << frobNormOfZ << "\n"
             << "    ||A Z - Z W||_1    = " << oneNormOfError << "\n"
             << "    ||A Z - Z W||_oo   = " << infNormOfError << "\n"
             << "    ||A Z - Z W||_F    = " << frobNormOfError << endl;
    }
}

void TestHermitianEigDouble
( bool testCorrectness, bool printMatrices,
  bool onlyEigenvalues, char range, bool clustered, UpperOrLower uplo, int m, 
  double vl, double vu, int il, int iu, const Grid& g )
{
    DistMatrix<double> A(g), AOrig(g), Z(g);
    DistMatrix<double,VR,STAR> w(g);

    if( clustered )
        Wilkinson( m/2, A );
    else
        HermitianUniformSpectrum( m, A, -10, 10 );
    if( testCorrectness && !onlyEigenvalues )
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
        cout << "  Starting Hermitian eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            HermitianEig( uplo, A, w );
        else if( range == 'I' )
            HermitianEig( uplo, A, w, il, iu );
        else
            HermitianEig( uplo, A, w, vl, vu );
    }
    else
    {
        if( range == 'A' )
            HermitianEig( uplo, A, w, Z );
        else if( range == 'I' )
            HermitianEig( uplo, A, w, Z, il, iu );
        else
            HermitianEig( uplo, A, w, Z, vl, vu );
    }
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << endl;
    }
    if( printMatrices )
    {
        w.Print("eigenvalues:");
        if( !onlyEigenvalues )
            Z.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigenvalues )
        TestCorrectness( printMatrices, uplo, A, w, Z, AOrig );
}
    
void TestHermitianEigDoubleComplex
( bool testCorrectness, bool printMatrices,
  bool onlyEigenvalues, char range, bool clustered, UpperOrLower uplo, int m, 
  double vl, double vu, int il, int iu, const Grid& g )
{
    DistMatrix<Complex<double> > A(g), AOrig(g), Z(g);
    DistMatrix<double,VR,STAR> w(g);

    if( clustered )
        Wilkinson( m/2, A );
    else
        HermitianUniformSpectrum( m, A, -10, 10 );
    if( testCorrectness && !onlyEigenvalues )
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
        cout << "  Starting Hermitian eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            HermitianEig( uplo, A, w );
        else if( range == 'I' )
            HermitianEig( uplo, A, w, il, iu );
        else
            HermitianEig( uplo, A, w, vl, vu );
    }
    else
    {
        if( range == 'A' )
            HermitianEig( uplo, A, w, Z );
        else if( range == 'I' )
            HermitianEig( uplo, A, w, Z, il, iu );
        else
            HermitianEig( uplo, A, w, Z, vl, vu );
    }
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << endl;
    }
    if( printMatrices )
    {
        w.Print("eigenvalues:");
        if( !onlyEigenvalues )
            Z.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigenvalues )
        TestCorrectness( printMatrices, uplo, A, w, Z, AOrig );
}

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
        const bool clustered = atoi(argv[++argNum]);
        const UpperOrLower uplo = CharToUpperOrLower(*argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int nb = atoi(argv[++argNum]);
        const int nbLocalSymv = atoi(argv[++argNum]);
        const bool testCorrectness = atoi(argv[++argNum]);
        const bool printMatrices = atoi(argv[++argNum]);

        if( onlyEigenvalues && testCorrectness && rank==0 )
        {
            cout << "Cannot test correctness with only the eigenvalues." 
                 << endl;
        }

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
        SetLocalSymvBlocksize<double>( nbLocalSymv );
        SetLocalHemvBlocksize<Complex<double> >( nbLocalSymv );

        if( rank == 0 )
        {
            cout << "Will test " << ( uplo==LOWER ? "lower" : "upper" )
                 << " HermitianEig." << endl;
        }

        if( rank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision normal tridiag algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestHermitianEigDouble
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, clustered, uplo, m, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision square tridiag algorithm,\n"
                 << "row-major grid:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianEigDouble
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, clustered, uplo, m, vl, vu, il, iu, g );
 
        if( rank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision square tridiag algorithm,\n"
                 << "col-major grid:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianEigDouble
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, clustered, uplo, m, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision complex normal tridiag algorithm:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestHermitianEigDoubleComplex
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, clustered, uplo, m, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "row-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianEigDoubleComplex
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, clustered, uplo, m, vl, vu, il, iu, g );

        if( rank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "col-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianEigDoubleComplex
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, clustered, uplo, m, vl, vu, il, iu, g );
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

