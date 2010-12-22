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
    cout << "Generates random Hermitian matrix then solves for its eigenpairs."
	 << "\n\n"
         << "  HermitianEig <r> <c> <only eigenvalues?> <range> <a> <b> "
            "<highAccuracy?> <shape> <m> <nb> <correctness?> <print?>\n\n"
         << "  r: number of process rows\n"
         << "  c: number of process cols\n"
         << "  only eigenvalues?: 0/1\n"
         << "  range: 'A' for all, 'I' for index range, "
            "'V' for floating-point range\n"
         << "  a: if range=='I', 0-indexed first eigenpair to compute\n"
            "     if range=='V', lower-bound on eigenvalues\n"
         << "  b: if range=='I', 0-indexed last eigenpair to compute\n"
            "     if range=='V', upper-bound on eigenvalues\n"
         << "  highAccuracy?: try for high acc. iff != 0\n"
         << "  shape: L/U\n"
         << "  m: height of matrix\n"
         << "  nb: algorithmic blocksize\n"
         << "  test correctness?: false iff 0\n"
         << "  print matrices?: false iff 0\n" << endl;
}

void TestCorrectnessDouble
( bool printMatrices,
  Shape shape,
  const DistMatrix<double,MC,  MR>& A,
  const DistMatrix<double,Star,VR>& w,
  const DistMatrix<double,MC,  MR>& Z,
  const DistMatrix<double,MC  ,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int n = Z.Height();
    const int k = Z.Width();

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

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of eigenvectors..." << endl;
    DistMatrix<double,MC,MR> X(k,k,g);
    X.SetToIdentity();
    blas::Herk( shape, ConjugateTranspose, (double)-1, Z, (double)1, X );
    double oneNormOfError = lapack::OneNorm( X );
    double infNormOfError = lapack::InfinityNorm( X );
    double frobNormOfError = lapack::FrobeniusNorm( X );
    if( g.VCRank() == 0 )
    {
        cout << "    ||Z^H Z - I||_1  = " << oneNormOfError << "\n"
             << "    ||Z^H Z - I||_oo = " << infNormOfError << "\n"
             << "    ||Z^H Z - I||_F  = " << frobNormOfError << endl;
    }

    if( g.VCRank() == 0 )
        cout << "  Testing for deviation of AZ from ZW..." << endl;
    // Set X := AZ
    X.AlignWith( Z );
    X.ResizeTo( n, k );
    blas::Hemm( Left, shape, (double)1, AOrig, Z, (double)0, X );
    // Set X := X - ZW = AZ - ZW
    for( int j=0; j<X.LocalWidth(); ++j )
    {
        double omega = w_Star_MR.GetLocalEntry(0,j);
        for( int i=0; i<X.LocalHeight(); ++i )
            X.SetLocalEntry(i,j,
                X.GetLocalEntry(i,j)-omega*Z.GetLocalEntry(i,j));
    }
    // Find the infinity norms of A, Z, and AZ-ZW
    double infNormOfA = lapack::HermitianInfinityNorm( shape, AOrig );
    double frobNormOfA = lapack::HermitianFrobeniusNorm( shape, AOrig );
    double oneNormOfZ = lapack::OneNorm( Z );
    double infNormOfZ = lapack::InfinityNorm( Z );
    double frobNormOfZ = lapack::FrobeniusNorm( Z );
    oneNormOfError = lapack::OneNorm( X );
    infNormOfError = lapack::InfinityNorm( X );
    frobNormOfError = lapack::FrobeniusNorm( X );
    if( g.VCRank() == 0 )
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

#ifndef WITHOUT_COMPLEX
void TestCorrectnessDoubleComplex
( bool printMatrices,
  Shape shape,
  const DistMatrix<std::complex<double>,MC,  MR>& A,
  const DistMatrix<             double, Star,VR>& w,
  const DistMatrix<std::complex<double>,MC,  MR>& Z,
  const DistMatrix<std::complex<double>,MC  ,MR>& AOrig )
{
    const Grid& g = A.Grid();
    const int n = Z.Height();
    const int k = Z.Width();

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<double,Star,MR> w_Star_MR(true,Z.RowAlignment(),g); 
    w_Star_MR = w;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of eigenvectors..." << endl;
    DistMatrix<std::complex<double>,MC,MR> X( k, k, g );
    X.SetToIdentity();
    blas::Herk
    ( shape, ConjugateTranspose, 
      std::complex<double>(-1), Z, std::complex<double>(1), X );
    double oneNormOfError = lapack::OneNorm( X );
    double infNormOfError = lapack::InfinityNorm( X );
    double frobNormOfError = lapack::FrobeniusNorm( X );
    if( g.VCRank() == 0 )
    {
        cout << "    ||Z^H Z - I||_1  = " << oneNormOfError << "\n"
             << "    ||Z^H Z - I||_oo = " << infNormOfError << "\n"
             << "    ||Z^H Z - I||_F  = " << frobNormOfError << endl;
    }

    if( g.VCRank() == 0 )
        cout << "  Testing for deviation of AZ from ZW..." << endl;
    // X := AZ
    X.AlignWith( Z );
    X.ResizeTo( n, k );
    blas::Hemm
    ( Left, shape, std::complex<double>(1), AOrig, Z, 
      std::complex<double>(0), X );
    // Find the residual ||X-ZW||_oo = ||AZ-ZW||_oo
    for( int j=0; j<X.LocalWidth(); ++j )
    {
        double omega = w_Star_MR.GetLocalEntry(0,j);
        for( int i=0; i<X.LocalHeight(); ++i )
            X.SetLocalEntry(i,j,
                X.GetLocalEntry(i,j)-omega*Z.GetLocalEntry(i,j));
    }
    // Find the infinity norms of A, Z, and AZ-ZW
    double infNormOfA = lapack::HermitianInfinityNorm( shape, AOrig );
    double frobNormOfA = lapack::HermitianFrobeniusNorm( shape, AOrig );
    double oneNormOfZ = lapack::OneNorm( Z );
    double infNormOfZ = lapack::InfinityNorm( Z );
    double frobNormOfZ = lapack::FrobeniusNorm( Z );
    oneNormOfError = lapack::OneNorm( X );
    infNormOfError = lapack::InfinityNorm( X );
    frobNormOfError = lapack::FrobeniusNorm( X );
    if( g.VCRank() == 0 )
    {
        cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
             << "    ||A||_F            = " << frobNormOfA << "\n"
             << "    ||Z||_1            = " << oneNormOfZ << "\n"
             << "    ||Z||_oo           = " << infNormOfZ << "\n"
             << "    ||Z||_F            = " << frobNormOfZ << "\n"
             << "    ||A Z - Z W||_oo   = " << infNormOfError << "\n"
             << "    ||A Z - Z W||_F    = " << frobNormOfError << endl;
    }
}
#endif // WITHOUT_COMPLEX

void TestHermitianEigDouble
( bool testCorrectness, bool printMatrices,
  bool onlyEigenvalues, char range, Shape shape, int m, 
  double vl, double vu, int il, int iu, bool tryForHighAccuracy, 
  const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<double,MC,MR> A(m,m,g);
    DistMatrix<double,MC,MR> AOrig(g);
    DistMatrix<double,Star,VR> w(g);
    DistMatrix<double,MC,MR> Z(g);

    A.SetToRandomHPD();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting Hermitian eigensolver...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            lapack::HermitianEig( shape, A, w, tryForHighAccuracy );
        else if( range == 'I' )
            lapack::HermitianEig( shape, A, w, il, iu, tryForHighAccuracy );
        else
            lapack::HermitianEig( shape, A, w, vl, vu, tryForHighAccuracy );
    }
    else
    {
        if( range == 'A' )
            lapack::HermitianEig( shape, A, w, Z, tryForHighAccuracy );
        else if( range == 'I' )
            lapack::HermitianEig( shape, A, w, Z, il, iu, tryForHighAccuracy );
        else
            lapack::HermitianEig( shape, A, w, Z, vl, vu, tryForHighAccuracy );
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
            Z.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigenvalues )
    {
        TestCorrectnessDouble( printMatrices, shape, A, w, Z, AOrig );
    }
}
    
#ifndef WITHOUT_COMPLEX
void TestHermitianEigDoubleComplex
( bool testCorrectness, bool printMatrices,
  bool onlyEigenvalues, char range, Shape shape, int m, 
  double vl, double vu, int il, int iu, bool tryForHighAccuracy, 
  const Grid& g )
{
    double startTime, endTime, runTime;
    DistMatrix<std::complex<double>,MC,  MR> A(m,m,g);
    DistMatrix<std::complex<double>,MC,  MR> AOrig(g);
    DistMatrix<             double, Star,VR> w(g);
    DistMatrix<std::complex<double>,MC,  MR> Z(g);

    A.SetToRandomHPD();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting Hermitian eigensolver...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    if( onlyEigenvalues )
    {
        if( range == 'A' )
            lapack::HermitianEig( shape, A, w, tryForHighAccuracy );
        else if( range == 'I' )
            lapack::HermitianEig( shape, A, w, il, iu, tryForHighAccuracy );
        else
            lapack::HermitianEig( shape, A, w, vl, vu, tryForHighAccuracy );
    }
    else
    {
        if( range == 'A' )
            lapack::HermitianEig( shape, A, w, Z, tryForHighAccuracy );
        else if( range == 'I' )
            lapack::HermitianEig( shape, A, w, Z, il, iu, tryForHighAccuracy );
        else
            lapack::HermitianEig( shape, A, w, Z, vl, vu, tryForHighAccuracy );
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
            Z.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigenvalues )
    {
        TestCorrectnessDoubleComplex( printMatrices, shape, A, w, Z, AOrig );
    }
}
#endif // WITHOUT_COMPLEX

int main( int argc, char* argv[] )
{
    int rank;
    Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 13 )
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
        TestHermitianEigDouble
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, shape, m, vl, vu, il, iu, 
          tryForHighAccuracy, g );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestHermitianEigDoubleComplex
        ( testCorrectness, printMatrices, 
          onlyEigenvalues, range, shape, m, vl, vu, il, iu, 
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

