/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"
#include "elemental/lapack-like/HermitianNorm.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
#include "elemental/matrices/Wilkinson.hpp"
using namespace std;
using namespace elem;

void TestCorrectness
( bool print,
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
( bool print,
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
( bool testCorrectness, bool print,
  bool onlyEigvals, char range, bool clustered, UpperOrLower uplo, int m, 
  double vl, double vu, int il, int iu, const Grid& g )
{
    DistMatrix<double> A(g), AOrig(g), Z(g);
    DistMatrix<double,VR,STAR> w(g);

    if( clustered )
        Wilkinson( m/2, A );
    else
        HermitianUniformSpectrum( m, A, -10, 10 );
    if( testCorrectness && !onlyEigvals )
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
    if( print )
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting Hermitian eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( onlyEigvals )
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
    if( print )
    {
        w.Print("eigenvalues:");
        if( !onlyEigvals )
            Z.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigvals )
        TestCorrectness( print, uplo, A, w, Z, AOrig );
}
    
void TestHermitianEigDoubleComplex
( bool testCorrectness, bool print,
  bool onlyEigvals, char range, bool clustered, UpperOrLower uplo, int m, 
  double vl, double vu, int il, int iu, const Grid& g )
{
    DistMatrix<Complex<double> > A(g), AOrig(g), Z(g);
    DistMatrix<double,VR,STAR> w(g);

    if( clustered )
        Wilkinson( m/2, A );
    else
        HermitianUniformSpectrum( m, A, -10, 10 );
    if( testCorrectness && !onlyEigvals )
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
    if( print )
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting Hermitian eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( onlyEigvals )
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
    if( print )
    {
        w.Print("eigenvalues:");
        if( !onlyEigvals )
            Z.Print("eigenvectors:");
    }
    if( testCorrectness && !onlyEigvals )
        TestCorrectness( print, uplo, A, w, Z, AOrig );
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
        const bool clustered = Input
            ("--cluster","force clustered eigenvalues?",false);
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
        SetLocalSymvBlocksize<Complex<double> >( nbLocal );
        if( range != 'A' && range != 'I' && range != 'V' )
            throw logic_error("'range' must be 'A', 'I', or 'V'");
        if( onlyEigvals && testCorrectness && commRank==0 )
            cout << "Cannot test correctness with only eigenvalues." << endl;
#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test HermitianEig " << uploChar << endl;

        if( commRank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision normal tridiag algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestHermitianEigDouble
        ( testCorrectness, print, 
          onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision square tridiag algorithm,\n"
                 << "row-major grid:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianEigDouble
        ( testCorrectness, print, 
          onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, g );
 
        if( commRank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision square tridiag algorithm,\n"
                 << "col-major grid:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianEigDouble
        ( testCorrectness, print, 
          onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision complex normal tridiag algorithm:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestHermitianEigDoubleComplex
        ( testCorrectness, print, 
          onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, g );

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
        TestHermitianEigDoubleComplex
        ( testCorrectness, print, 
          onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, g );

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
        TestHermitianEigDoubleComplex
        ( testCorrectness, print, 
          onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, g );
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
