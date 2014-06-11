/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_HERMITIANUNIFORMSPECTRUM_INC
#include EL_IDENTITY_INC
#include EL_WILKINSON_INC
#include EL_ZEROS_INC
using namespace std;
using namespace El;

template<typename F>
void TestCorrectness
( bool print,
  UpperOrLower uplo,
  const DistMatrix<F>& A,
  const DistMatrix<Base<F>,VR,STAR>& w,
  const DistMatrix<F>& Z,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = Z.Height();
    const Int k = Z.Width();

    if( g.Rank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<Real,MR,STAR> w_MR_STAR(true,Z.RowAlign(),g); 
    w_MR_STAR = w;
    if( g.Rank() == 0 )
        cout << "DONE" << endl;

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of eigenvectors..." << endl;
    DistMatrix<F> X(g);
    Identity( X, k, k );
    Herk( uplo, ADJOINT, F(-1), Z, F(1), X );
    Real oneNormOfError = OneNorm( X );
    Real infNormOfError = InfinityNorm( X );
    Real frobNormOfError = FrobeniusNorm( X );
    if( g.Rank() == 0 )
    {
        cout << "    ||Z^H Z - I||_1  = " << oneNormOfError << "\n"
             << "    ||Z^H Z - I||_oo = " << infNormOfError << "\n"
             << "    ||Z^H Z - I||_F  = " << frobNormOfError << "\n\n"
             << "  Testing for deviation of AZ from ZW..." << endl;
    }
    // X := AZ
    X.AlignWith( Z );
    Zeros( X, n, k );
    Hemm( LEFT, uplo, F(1), AOrig, Z, F(0), X );
    // Find the residual ||X-ZW||_oo = ||AZ-ZW||_oo
    for( Int jLoc=0; jLoc<X.LocalWidth(); ++jLoc )
    {
        const Real omega = w_MR_STAR.GetLocal(jLoc,0);
        for( Int iLoc=0; iLoc<X.LocalHeight(); ++iLoc )
        {
            const F chi = X.GetLocal(iLoc,jLoc);
            const F zeta = Z.GetLocal(iLoc,jLoc);
            X.SetLocal(iLoc,jLoc,chi-omega*zeta);
        }
    }
    // Find the infinity norms of A, Z, and AZ-ZW
    Real infNormOfA = HermitianInfinityNorm( uplo, AOrig );
    Real frobNormOfA = HermitianFrobeniusNorm( uplo, AOrig );
    Real oneNormOfZ = OneNorm( Z );
    Real infNormOfZ = InfinityNorm( Z );
    Real frobNormOfZ = FrobeniusNorm( Z );
    oneNormOfError = OneNorm( X );
    infNormOfError = InfinityNorm( X );
    frobNormOfError = FrobeniusNorm( X );
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

template<typename F>
void TestHermitianEig
( bool testCorrectness, bool print,
  bool onlyEigvals, char range, bool clustered, UpperOrLower uplo, Int m, 
  Base<F> vl, Base<F> vu, Int il, Int iu, SortType sort, const Grid& g,
  const HermitianEigCtrl<Base<F>> ctrl )
{
    typedef Base<F> Real;
    DistMatrix<F> A(g), AOrig(g), Z(g);
    DistMatrix<Real,VR,STAR> w(g);

    if( clustered )
        Wilkinson( A, m/2 );
    else
        HermitianUniformSpectrum( A, m, -10, 10 );
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
        Print( A, "A" );

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
            HermitianEig( uplo, A, w, sort, ctrl );
        else if( range == 'I' )
            HermitianEig( uplo, A, w, il, iu, sort, ctrl );
        else
            HermitianEig( uplo, A, w, vl, vu, sort, ctrl );
    }
    else
    {
        if( range == 'A' )
            HermitianEig( uplo, A, w, Z, sort, ctrl );
        else if( range == 'I' )
            HermitianEig( uplo, A, w, Z, il, iu, sort, ctrl );
        else
            HermitianEig( uplo, A, w, Z, vl, vu, sort, ctrl );
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
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( Z, "eigenvectors:" );
    }
    if( testCorrectness && !onlyEigvals )
        TestCorrectness( print, uplo, A, w, Z, AOrig );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const bool onlyEigvals = Input
            ("--onlyEigvals","only compute eigenvalues?",false);
        const char range = Input
            ("--range",
             "range of eigenpairs: 'A' for all, 'I' for index range, "
             "'V' for value range",'A');
        const Int il = Input("--il","lower bound of index range",0);
        const Int iu = Input("--iu","upper bound of index range",100);
        const double vl = Input("--vl","lower bound of value range",0.);
        const double vu = Input("--vu","upper bound of value range",100.);
        const Int sortInt = Input("--sort","sort type",0);
        const bool clustered = Input
            ("--cluster","force clustered eigenvalues?",false);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalSymvBlocksize<double>( nbLocal );
        SetLocalSymvBlocksize<Complex<double>>( nbLocal );
        if( range != 'A' && range != 'I' && range != 'V' )
            LogicError("'range' must be 'A', 'I', or 'V'");
        const SortType sort = static_cast<SortType>(sortInt);
        if( onlyEigvals && testCorrectness && commRank==0 )
            cout << "Cannot test correctness with only eigenvalues." << endl;
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test HermitianEig " << uploChar << endl;

        HermitianEigCtrl<double> ctrl;

        if( commRank == 0 )
            cout << "Normal tridiag algorithms:" << endl;
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
        if( testReal )
            TestHermitianEig<double>
            ( testCorrectness, print, 
              onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, sort, g,
              ctrl );
        if( testCpx )
            TestHermitianEig<Complex<double>>
            ( testCorrectness, print, 
              onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, sort, g,
              ctrl );

        if( commRank == 0 )
            cout << "Square row-major tridiag algorithms:" << endl;
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl.tridiagCtrl.order = ROW_MAJOR;
        if( testReal )
            TestHermitianEig<double>
            ( testCorrectness, print, 
              onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, sort, g,
              ctrl );
        if( testCpx )
            TestHermitianEig<Complex<double>>
            ( testCorrectness, print, 
              onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, sort, g,
              ctrl );

        if( commRank == 0 )
            cout << "Square column-major tridiag algorithms:" << endl;
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl.tridiagCtrl.order = COLUMN_MAJOR;
        if( testReal )
            TestHermitianEig<double>
            ( testCorrectness, print, 
              onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, sort, g,
              ctrl );
        if( testCpx )
            TestHermitianEig<Complex<double>>
            ( testCorrectness, print, 
              onlyEigvals, range, clustered, uplo, m, vl, vu, il, iu, sort, g,
              ctrl );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
