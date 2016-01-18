/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F>
void TestCorrectness
( bool print,
  UpperOrLower uplo,
  const ElementalMatrix<F>& AOrig,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<Base<F>>& w,
  const ElementalMatrix<F>& Z )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = Z.Height();
    const Int k = Z.Width();

    DistMatrix<F> X(g);
    Identity( X, k, k );
    Herk( uplo, ADJOINT, Real(-1), Z, Real(1), X );
    Real frobNormE = FrobeniusNorm( X );
    if( g.Rank() == 0 )
        Output("    ||Z^H Z - I||_F  = ",frobNormE);
    // X := AZ
    X.AlignWith( Z );
    Zeros( X, n, k );
    Hemm( LEFT, uplo, F(1), AOrig, Z, F(0), X );
    // Find the residual ||X-ZW||_oo = ||AZ-ZW||_oo
    DistMatrix<F> ZW( Z );
    DiagonalScale( RIGHT, NORMAL, w, ZW );
    X -= ZW;
    // Find the infinity norms of A, Z, and AZ-ZW
    Real frobNormA = HermitianFrobeniusNorm( uplo, AOrig );
    frobNormE = FrobeniusNorm( X );
    if( g.Rank() == 0 )
        Output("    ||A Z - Z W||_F / ||A||_F = ",frobNormE/frobNormA);
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=MC>
void TestHermitianEig
( bool testCorrectness,
  bool print,
  bool onlyEigvals,
  bool clustered,
  UpperOrLower uplo,
  Int m, 
  SortType sort,
  const Grid& g,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl,
  bool scalapack )
{
    typedef Base<F> Real;
    DistMatrix<F,U,V> A(g), AOrig(g), Z(g);
    DistMatrix<Real,S,STAR> w(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    if( clustered )
        Wilkinson( A, m/2 );
    else
        HermitianUniformSpectrum( A, m, -10, 10 );
    if( testCorrectness && !onlyEigvals )
        AOrig = A;
    if( print )
        Print( A, "A" );

    if( scalapack && U == MC && V == MR )
    {
        if( onlyEigvals )
        {
            DistMatrix<F,MC,MR,BLOCK> ABlock( A );
            Matrix<Base<F>> wBlock;
            const double startTime = mpi::Time();
            HermitianEig( uplo, ABlock, wBlock, subset );
            const double runTime = mpi::Time() - startTime;
            if( g.Rank() == 0 )
                Output("  ScaLAPACK HermitianEig time: ",runTime);
        }
        else
        {
            DistMatrix<F,MC,MR,BLOCK> ABlock( A ), ZBlock(g);
            Matrix<Base<F>> wBlock;
            const double startTime = mpi::Time();
            HermitianEig( uplo, ABlock, wBlock, ZBlock, subset );
            const double runTime = mpi::Time() - startTime;
            if( g.Rank() == 0 )
                Output("  ScaLAPACK HermitianEig time: ",runTime);
        }
    }

    if( g.Rank() == 0 )
        Output("  Starting Hermitian eigensolver...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( onlyEigvals )
        HermitianEig( uplo, A, w, sort, subset, ctrl );
    else
        HermitianEig( uplo, A, w, Z, sort, subset, ctrl );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
        Output("  Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( Z, "eigenvectors:" );
    }
    if( testCorrectness && !onlyEigvals )
        TestCorrectness( print, uplo, AOrig, A, w, Z );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );
    const int commSize = mpi::Size( comm );

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
        const bool avoidTrmv = 
            Input("--avoidTrmv","avoid Trmv based Symv",true);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
#else
        const bool scalapack = false;
#endif
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
        const bool timeStages = Input("--timeStages","time stages?",true);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        if( range != 'A' && range != 'I' && range != 'V' )
            LogicError("'range' must be 'A', 'I', or 'V'");
        const SortType sort = static_cast<SortType>(sortInt);
        if( onlyEigvals && testCorrectness && commRank==0 )
            Output("Cannot test correctness with only eigenvalues.");
        ComplainIfDebug();

        HermitianEigSubset<double> subset;
        if( range == 'I' )
        {
            subset.indexSubset = true;
            subset.lowerIndex = il;
            subset.upperIndex = iu;
        }
        else if( range == 'V' )
        {
            subset.rangeSubset = true;
            subset.lowerBound = vl;
            subset.upperBound = vu;
        }

        HermitianEigCtrl<double> ctrl_d;
        ctrl_d.timeStages = timeStages;
        ctrl_d.tridiagCtrl.symvCtrl.bsize = nbLocal;
        ctrl_d.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;

        HermitianEigCtrl<Complex<double>> ctrl_z;
        ctrl_z.timeStages = timeStages;
        ctrl_z.tridiagCtrl.symvCtrl.bsize = nbLocal;
        ctrl_z.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;

        if( commRank == 0 )
            Output("Normal tridiag algorithms:");
        ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
        ctrl_z.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
        if( testReal )
            TestHermitianEig<double>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_d, scalapack );
        if( testCpx )
            TestHermitianEig<Complex<double>>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_z, scalapack );

        if( commRank == 0 )
            Output("Square row-major tridiag algorithms:");
        ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_z.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_d.tridiagCtrl.order = ROW_MAJOR;
        ctrl_z.tridiagCtrl.order = ROW_MAJOR;
        if( testReal )
            TestHermitianEig<double>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_d, scalapack );
        if( testCpx )
            TestHermitianEig<Complex<double>>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_z, scalapack );

        if( commRank == 0 )
            Output("Square column-major tridiag algorithms:");
        ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_z.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_d.tridiagCtrl.order = COLUMN_MAJOR;
        ctrl_z.tridiagCtrl.order = COLUMN_MAJOR;
        if( testReal )
            TestHermitianEig<double>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_d, scalapack );
        if( testCpx )
            TestHermitianEig<Complex<double>>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_z, scalapack );

        // Also test with non-standard distributions
        if( commRank == 0 )
            Output("Nonstandard distributions:");
        if( testReal )
            TestHermitianEig<double,MR,MC,MC>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_d, scalapack );
        if( testCpx )
            TestHermitianEig<Complex<double>,MR,MC,MC>
            ( testCorrectness, print, onlyEigvals, clustered, 
              uplo, m, sort, g, subset, ctrl_z, scalapack );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
