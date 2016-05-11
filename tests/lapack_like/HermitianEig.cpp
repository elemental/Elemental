/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
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
    const Real eps = limits::Epsilon<Real>();

    DistMatrix<F> X(g);
    Identity( X, k, k );
    Herk( uplo, ADJOINT, Real(-1), Z, Real(1), X );
    const Real infOrthogError = HermitianInfinityNorm( uplo, X );
    const Real relOrthogError = infOrthogError / (eps*n);
    OutputFromRoot(g.Comm(),"||Z^H Z - I||_oo / (eps n) = ",relOrthogError);

    // X := AZ
    X.AlignWith( Z );
    Zeros( X, n, k );
    Hemm( LEFT, uplo, F(1), AOrig, Z, F(0), X );
    // Find the residual ||X-ZW||_oo = ||AZ-ZW||_oo
    DistMatrix<F> ZW( Z );
    DiagonalScale( RIGHT, NORMAL, w, ZW );
    X -= ZW;
    const Real oneNormA = HermitianOneNorm( uplo, AOrig );
    if( oneNormA == Real(0) )
        LogicError("Tried to test relative accuracy on zero matrix...");
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (n*eps*oneNormA);
    OutputFromRoot(g.Comm(),"||A Z - Z W||_oo / (eps n ||A||_1) = ",relError);

    // TODO: More refined failure conditions
    if( relOrthogError > Real(200) ) // yes, really
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=MC>
void TestHermitianEig
( Int m,
  UpperOrLower uplo,
  bool testCorrectness,
  bool print,
  bool onlyEigvals,
  bool clustered,
  SortType sort,
  const Grid& g,
  const HermitianEigSubset<Base<F>> subset,
  const HermitianEigCtrl<F>& ctrl,
  bool scalapack )
{
    typedef Base<F> Real;
    DistMatrix<F,U,V> A(g), AOrig(g), Z(g);
    DistMatrix<Real,S,STAR> w(g);
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    if( clustered )
        Wilkinson( A, m/2 );
    else
        HermitianUniformSpectrum( A, m, -10, 10 );
    if( testCorrectness && !onlyEigvals )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Timer timer;
    if( scalapack && U == MC && V == MR )
    {
        if( onlyEigvals )
        {
            DistMatrix<F,MC,MR,BLOCK> ABlock( A );
            Matrix<Base<F>> wBlock;
            timer.Start();
            HermitianEig( uplo, ABlock, wBlock, subset );
            const double runTime = timer.Stop();
            OutputFromRoot(g.Comm(),"ScaLAPACK HermitianEig time: ",runTime);
        }
        else
        {
            DistMatrix<F,MC,MR,BLOCK> ABlock( A ), ZBlock(g);
            Matrix<Base<F>> wBlock;
            timer.Start();
            HermitianEig( uplo, ABlock, wBlock, ZBlock, subset );
            const double runTime = timer.Stop();
            OutputFromRoot(g.Comm(),"ScaLAPACK HermitianEig time: ",runTime);
        }
    }

    OutputFromRoot(g.Comm(),"Starting Hermitian eigensolver...");
    mpi::Barrier( g.Comm() );
    timer.Start();
    if( onlyEigvals )
        HermitianEig( uplo, A, w, sort, subset, ctrl );
    else
        HermitianEig( uplo, A, w, Z, sort, subset, ctrl );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( Z, "eigenvectors:" );
    }
    if( testCorrectness && !onlyEigvals )
        TestCorrectness( print, uplo, AOrig, A, w, Z );
    PopIndent();
}


template<typename F>
void TestSuite
( Int m,
  UpperOrLower uplo,
  char range,
  Int il, Int iu,
  Base<F> vl, Base<F> vu,
  bool timeStages,
  Int nbLocal,
  bool avoidTrmv,
  bool testCorrectness,
  bool print,
  bool onlyEigvals,
  bool clustered,
  SortType sort,
  const Grid& g,
  bool scalapack )
{
    OutputFromRoot(g.Comm(),"Will test with ",TypeName<F>());
    PushIndent();

    typedef Base<F> Real;
    HermitianEigSubset<Real> subset;
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

    HermitianEigCtrl<F> ctrl;
    ctrl.timeStages = timeStages;
    ctrl.tridiagCtrl.symvCtrl.bsize = nbLocal;
    ctrl.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;

    OutputFromRoot(g.Comm(),"Normal tridiag algorithms:");
    ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
    TestHermitianEig<F>
    ( m, uplo, testCorrectness, print, onlyEigvals, clustered, 
      sort, g, subset, ctrl, scalapack );

    OutputFromRoot(g.Comm(),"Square row-major tridiag algorithms:");
    ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.tridiagCtrl.order = ROW_MAJOR;
    TestHermitianEig<F>
    ( m, uplo, testCorrectness, print, onlyEigvals, clustered, 
      sort, g, subset, ctrl, scalapack );

    OutputFromRoot(g.Comm(),"Square column-major tridiag algorithms:");
    ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.tridiagCtrl.order = COLUMN_MAJOR;
    TestHermitianEig<F>
    ( m, uplo, testCorrectness, print, onlyEigvals, clustered, 
      sort, g, subset, ctrl, scalapack );

    // Also test with non-standard distributions
    OutputFromRoot(g.Comm(),"Nonstandard distributions:");
    TestHermitianEig<F,MR,MC,MC>
    ( m, uplo, testCorrectness, print, onlyEigvals, clustered, 
      sort, g, subset, ctrl, scalapack );

    PopIndent();
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","height of process grid",0);
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

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        if( range != 'A' && range != 'I' && range != 'V' )
            LogicError("'range' must be 'A', 'I', or 'V'");
        const SortType sort = static_cast<SortType>(sortInt);
        if( onlyEigvals && testCorrectness )
            OutputFromRoot
            (g.Comm(),"Cannot test correctness with only eigenvalues.");
        ComplainIfDebug();

        if( testReal )
        {
            TestSuite<float>
            ( m, uplo, range, il, iu, float(vl), float(vu),
              timeStages, nbLocal, avoidTrmv,
              testCorrectness, print, onlyEigvals, clustered,
              sort, g, scalapack );

            TestSuite<double>
            ( m, uplo, range, il, iu, double(vl), double(vu),
              timeStages, nbLocal, avoidTrmv,
              testCorrectness, print, onlyEigvals, clustered,
              sort, g, scalapack );
         }
         if( testCpx )
         {
            TestSuite<Complex<float>>
            ( m, uplo, range, il, iu, float(vl), float(vu),
              timeStages, nbLocal, avoidTrmv,
              testCorrectness, print, onlyEigvals, clustered,
              sort, g, scalapack );

            TestSuite<Complex<double>>
            ( m, uplo, range, il, iu, double(vl), double(vu),
              timeStages, nbLocal, avoidTrmv,
              testCorrectness, print, onlyEigvals, clustered,
              sort, g, scalapack );
         }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
