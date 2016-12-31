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
  const Matrix<F>& AOrig,
  const Matrix<F>& A,
  const Matrix<Base<F>>& w,
  const Matrix<F>& Q )
{
    typedef Base<F> Real;
    const Int n = Q.Height();
    const Int k = Q.Width();
    const Real eps = limits::Epsilon<Real>();

    Matrix<F> X;
    Identity( X, k, k );
    Herk( uplo, ADJOINT, Real(-1), Q, Real(1), X );
    const Real infOrthogError = HermitianInfinityNorm( uplo, X );
    const Real relOrthogError = infOrthogError / (eps*n);
    Output("||Q^H Q - I||_oo / (eps n) = ",relOrthogError);

    // X := A Q
    Zeros( X, n, k );
    Hemm( LEFT, uplo, F(1), AOrig, Q, F(0), X );
    // Find the residual ||X-QW||_oo = ||AQ-QW||_oo
    Matrix<F> QW( Q );
    DiagonalScale( RIGHT, NORMAL, w, QW );
    X -= QW;
    const Real oneNormA = HermitianOneNorm( uplo, AOrig );
    if( oneNormA == Real(0) )
        LogicError("Tried to test relative accuracy on zero matrix...");
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (n*eps*oneNormA);
    Output("||A Q - Q W||_oo / (eps n ||A||_1) = ",relError);

    // TODO: More refined failure conditions
    if( relOrthogError > Real(200) ) // yes, really
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename F>
void TestCorrectness
( bool print,
  UpperOrLower uplo,
  const AbstractDistMatrix<F>& AOrig,
  const AbstractDistMatrix<F>& A,
  const AbstractDistMatrix<Base<F>>& w,
  const AbstractDistMatrix<F>& Q )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = Q.Height();
    const Int k = Q.Width();
    const Real eps = limits::Epsilon<Real>();

    DistMatrix<F> X(g);
    Identity( X, k, k );
    Herk( uplo, ADJOINT, Real(-1), Q, Real(1), X );
    const Real infOrthogError = HermitianInfinityNorm( uplo, X );
    const Real relOrthogError = infOrthogError / (eps*n);
    OutputFromRoot(g.Comm(),"||Q^H Q - I||_oo / (eps n) = ",relOrthogError);

    // X := A Q
    X.AlignWith( Q );
    Zeros( X, n, k );
    Hemm( LEFT, uplo, F(1), AOrig, Q, F(0), X );
    // Find the residual ||X-QW||_oo = ||AQ-QW||_oo
    DistMatrix<F> QW( Q );
    DiagonalScale( RIGHT, NORMAL, w, QW );
    X -= QW;
    const Real oneNormA = HermitianOneNorm( uplo, AOrig );
    if( oneNormA == Real(0) )
        LogicError("Tried to test relative accuracy on zero matrix...");
    const Real infError = InfinityNorm( X );
    const Real relError = infError / (n*eps*oneNormA);
    OutputFromRoot(g.Comm(),"||A Q - Q W||_oo / (eps n ||A||_1) = ",relError);

    // TODO: More refined failure conditions
    if( relOrthogError > Real(200) ) // yes, really
        LogicError("Relative orthogonality error was unacceptably large");
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename F>
void TestHermitianEigSequential
( Int m,
  UpperOrLower uplo,
  bool onlyEigvals,
  bool clustered,
  bool correctness,
  bool print,
  const HermitianEigCtrl<F>& ctrl )
{
    typedef Base<F> Real;
    Matrix<F> A, AOrig, Q;
    Matrix<Real> w;
    Output("Testing with ",TypeName<F>());
    PushIndent();

    if( clustered )
        Wilkinson( A, m/2 );
    else
        HermitianUniformSpectrum( A, m, -10, 10 );
    if( correctness && !onlyEigvals )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Timer timer;
    Output("Starting Hermitian eigensolver...");
    timer.Start();
    if( onlyEigvals )
        HermitianEig( uplo, A, w, ctrl );
    else
        HermitianEig( uplo, A, w, Q, ctrl );
    const double runTime = timer.Stop();
    Output("Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( Q, "eigenvectors:" );
    }
    if( correctness && !onlyEigvals )
        TestCorrectness( print, uplo, AOrig, A, w, Q );
    PopIndent();
}
template<typename F,Dist U=MC,Dist V=MR,Dist S=MC>
void TestHermitianEig
( Int m,
  UpperOrLower uplo,
  bool onlyEigvals,
  bool clustered,
  bool correctness,
  bool print,
  const Grid& g,
  const HermitianEigCtrl<F>& ctrl )
{
    typedef Base<F> Real;
    DistMatrix<F,U,V> A(g), AOrig(g), Q(g);
    DistMatrix<Real,S,STAR> w(g);
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    if( clustered )
        Wilkinson( A, m/2 );
    else
        HermitianUniformSpectrum( A, m, -10, 10 );
    if( correctness && !onlyEigvals )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Timer timer;
    OutputFromRoot(g.Comm(),"Starting Hermitian eigensolver...");
    mpi::Barrier( g.Comm() );
    timer.Start();
    if( onlyEigvals )
        HermitianEig( uplo, A, w, ctrl );
    else
        HermitianEig( uplo, A, w, Q, ctrl );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( Q, "eigenvectors:" );
    }
    if( correctness && !onlyEigvals )
        TestCorrectness( print, uplo, AOrig, A, w, Q );
    PopIndent();
}

template<typename F>
void TestSuite
( Int m,
  UpperOrLower uplo,
  bool onlyEigvals,
  bool clustered,
  bool sequential,
  bool distributed,
  bool correctness,
  bool print,
  const Grid& g,
  const HermitianEigCtrl<double>& ctrlDbl )
{
    typedef Base<F> Real;
    OutputFromRoot(g.Comm(),"Will test with ",TypeName<F>());
    PushIndent();

    auto subsetDbl = ctrlDbl.tridiagEigCtrl.subset;
    HermitianEigSubset<Real> subset;
    subset.indexSubset = subsetDbl.indexSubset;
    subset.lowerIndex = subsetDbl.lowerIndex;
    subset.upperIndex = subsetDbl.upperIndex;
    subset.rangeSubset = subsetDbl.rangeSubset;
    subset.lowerBound = subsetDbl.lowerBound;
    subset.upperBound = subsetDbl.upperBound;

    HermitianEigCtrl<F> ctrl;
    ctrl.timeStages = ctrlDbl.timeStages;
    ctrl.useScaLAPACK = ctrlDbl.useScaLAPACK;
    ctrl.tridiagCtrl.symvCtrl.bsize =
      ctrlDbl.tridiagCtrl.symvCtrl.bsize;
    ctrl.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv =
      ctrlDbl.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv;
    ctrl.tridiagEigCtrl.sort = ctrlDbl.tridiagEigCtrl.sort;
    ctrl.tridiagEigCtrl.alg = ctrlDbl.tridiagEigCtrl.alg;
    ctrl.tridiagEigCtrl.subset = subset;
    ctrl.tridiagEigCtrl.progress = ctrlDbl.tridiagEigCtrl.progress;

    if( sequential && g.Rank() == 0 )
    {
        TestHermitianEigSequential<F>
        ( m, uplo, onlyEigvals, clustered, correctness, print, ctrl );
    }
    if( distributed )
    {
        OutputFromRoot(g.Comm(),"Normal tridiag algorithms:");
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
        TestHermitianEig<F>
        ( m, uplo, onlyEigvals, clustered, correctness, print, g, ctrl );

        OutputFromRoot(g.Comm(),"Square row-major tridiag algorithms:");
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl.tridiagCtrl.order = ROW_MAJOR;
        TestHermitianEig<F>
        ( m, uplo, onlyEigvals, clustered, correctness, print, g, ctrl );

        OutputFromRoot(g.Comm(),"Square column-major tridiag algorithms:");
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl.tridiagCtrl.order = COLUMN_MAJOR;
        TestHermitianEig<F>
        ( m, uplo, onlyEigvals, clustered, correctness, print, g, ctrl );

        // Also test with non-standard distributions
        OutputFromRoot(g.Comm(),"Nonstandard distributions:");
        TestHermitianEig<F,MR,MC,MC>
        ( m, uplo, onlyEigvals, clustered, correctness, print, g, ctrl );
    }

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
        const bool onlyEigvals =
          Input("--onlyEigvals","only compute eigenvalues?",false);
        const char range =
          Input("--range",
           "range of eigenpairs: 'A' for all, 'I' for index range, "
           "'V' for value range",'A');
        const Int il = Input("--il","lower bound of index range",0);
        const Int iu = Input("--iu","upper bound of index range",100);
        const double vl = Input("--vl","lower bound of value range",0.);
        const double vu = Input("--vu","upper bound of value range",100.);
        const Int sortInt = Input("--sort","sort type",0);
        const bool clustered =
          Input("--cluster","force clustered eigenvalues?",false);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool avoidTrmv =
          Input("--avoidTrmv","avoid Trmv based Symv",true);
        const bool useScaLAPACK =
          Input("--useScaLAPACK","test ScaLAPACK?",false);
        const Int algInt = Input("--algInt","0: QR, 1: D&C, 2: MRRR",1);
        const bool sequential =
          Input("--sequential","test sequential?",true);
        const bool distributed =
          Input("--distributed","test distributed?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool progress = Input("--progress","print progress?",false);
        const bool print = Input("--print","print matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
        const bool timeStages = Input("--timeStages","time stages?",true);
        ProcessInput();
        PrintInputReport();

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = colMajor ? COLUMN_MAJOR : ROW_MAJOR;
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const auto alg = static_cast<HermitianTridiagEigAlg>(algInt);
        SetBlocksize( nb );
        if( range != 'A' && range != 'I' && range != 'V' )
            LogicError("'range' must be 'A', 'I', or 'V'");
        const SortType sort = static_cast<SortType>(sortInt);
        if( onlyEigvals && correctness )
            OutputFromRoot
            (g.Comm(),"Cannot test correctness with only eigenvalues.");
        ComplainIfDebug();

        // Convert an initial double-precision control structure into each of
        // the datatypes for simplicity
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

        HermitianEigCtrl<double> ctrl;
        ctrl.timeStages = timeStages;
        ctrl.useScaLAPACK = useScaLAPACK;
        ctrl.tridiagCtrl.symvCtrl.bsize = nbLocal;
        ctrl.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;
        ctrl.tridiagEigCtrl.sort = sort;
        ctrl.tridiagEigCtrl.alg = alg;
        ctrl.tridiagEigCtrl.subset = subset;
        ctrl.tridiagEigCtrl.progress = progress;

        if( testReal )
        {
            TestSuite<float>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<double>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );

#ifdef EL_HAVE_QD
            TestSuite<DoubleDouble>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<QuadDouble>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );
#endif

#ifdef EL_HAVE_QUAD
            TestSuite<Quad>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );
#endif
#ifdef EL_HAVE_MPC
            TestSuite<BigFloat>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );
#endif
         }
         if( testCpx )
         {
            TestSuite<Complex<float>>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<Complex<double>>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );

#ifdef EL_HAVE_QD
            TestSuite<Complex<DoubleDouble>>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<Complex<QuadDouble>>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );
#endif

#ifdef EL_HAVE_QUAD
            TestSuite<Complex<Quad>>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );
#endif

#ifdef EL_HAVE_MPC
            TestSuite<Complex<BigFloat>>
            ( m, uplo, onlyEigvals, clustered,
              sequential, distributed, correctness, print, g, ctrl );
#endif
         }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
