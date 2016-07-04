/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace std;
using namespace El;

template<typename F>
void TestCorrectness
( bool print,
  Pencil pencil, UpperOrLower uplo,
  const Matrix<F>& AOrig,
  const Matrix<F>& BOrig,
  const Matrix<F>& A,
  const Matrix<F>& B,
  const Matrix<Base<F>>& w,
  const Matrix<F>& X )
{
    typedef Base<F> Real;
    const Int n = X.Height();
    const Int k = X.Width();
    const Real eps = limits::Epsilon<Real>();

    const Real oneNormA = HermitianOneNorm( uplo, AOrig );
    const Real oneNormB = HermitianOneNorm( uplo, BOrig );

    if( pencil == AXBX )
    {
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        Matrix<F> Y;
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        DiagonalScale( RIGHT, NORMAL, w, Y );
        // Y := Y - AX = BXW - AX
        Hemm( LEFT, uplo, F(-1), AOrig, X, F(1), Y );
        const Real infError = InfinityNorm( Y );
        const Real relError = infError / (eps*n*Max(oneNormA,oneNormB));
        Output("||A X - B X W||_oo / (eps n max(||A||_1,||B||_1) = ", relError);

        Matrix<F> Z;
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        const Real infOrthogError = HermitianInfinityNorm( uplo, Y );
        const Real relOrthogError = infOrthogError / (eps*n*oneNormB);
        Output("||X^H B X - I||_oo  / (eps n ||B||_1) = ", relOrthogError);

        // TODO: More refined failure conditions
        if( relOrthogError > Real(200) )
            LogicError("Relative orthogonality error was unacceptable");
        if( relError > Real(100) )
            LogicError("Relative error was unacceptably large");
    }
    else if( pencil == ABX )
    {
        // Set Y := BX
        Matrix<F> Y;
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        // Set Z := AY = ABX
        Matrix<F> Z( n, k );
        Hemm( LEFT, uplo, F(1), AOrig, Y, F(0), Z );
        // Set Z := Z - XW = ABX - XW
        for( Int j=0; j<Z.Width(); ++j )
        {
            const Real omega = w(j); 
            for( Int i=0; i<Z.Height(); ++i )
            {
                Z(i,j) -= omega*X(i,j);
            }
        }
        const Real infError = InfinityNorm( Z );
        const Real relError = infError / (eps*n*Max(oneNormA,oneNormB));
        Output
        ("||A B X - X W||_oo  / (eps n Max(||A||_1,||B||_1)) = ",relError);

        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        const Real infOrthogError = HermitianInfinityNorm( uplo, Y );
        const Real relOrthogError = infOrthogError / (eps*n*oneNormB);
        Output("||X^H B X - I||_oo / (eps n ||B||_1) = ",relOrthogError);

        // TODO: More refined failure conditions
        if( relOrthogError > Real(200) )
            LogicError("Relative orthogonality error was unacceptable");
        if( relError > Real(100) )
            LogicError("Relative error was unacceptably large");
    }
    else /* pencil == BAX */
    {
        // Set Y := AX
        Matrix<F> Y;
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
        // Set Z := BY = BAX
        Matrix<F> Z( n, k );
        Hemm( LEFT, uplo, F(1), BOrig, Y, F(0), Z );
        // Set Z := Z - XW = BAX-XW
        for( Int j=0; j<Z.Width(); ++j )
        {
            const Real omega = w(j); 
            for( Int i=0; i<Z.Height(); ++i )
            {
                Z(i,j) -= omega*X(i,j);
            }
        }
        const Real infError = InfinityNorm( Z );
        const Real relError = infError / (eps*n*Max(oneNormA,oneNormB));
        Output
        ("||B A X - X W||_oo / (eps n Max(||A||_1,||B||_1)) = ",relError);

        Z = X;
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), B, Z );
        else
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        const Real infOrthogError = HermitianInfinityNorm( uplo, Y );
        const Real relOrthogError = infOrthogError / (eps*n*oneNormB);
        Output("||X^H B^-1 X - I||_oo  / (eps n ||B||_1) = ", relOrthogError);

        // TODO: More refined failure conditions
        if( relOrthogError > Real(200) )
            LogicError("Relative orthogonality error was unacceptable");
        if( relError > Real(100) )
            LogicError("Relative error was unacceptably large");
    }
}

template<typename F>
void TestCorrectness
( bool print,
  Pencil pencil, UpperOrLower uplo,
  const ElementalMatrix<F>& AOrig,
  const ElementalMatrix<F>& BOrig,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& B,
  const ElementalMatrix<Base<F>>& w,
  const ElementalMatrix<F>& X )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = X.Height();
    const Int k = X.Width();
    const Real eps = limits::Epsilon<Real>();

    const Real oneNormA = HermitianOneNorm( uplo, AOrig );
    const Real oneNormB = HermitianOneNorm( uplo, BOrig );

    DistMatrix<Real,MR,STAR> w_MR_STAR(true,X.RowAlign(),g); 
    w_MR_STAR = w;

    if( pencil == AXBX )
    {
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        DiagonalScale( RIGHT, NORMAL, w_MR_STAR, Y );
        // Y := Y - AX = BXW - AX
        Hemm( LEFT, uplo, F(-1), AOrig, X, F(1), Y );
        const Real infError = InfinityNorm( Y );
        const Real relError = infError / (eps*n*Max(oneNormA,oneNormB));
        OutputFromRoot
        (g.Comm(),
         "||A X - B X W||_oo / (eps n max(||A||_1,||B||_1) = ", relError);

        DistMatrix<F> Z(g);
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        const Real infOrthogError = HermitianInfinityNorm( uplo, Y );
        const Real relOrthogError = infOrthogError / (eps*n*oneNormB);
        OutputFromRoot
        (g.Comm(),
         "||X^H B X - I||_oo  / (eps n ||B||_1) = ", relOrthogError);

        // TODO: More refined failure conditions
        if( relOrthogError > Real(200) )
            LogicError("Relative orthogonality error was unacceptable");
        if( relError > Real(100) )
            LogicError("Relative error was unacceptably large");
    }
    else if( pencil == ABX )
    {
        // Set Y := BX
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        // Set Z := AY = ABX
        DistMatrix<F> Z( n, k, g );
        Hemm( LEFT, uplo, F(1), AOrig, Y, F(0), Z );
        // Set Z := Z - XW = ABX - XW
        const Matrix<Real>& w_MR_STARLoc = w_MR_STAR.LockedMatrix();
        const Matrix<F>& XLoc = X.LockedMatrix();
              Matrix<F>& ZLoc = Z.Matrix();
        for( Int jLoc=0; jLoc<ZLoc.Width(); ++jLoc )
        {
            const Real omega = w_MR_STARLoc(jLoc); 
            for( Int iLoc=0; iLoc<ZLoc.Height(); ++iLoc )
            {
                ZLoc(iLoc,jLoc) -= omega*XLoc(iLoc,jLoc);
            }
        }
        const Real infError = InfinityNorm( Z );
        const Real relError = infError / (eps*n*Max(oneNormA,oneNormB));
        OutputFromRoot
        (g.Comm(),
         "||A B X - X W||_oo  / (eps n Max(||A||_1,||B||_1)) = ",relError);

        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        const Real infOrthogError = HermitianInfinityNorm( uplo, Y );
        const Real relOrthogError = infOrthogError / (eps*n*oneNormB);
        OutputFromRoot
        (g.Comm(),
         "||X^H B X - I||_oo / (eps n ||B||_1) = ",relOrthogError);

        // TODO: More refined failure conditions
        if( relOrthogError > Real(200) )
            LogicError("Relative orthogonality error was unacceptable");
        if( relError > Real(100) )
            LogicError("Relative error was unacceptably large");
    }
    else /* pencil == BAX */
    {
        // Set Y := AX
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
        // Set Z := BY = BAX
        DistMatrix<F> Z( n, k, g );
        Hemm( LEFT, uplo, F(1), BOrig, Y, F(0), Z );
        // Set Z := Z - XW = BAX-XW
        const Matrix<Real>& w_MR_STARLoc = w_MR_STAR.LockedMatrix();
        const Matrix<F>& XLoc = X.LockedMatrix();
              Matrix<F>& ZLoc = Z.Matrix();
        for( Int jLoc=0; jLoc<ZLoc.Width(); ++jLoc )
        {
            const Real omega = w_MR_STARLoc(jLoc); 
            for( Int iLoc=0; iLoc<ZLoc.Height(); ++iLoc )
            {
                ZLoc(iLoc,jLoc) -= omega*XLoc(iLoc,jLoc);
            }
        }
        const Real infError = InfinityNorm( Z );
        const Real relError = infError / (eps*n*Max(oneNormA,oneNormB));
        OutputFromRoot
        (g.Comm(),
         "||B A X - X W||_oo / (eps n Max(||A||_1,||B||_1)) = ",relError);

        Z = X;
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), B, Z );
        else
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        const Real infOrthogError = HermitianInfinityNorm( uplo, Y );
        const Real relOrthogError = infOrthogError / (eps*n*oneNormB);
        OutputFromRoot
        (g.Comm(),
         "||X^H B^-1 X - I||_oo  / (eps n ||B||_1) = ", relOrthogError);

        // TODO: More refined failure conditions
        if( relOrthogError > Real(200) )
            LogicError("Relative orthogonality error was unacceptable");
        if( relError > Real(100) )
            LogicError("Relative error was unacceptably large");
    }
}

template<typename F>
void TestHermitianGenDefEigSequential
( Int m,
  UpperOrLower uplo,
  Pencil pencil,
  bool onlyEigvals,
  bool correctness,
  bool print,
  const HermitianEigCtrl<F>& ctrl )
{
    typedef Base<F> Real;
    Matrix<F> A, B, AOrig, BOrig;
    Matrix<Real> w;
    Matrix<F> X;
    Output("Testing with ",TypeName<F>());
    PushIndent();

    HermitianUniformSpectrum( A, m, 1, 10 );
    if( pencil == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        Zeros( B, m, m );
        Matrix<F> C;
        Uniform( C, m, m );
        Herk( uplo, ADJOINT, Real(1), C, Real(0), B );
    }
    else
        HermitianUniformSpectrum( B, m, 1, 10 );

    if( correctness && !onlyEigvals )
    {
        AOrig = A;
        BOrig = B;
    }
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
    }

    Output("Starting Hermitian Generalized-Definite Eigensolver...");
    Timer timer;
    timer.Start();
    if( onlyEigvals )
        HermitianGenDefEig( pencil, uplo, A, B, w, ctrl );
    else
        HermitianGenDefEig( pencil, uplo, A, B, w, X, ctrl );
    const double runTime = timer.Stop();
    Output("Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( X, "eigenvectors:" );
    }
    if( correctness && !onlyEigvals )
        TestCorrectness( print, pencil, uplo, AOrig, BOrig, A, B, w, X );
    PopIndent();
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=VR>
void TestHermitianGenDefEig
( Int m,
  UpperOrLower uplo,
  Pencil pencil,
  bool onlyEigvals,
  bool correctness,
  bool print,
  const Grid& g, 
  const HermitianEigCtrl<F>& ctrl )
{
    typedef Base<F> Real;
    DistMatrix<F,U,V> A(g), B(g), AOrig(g), BOrig(g);
    DistMatrix<Real,S,STAR> w(g);
    DistMatrix<F> X(g);
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

    HermitianUniformSpectrum( A, m, 1, 10 );
    if( pencil == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        Zeros( B, m, m );
        DistMatrix<F> C(g);
        Uniform( C, m, m );
        Herk( uplo, ADJOINT, Real(1), C, Real(0), B );
    }
    else
        HermitianUniformSpectrum( B, m, 1, 10 );

    if( correctness && !onlyEigvals )
    {
        AOrig = A;
        BOrig = B;
    }
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
    }

    OutputFromRoot
    (g.Comm(),"Starting Hermitian Generalized-Definite Eigensolver...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    if( onlyEigvals )
        HermitianGenDefEig( pencil, uplo, A, B, w, ctrl );
    else
        HermitianGenDefEig( pencil, uplo, A, B, w, X, ctrl );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( X, "eigenvectors:" );
    }
    if( correctness && !onlyEigvals )
        TestCorrectness( print, pencil, uplo, AOrig, BOrig, A, B, w, X );
    PopIndent();
}

template<typename F>
void TestSuite
( Int m,
  UpperOrLower uplo,
  Pencil pencil,
  bool onlyEigvals,
  bool sequential,
  bool distributed,
  bool correctness,
  bool print,
  const Grid& g,
  const HermitianEigCtrl<double>& ctrlDbl )
{
    typedef Base<F> Real;
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
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
    ctrl.tridiagEigCtrl.useQR = ctrlDbl.tridiagEigCtrl.useQR;
    ctrl.tridiagEigCtrl.subset = subset;

    if( sequential )
    {
        TestHermitianGenDefEigSequential<F>
        ( m, uplo, pencil, onlyEigvals, correctness, print, ctrl );
    }
    // Distributed eigensolvers are currently only supported for {float,double}
    const bool canSolve = IsBlasScalar<Real>::value || (g.Size()==1);
    if( distributed && canSolve )
    {
        OutputFromRoot(g.Comm(),"Normal tridiag algorithms:");
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
        TestHermitianGenDefEig<F>
        ( m, uplo, pencil, onlyEigvals, correctness, print, g, ctrl );

        OutputFromRoot(g.Comm(),"Square row-major algorithms:");
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl.tridiagCtrl.order = ROW_MAJOR;
        TestHermitianGenDefEig<F>
        ( m, uplo, pencil, onlyEigvals, correctness, print, g, ctrl );

        OutputFromRoot(g.Comm(),"Square column-major algorithms:");
        ctrl.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl.tridiagCtrl.order = COLUMN_MAJOR;
        TestHermitianGenDefEig<F>
        ( m, uplo, pencil, onlyEigvals, correctness, print, g, ctrl );

        // Also test with non-standard distributions
        OutputFromRoot(g.Comm(),"Nonstandard distributions:");
        TestHermitianGenDefEig<F,MR,MC,MC>
        ( m, uplo, pencil, onlyEigvals, correctness, print, g, ctrl );
    }

    PopIndent();
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int pencilInt = Input("--pencil",
             "1 is A x = lambda B x, "
             "2 is A B x = lambda x, "
             "3 is B A x = lambda x",1);
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
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool avoidTrmv = 
            Input("--avoidTrmv","avoid Trmv based Symv",true);
        const bool useScaLAPACK =
          Input("--useScaLAPACK","test ScaLAPACK?",false);
        const bool useQR = Input("--useQR","use QR algorithm?",false);
        const bool sequential =
          Input("--sequential","test sequential?",true);
        const bool distributed =
          Input("--distributed","test distributed?",true);
        const bool correctness = Input
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
            throw runtime_error("'range' must be 'A', 'I', or 'V'");
        const SortType sort = static_cast<SortType>(sortInt);
        if( correctness && onlyEigvals )
            OutputFromRoot
            (g.Comm(),"Cannot test correctness with only eigenvalues.");
        Pencil pencil;
        string pencilString;
        if( pencilInt == 1 )
        {
            pencil = AXBX;
            pencilString = "AXBX";
        }
        else if( pencilInt == 2 )
        {
            pencil = ABX;
            pencilString = "ABX";
        }
        else if( pencilInt == 3 )
        {
            pencil = BAX;
            pencilString = "BAX";
        }
        else
            LogicError("Invalid pencil integer");
        ComplainIfDebug();
        OutputFromRoot
        (g.Comm(),"Will test ",( uplo==LOWER ? "lower" : "upper" )," ",
         pencilString," HermitianGenDefEig.");

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
        ctrl.tridiagEigCtrl.useQR = useQR;
        ctrl.tridiagEigCtrl.subset = subset;

        if( testReal )
        {
            TestSuite<float>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<double>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );

#ifdef EL_HAVE_QD
            TestSuite<DoubleDouble>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<QuadDouble>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );
#endif

#ifdef EL_HAVE_QUAD
            TestSuite<Quad>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );
#endif

#ifdef EL_HAVE_MPC
            TestSuite<BigFloat>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );
#endif
        }
        if( testCpx )
        {
            TestSuite<Complex<float>>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<Complex<double>>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );

#ifdef EL_HAVE_QD
            TestSuite<Complex<DoubleDouble>>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );

            TestSuite<Complex<QuadDouble>>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );
#endif

#ifdef EL_HAVE_QUAD
            TestSuite<Complex<Quad>>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );
#endif

#ifdef EL_HAVE_MPC
            TestSuite<Complex<BigFloat>>
            ( m, uplo, pencil, onlyEigvals,
              sequential, distributed, correctness, print, g, ctrl );
#endif
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
