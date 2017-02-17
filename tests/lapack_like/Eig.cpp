/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2016, Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Real>
void TestCorrectness
( bool print,
  const Matrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& w,
  const Matrix<Complex<Real>>& V )
{
    const Int n = A.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( A );

    // Find the residual R = AV-VW
    Matrix<Complex<Real>> R( V.Height(), V.Width() );
    Gemm( NORMAL, NORMAL, Complex<Real>(1), A, V, R );
    Matrix<Complex<Real>> VW( V );
    DiagonalScale( RIGHT, NORMAL, w, VW );
    R -= VW;

    const Real infError = InfinityNorm( R );
    const Real relError = infError / (eps*n*oneNormA);
    Output("|| A V - V W ||_oo / (eps n || A ||_1) = ",relError);

    // TODO(poulson): A more refined failure condition
    if( relError > Real(100) )
        LogicError("Relative error was unacceptably large");
}

template<typename Real>
void TestCorrectness
( bool print,
  const ElementalMatrix<Complex<Real>>& A,
  const ElementalMatrix<Complex<Real>>& w,
  const ElementalMatrix<Complex<Real>>& V )
{
    const Int n = A.Height();
    const Real eps = limits::Epsilon<Real>();
    const Real oneNormA = OneNorm( A );

    // Find the residual R = AV-VW
    DistMatrix<Complex<Real>> R( V.Height(), V.Width(), A.Grid() );
    Gemm( NORMAL, NORMAL, Complex<Real>(1), A, V, R );
    DistMatrix<Complex<Real>> VW( V );
    DiagonalScale( RIGHT, NORMAL, w, VW );
    R -= VW;

    const Real infError = InfinityNorm( R );
    const Real relError = infError / (eps*n*oneNormA);
    OutputFromRoot
    (A.Grid().Comm(),"|| A V - V W ||_oo / (eps n || A ||_1) = ",relError);

    // TODO(poulson): A more refined failure condition
    if( relError > Real(100) )
        LogicError("Relative error was unacceptably large");
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void LAPACKHelper
( const Matrix<Complex<Real>>& AOrig,
  bool correctness,
  bool print )
{
    Output( "Testing LAPACK with ", TypeName<Complex<Real>>() );
    const Int m = AOrig.Height();
    Timer timer;

    auto A( AOrig );
    Matrix<Complex<Real>> w(m,1), V(m,m), X, tau;

    // Compute eigenvectors with LAPACK (GEHRD, HSEQR, TREVC, TRMM)
    Output("LAPACK (GEHRD, UNGHR, HSEQR, TREVC)");
    PushIndent();
    A = AOrig;
    tau.Resize( m, 1 );
    timer.Reset();
    Output("Transforming to upper Hessenberg form...");
    PushIndent();
    timer.Start();
    lapack::Hessenberg( m, A.Buffer(), A.LDim(), tau.Buffer() );
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    Output("Obtaining orthogonal matrix...");
    PushIndent();
    timer.Start();
    V = A;
    lapack::HessenbergGenerateUnitary( m, V.Buffer(), V.LDim(), tau.Buffer() );
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    Output("Schur decomposition...");
    PushIndent();
    timer.Start();
    {
        bool fullTriangle=true;
        bool multiplyQ=true;
        lapack::HessenbergSchur
        ( m,
          A.Buffer(), A.LDim(),
          w.Buffer(),
          V.Buffer(), V.LDim(),
          fullTriangle, multiplyQ );
    }
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    Output("Triangular eigensolver...");
    PushIndent();
    timer.Start();
    {
        bool accumulate=true;
        lapack::TriangEig
        ( m, A.Buffer(), A.LDim(), V.Buffer(), V.LDim(), accumulate );
    }
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    Output("Total Time = ",timer.Total()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( print, AOrig, w, V );
    PopIndent();

    // Compute eigenvectors with LAPACK (GEEV)
    Output("LAPACK (GEEV)");
    PushIndent();
    A = AOrig;
    timer.Reset();
    timer.Start();
    lapack::Eig
    ( m,
      A.Buffer(), A.LDim(),
      w.Buffer(),
      V.Buffer(), V.LDim() );
    Output("Total Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( print, AOrig, w, V );
    PopIndent();
}

template<typename Real>
void ElementalHelper
( const Matrix<Complex<Real>>& AOrig,
  bool correctness,
  bool print )
{
    Output( "Testing Elemental with ", TypeName<Complex<Real>>() );
    Timer timer;

    auto A( AOrig );
    Matrix<Complex<Real>> w, V, X;

    PushIndent();
    Output("Schur decomposition...");
    PushIndent();
    timer.Start();
    SchurCtrl<Real> schurCtrl;
    schurCtrl.time = true;
    //schurCtrl.hessSchurCtrl.progress = true;
    Schur( A, w, V, schurCtrl );
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    if( print )
    {
        Print( A, "T" );
        Print( w, "w" );
        Print( V, "Q" );
    }
    Output("Triangular eigensolver...");
    PushIndent();
    timer.Start();
    TriangEig( A, X );
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    if( print )
        Print( X, "X" );
    Output("Transforming to get eigenvectors...");
    PushIndent();
    timer.Start();
    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, Complex<Real>(1), X, V );
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    Output("Total Time = ",timer.Total()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( print, AOrig, w, V );
    PopIndent();
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void EigBenchmark
( Int m,
  bool correctness,
  bool print,
  Int whichMatrix )
{
    Output( "Testing with ", TypeName<Complex<Real>>() );
    Matrix<Complex<Real>> A(m,m), AOrig(m,m);
    Matrix<Complex<Real>> w(m,1), V(m,m);
    Matrix<Complex<Real>> X, tau;

    const Real foxLiOmega = 16*M_PI;

    // Generate test matrix
    switch( whichMatrix )
    {
    case 0: Gaussian( AOrig, m, m );       break;
    case 1: FoxLi( AOrig, m, foxLiOmega ); break;
    case 2: Grcar( AOrig, m );             break;
    case 3: Jordan( AOrig, m, Complex<Real>(7) ); break;
    default: LogicError("Unknown test matrix");
    }
    if( print )
        Print( AOrig, "A" );

    ElementalHelper( AOrig, correctness, print );
    LAPACKHelper( AOrig, correctness, print );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void EigBenchmark
( Int m,
  bool correctness,
  bool print,
  Int whichMatrix )
{
    Output( "Testing with ", TypeName<Complex<Real>>() );
    Matrix<Complex<Real>> A(m,m), AOrig(m,m);
    Matrix<Complex<Real>> w(m,1), V(m,m);
    Matrix<Complex<Real>> X, tau;

    const Real foxLiOmega = 16*M_PI;

    // Generate test matrix
    switch( whichMatrix )
    {
    case 0: Gaussian( AOrig, m, m );       break;
    case 1: FoxLi( AOrig, m, foxLiOmega ); break;
    case 2: Grcar( AOrig, m );             break;
    case 3: Jordan( AOrig, m, Complex<Real>(7) ); break;
    default: LogicError("Unknown test matrix");
    }
    if( print )
        Print( AOrig, "A" );

    ElementalHelper( AOrig, correctness, print );
}

template<typename Real>
void EigBenchmark
( const Grid& g,
  Int m,
  bool correctness,
  bool print,
  Int whichMatrix,
  bool scalapack,
  Int blockHeight )
{
    OutputFromRoot( g.Comm(), "Testing with ", TypeName<Complex<Real>>() );
    // TODO(poulson): Convert to distributed analogue
    DistMatrix<Complex<Real>> A(m,m,g), AOrig(m,m,g);
    DistMatrix<Complex<Real>> w(m,1,g), V(m,m,g), X(g);

    const Real foxLiOmega = 16*M_PI;

    // Generate test matrix
    switch( whichMatrix )
    {
    case 0: Gaussian( AOrig, m, m );       break;
    case 1: FoxLi( AOrig, m, foxLiOmega ); break;
    case 2: Grcar( AOrig, m );             break;
    case 3: Jordan( AOrig, m, Complex<Real>(7) ); break;
    default: LogicError("Unknown test matrix");
    }
    if( print )
        Print( AOrig, "A" );

    Timer timer;

    SchurCtrl<Real> schurCtrl;
    schurCtrl.hessSchurCtrl.fullTriangle = true;
    schurCtrl.hessSchurCtrl.scalapack = scalapack;
    schurCtrl.hessSchurCtrl.blockHeight = blockHeight;
    schurCtrl.time = true;

    // Compute eigenvectors with Elemental
    OutputFromRoot(g.Comm(),"Elemental");
    PushIndent();
    A = AOrig;
    OutputFromRoot(g.Comm(),"Schur decomposition...");
    PushIndent();
    timer.Start();
    Schur( A, w, V, schurCtrl );
    OutputFromRoot(g.Comm(),"Time = ",timer.Stop()," seconds");
    PopIndent();
    if( print )
    {
        Print( A, "T" );
        Print( V, "Q" );
    }
    OutputFromRoot(g.Comm(),"Triangular eigensolver...");
    PushIndent();
    timer.Start();
    TriangEig( A, X );
    OutputFromRoot(g.Comm(),"Time = ",timer.Stop()," seconds");
    PopIndent();
    if( print )
        Print( X, "X" );
    OutputFromRoot(g.Comm(),"Transforming to get eigenvectors...");
    PushIndent();
    timer.Start();
    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, Complex<Real>(1), X, V );
    OutputFromRoot(g.Comm(),"Time = ",timer.Stop()," seconds");
    PopIndent();
    OutputFromRoot(g.Comm(),"Total Time = ",timer.Total()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( correctness )
        TestCorrectness( print, AOrig, w, V );
    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        // Parse command line arguments
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int blockHeight =
          Input("--blockHeight","ScaLAPACK block height",32);
        // NOTE: Distributed AED is not supported by ScaLAPACK for complex :-(
        const bool scalapack =
          Input
          ("--scalapack","Use ScaLAPACK for the distributed solver?",false);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool distributed =
          Input("--distributed","test distributed?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const Int whichMatrix =
          Input("--whichMatrix","(0=Gaussian,1=Fox-Li,2=Grcar,3=Jordan)",0);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        ComplainIfDebug();

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = colMajor ? COLUMN_MAJOR : ROW_MAJOR;
        const Grid grid( comm, gridHeight, order );

        if( sequential && commRank == 0 )
        {
            EigBenchmark<float>( n, correctness, print, whichMatrix );
            EigBenchmark<double>( n, correctness, print, whichMatrix );
#ifdef EL_HAVE_QD
            EigBenchmark<DoubleDouble>( n, correctness, print, whichMatrix );
            EigBenchmark<QuadDouble>( n, correctness, print, whichMatrix );
#endif
#ifdef EL_HAVE_QUAD
            EigBenchmark<Quad>( n, correctness, print, whichMatrix );
#endif
#ifdef EL_HAVE_MPC
            EigBenchmark<BigFloat>( n, correctness, print, whichMatrix );
#endif
        }
        if( distributed )
        {
            EigBenchmark<float>
            ( grid, n, correctness, print, whichMatrix, scalapack,
              blockHeight );
            EigBenchmark<double>
            ( grid, n, correctness, print, whichMatrix, scalapack,
              blockHeight );
#ifdef EL_HAVE_QD
            EigBenchmark<DoubleDouble>
            ( grid, n, correctness, print, whichMatrix, scalapack,
              blockHeight );
            EigBenchmark<QuadDouble>
            ( grid, n, correctness, print, whichMatrix, scalapack,
              blockHeight );
#endif
#ifdef EL_HAVE_QUAD
            EigBenchmark<Quad>
            ( grid, n, correctness, print, whichMatrix, scalapack,
              blockHeight );
#endif
#ifdef EL_HAVE_MPC
            EigBenchmark<BigFloat>
            ( grid, n, correctness, print, whichMatrix, scalapack,
              blockHeight );
#endif
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
