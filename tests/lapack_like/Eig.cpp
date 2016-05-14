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

    // Find the residual R = AV-VW
    Matrix<Complex<Real>> R( V.Height(), V.Width() );
    Gemm
    ( NORMAL, NORMAL,
      Complex<Real>(1), A, V,
      Complex<Real>(0), R);
    Matrix<Complex<Real>> VW( V );
    DiagonalScale( RIGHT, NORMAL, w, VW );
    R -= VW;
    
    // Find the Frobenius norms of A and AV-VW
    const Real frobNormA = FrobeniusNorm( A );
    const Real frobNormR = FrobeniusNorm( R );
    const Real relError = frobNormR / (eps*n*frobNormA); 
    Output("|| A V - V W ||_F / (eps n || A ||_F) = ",relError);

    // TODO: A more refined failure condition
    if( relError > Real(10) )
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

    // Find the residual R = AV-VW
    DistMatrix<Complex<Real>> R( V.Height(), V.Width(), A.Grid() );
    Gemm
    ( NORMAL, NORMAL,
      Complex<Real>(1), A, V,
      Complex<Real>(0), R);
    DistMatrix<Complex<Real>> VW( V );
    DiagonalScale( RIGHT, NORMAL, w, VW );
    R -= VW;
    
    // Find the Frobenius norms of A and AV-VW
    const Real frobNormA = FrobeniusNorm( A );
    const Real frobNormR = FrobeniusNorm( R );
    const Real relError = frobNormR / (eps*n*frobNormA); 
    OutputFromRoot
    (A.Grid().Comm(),"|| A V - V W ||_F / (eps || A ||_F) = ",relError);

    // TODO: A more refined failure condition
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename Real>
void SequentialEigBenchmark
( bool testCorrectness,
  bool print,
  Int m,
  Int testMatrix )
{
    Output( "Testing with ", TypeName<Real>() );
    Matrix<Complex<Real>> A(m,m), AOrig(m,m);
    Matrix<Complex<Real>> w(m,1), V(m,m);
    Matrix<Complex<Real>> X, tau;
    
    const Real foxLiOmega = 16*M_PI;

    // Generate test matrix
    switch( testMatrix )
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
    schurCtrl.time = true;

    // Compute eigenvectors with Elemental
    Output("Elemental");
    PushIndent();
    A = AOrig;
    Output("Schur decomposition...");
    PushIndent();
    timer.Start();
    Schur( A, w, V, true, schurCtrl );
    Output("Time = ",timer.Stop()," seconds");
    PopIndent();
    if( print )
    {
        Print( A, "T" );
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
    if( testCorrectness )
        TestCorrectness( print, AOrig, w, V );
    PopIndent();
    
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
    if( testCorrectness )
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
    if( testCorrectness )
        TestCorrectness( print, AOrig, w, V );
    PopIndent();
}

template<typename Real>
void EigBenchmark
( bool testCorrectness,
  bool print,
  Int m,
  Int testMatrix,
  bool distAED,
  Int blockHeight,
  const Grid& g )
{
    OutputFromRoot( g.Comm(), "Testing with ", TypeName<Real>() );
    // TODO: Convert to distributed analogue
    DistMatrix<Complex<Real>> A(m,m,g), AOrig(m,m,g);
    DistMatrix<Complex<Real>> w(m,1,g), V(m,m,g), X(g);

    const Real foxLiOmega = 16*M_PI;
    
    // Generate test matrix
    switch( testMatrix )
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
    schurCtrl.qrCtrl.distAED = distAED;
    schurCtrl.qrCtrl.blockHeight = blockHeight;
    schurCtrl.time = true;

    // Compute eigenvectors with Elemental
    OutputFromRoot(g.Comm(),"Elemental");
    PushIndent();
    A = AOrig;
    OutputFromRoot(g.Comm(),"Schur decomposition...");
    PushIndent();
    timer.Start();
    Schur( A, w, V, true, schurCtrl );
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
    if( testCorrectness )
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
        const bool distAED =
          Input
          ("--distAED",
           "Distributed Aggressive Early Deflation? (it can be buggy...)",
           false);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const Int testMatrix =
          Input("--testMatrix","(0=Gaussian,1=Fox-Li,2=Grcar,3=Jordan)",0);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        ComplainIfDebug();

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR ); 
        const Grid grid( comm, gridHeight, order );
       
        if( commRank == 0 )
        {
            SequentialEigBenchmark<float>
            ( testCorrectness, print, n, testMatrix );
            SequentialEigBenchmark<double>
            ( testCorrectness, print, n, testMatrix );
        }

        EigBenchmark<float>
        ( testCorrectness, print, n, testMatrix,
          distAED, blockHeight, grid );
        EigBenchmark<double>
        ( testCorrectness, print, n, testMatrix,
          distAED, blockHeight, grid );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
