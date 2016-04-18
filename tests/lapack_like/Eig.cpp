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

void TestCorrectness
( bool print,
  const Matrix<Complex<double>>& A,
  const Matrix<Complex<double>>& w,
  const Matrix<Complex<double>>& V )
{
    // Find the residual R = AV-VW
    Matrix<Complex<double>> R( V.Height(), V.Width() );
    Gemm
    ( NORMAL, NORMAL,
      Complex<double>(1), A, V,
      Complex<double>(0), R);
    Matrix<Complex<double>> VW( V );
    DiagonalScale( RIGHT, NORMAL, w, VW );
    R -= VW;
    
    // Find the Frobenius norms of A and AV-VW
    double frobNormA = FrobeniusNorm( A );
    double frobNormR = FrobeniusNorm( R );

    Output("  ||A V - V W||_F / ||A||_F = ",frobNormR/frobNormA);
}

void TestCorrectness
( bool print,
  const ElementalMatrix<Complex<double>>& A,
  const ElementalMatrix<Complex<double>>& w,
  const ElementalMatrix<Complex<double>>& V )
{
    // Find the residual R = AV-VW
    DistMatrix<Complex<double>> R( V.Height(), V.Width(), A.Grid() );
    Gemm
    ( NORMAL, NORMAL,
      Complex<double>(1), A, V,
      Complex<double>(0), R);
    DistMatrix<Complex<double>> VW( V );
    DiagonalScale( RIGHT, NORMAL, w, VW );
    R -= VW;
    
    // Find the Frobenius norms of A and AV-VW
    double frobNormA = FrobeniusNorm( A );
    double frobNormR = FrobeniusNorm( R );

    if( A.Grid().Rank() == 0 )
        Output("  ||A V - V W||_F / ||A||_F = ",frobNormR/frobNormA);
}

void SequentialEigBenchmark
( bool testCorrectness,
  bool print,
  Int m,
  Int testMatrix )
{
    Matrix<Complex<double>> A(m,m), AOrig(m,m);
    Matrix<Complex<double>> w(m,1), V(m,m);
    Matrix<Complex<double>> X, tau;
    
    //double foxLiOmega = -0.179;
    double foxLiOmega = 16*M_PI;

    // Generate test matrix
    switch( testMatrix )
    {
    case 0: Gaussian( AOrig, m, m );       break;
    case 1: FoxLi( AOrig, m, foxLiOmega ); break;
    case 2: Grcar( AOrig, m );             break;
    default: LogicError("Unknown test matrix");
    }
    if( print )
        Print( AOrig, "A" );
 
    Timer timer;

    SchurCtrl<double> schurCtrl;
    schurCtrl.time = true;

    // Compute eigenvectors with Elemental
    Output("Elemental");
    A = AOrig;
    Output("  Schur decomposition...");
    timer.Start();
    Schur( A, w, V, true, schurCtrl );
    Output("    Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( A, "T" );
        Print( V, "Q" );
    }
    Output("  Triangular eigensolver...");
    timer.Start();
    TriangEig( A, X );
    Output("    Time = ",timer.Stop()," seconds");
    if( print ) 
        Print( X, "X" );
    Output("  Transforming to get eigenvectors...");
    timer.Start();
    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, Complex<double>(1), X, V );
    Output("    Time = ",timer.Stop()," seconds");
    Output("  Total Time = ",timer.Total()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( testCorrectness )
        TestCorrectness( print, AOrig, w, V );
    
    // Compute eigenvectors with LAPACK (GEHRD, HSEQR, TREVC, TRMM)
    Output("LAPACK (GEHRD, UNGHR, HSEQR, TREVC)");
    A = AOrig;
    tau.Resize( m, 1 );
    timer.Reset();
    Output("  Transforming to upper Hessenberg form...");
    timer.Start();
    lapack::Hessenberg( m, A.Buffer(), A.LDim(), tau.Buffer() );
    Output("    Time = ",timer.Stop()," seconds");
    Output("  Obtaining orthogonal matrix...");
    timer.Start();
    V = A;
    lapack::HessenbergGenerateUnitary( m, V.Buffer(), V.LDim(), tau.Buffer() );
    Output("    Time = ",timer.Stop()," seconds");
    Output("  Schur decomposition...");
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
    Output("    Time = ",timer.Stop()," seconds");
    Output("  Triangular eigensolver...");
    timer.Start();
    {
        bool accumulate=true;
        lapack::TriangEig
        ( m, A.Buffer(), A.LDim(), V.Buffer(), V.LDim(), accumulate );
    }
    Output("    Time = ",timer.Stop()," seconds");
    Output("  Total Time = ",timer.Total()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( testCorrectness )
        TestCorrectness( print, AOrig, w, V );

    // Compute eigenvectors with LAPACK (GEEV)
    Output("LAPACK (GEEV)");
    A = AOrig;
    timer.Reset();
    timer.Start();
    lapack::Eig
    ( m, 
      A.Buffer(), A.LDim(),
      w.Buffer(),
      V.Buffer(), V.LDim() );
    Output("  Total Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( testCorrectness )
        TestCorrectness( print, AOrig, w, V );
}

void EigBenchmark
( bool testCorrectness,
  bool print,
  Int m,
  Int testMatrix,
  bool distAED,
  Int blockHeight,
  const Grid& g )
{
    const int gridRank = g.Rank();

    // TODO: Convert to distributed analogue
    DistMatrix<Complex<double>> A(m,m,g), AOrig(m,m,g);
    DistMatrix<Complex<double>> w(m,1,g), V(m,m,g), X(g);

    //double foxLiOmega = -0.179;
    double foxLiOmega = 16*M_PI;
    
    // Generate test matrix
    switch( testMatrix )
    {
    case 0: Gaussian( AOrig, m, m );       break;
    case 1: FoxLi( AOrig, m, foxLiOmega ); break;
    case 2: Grcar( AOrig, m );             break;
    default: LogicError("Unknown test matrix");
    }
    if( print )
        Print( AOrig, "A" );

    Timer timer;
 
    SchurCtrl<double> schurCtrl;
    schurCtrl.qrCtrl.distAED = distAED;
    schurCtrl.qrCtrl.blockHeight = blockHeight;
    schurCtrl.time = true;

    // Compute eigenvectors with Elemental
    if( gridRank == 0 )
        Output("Elemental");
    A = AOrig;
    if( gridRank == 0 )
    {
        Output("  Schur decomposition...");
        timer.Start();
    }
    Schur( A, w, V, true, schurCtrl );
    if( gridRank == 0 )
        Output("    Time = ",timer.Stop()," seconds");
    if( print )
    {
        Print( A, "T" );
        Print( V, "Q" );
    }
    if( gridRank == 0 )
    {
        Output("  Triangular eigensolver...");
        timer.Start();
    }
    TriangEig( A, X );
    if( gridRank == 0 )
        Output("    Time = ",timer.Stop()," seconds");
    if( print )
        Print( X, "X" );
    if( gridRank == 0 )
    {
        Output("  Transforming to get eigenvectors...");
        timer.Start();
    }
    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT, Complex<double>(1), X, V );
    if( gridRank == 0 )
    {
        Output("    Time = ",timer.Stop()," seconds");
        Output("  Total Time = ",timer.Total()," seconds");
    }
    if( print )
    {
        Print( w, "eigenvalues:" );
        Print( V, "eigenvectors:" );
    }
    if( testCorrectness )
        TestCorrectness( print, AOrig, w, V );
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
        // Parse command line arguments
        int r = Input("--gridHeight","height of process grid",0);
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
           true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",false);
        const bool print = Input("--print","print matrices?",false);
        const Int testMatrix =
          Input("--testMatrix","test matrix (0=Gaussian,1=Fox-Li,2=Grcar)",0);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        ComplainIfDebug();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR ); 
        const Grid grid( comm, r, order );
       
        if( commRank == 0 )
            SequentialEigBenchmark( testCorrectness, print, n, testMatrix );

        EigBenchmark
        ( testCorrectness, print, n, testMatrix,
          distAED, blockHeight, grid );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
