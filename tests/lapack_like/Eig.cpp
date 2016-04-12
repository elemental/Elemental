/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2016, Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

void TestCorrectness
( bool print,
  const Matrix<Complex<double>>& A,
  const Matrix<Complex<double>>& w,
  const Matrix<Complex<double>>& V)
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

    // Find condition number
    Output("  ||A V - V W||_F / ||A||_F = ",frobNormR/frobNormA);
}

void EigBenchmark
( bool testCorrectness,
  bool print,
  Int m,
  Int testMatrix )
{
    Matrix<Complex<double>> A(m,m), AOrig(m,m);
    Matrix<Complex<double>> w(m,1), V(m,m);
    Matrix<Complex<double>> work, tau;
    
    double time1, time2, time3, time4;
    
    // Generate test matrix
    switch( testMatrix )
    {

    case 0:
        // Gaussian matrix
        Gaussian( AOrig, m, m );
        break;

    case 1:
        // Fox-Li matrix
        FoxLi( AOrig, m, -0.179 );
        break;
        
    case 2:
        // Grcar matrix
        Grcar( AOrig, m );
        break;
        
    default:
        LogicError("Unknown test matrix");
        break;

    }
    if( print )
        Print( AOrig, "A" );
 
    SchurCtrl<double> schurCtrl;
    schurCtrl.time = true;

    // Compute eigenvectors with Elemental
    Output("Elemental");
    A = AOrig;
    time1 = mpi::Time();
    Output("  Schur decomposition...");
    time3 = mpi::Time();
    Schur( A, w, V, true, schurCtrl );
    time4 = mpi::Time();
    Output("    Time = ",time4-time3," seconds");
    Output("  Triangular eigensolver...");
    time3 = mpi::Time();
    TriangEig( A, work );
    time4 = mpi::Time();
    Output("    Time = ",time4-time3," seconds");
    Output("  Transforming to get eigenvectors...");
    time3 = mpi::Time();
    Trmm( RIGHT, UPPER, NORMAL, NON_UNIT,
          Complex<double>(1), work, V );
    time4 = mpi::Time();
    Output("    Time = ",time4-time3," seconds");
    time2 = mpi::Time();
    Output("  Total Time = ",time2-time1," seconds");
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
    time1 = mpi::Time();
    Output("  Transforming to upper Hessenberg form...");
    time3 = mpi::Time();
    lapack::Hessenberg( m, A.Buffer(), A.LDim(), tau.Buffer() );
    time4 = mpi::Time();
    Output("    Time = ",time4-time3," seconds");
    Output("  Obtaining orthogonal matrix...");
    time3 = mpi::Time();
    V = A;
    lapack::HessenbergGenerateUnitary( m, V.Buffer(), V.LDim(), tau.Buffer() );
    time4 = mpi::Time();
    Output("    Time = ",time4-time3," seconds");
    Output("  Schur decomposition...");
    time3 = mpi::Time();
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
    time4 = mpi::Time();
    Output("    Time = ",time4-time3," seconds");
    Output("  Triangular eigensolver...");
    time3 = mpi::Time();
    {
        bool accumulate=true;
        lapack::TriangEig
        ( m, A.Buffer(), A.LDim(), V.Buffer(), V.LDim(), accumulate );
    }
    time4 = mpi::Time();
    Output("    Time = ",time4-time3," seconds");
    time2 = mpi::Time();
    Output("  Total Time = ",time2-time1," seconds");
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
    time1 = mpi::Time();
    lapack::Eig
    ( m, 
      A.Buffer(), A.LDim(),
      w.Buffer(),
      V.Buffer(), V.LDim() );
    time2 = mpi::Time();
    Output("  Total Time = ",time2-time1," seconds");
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

    try
    {
        // Parse command line arguments
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",false);
        const bool print = Input("--print","print matrices?",false);
        const Int testMatrix = Input("--testMatrix","test matrix (0=Gaussian,1=Fox-Li,2=Grcar)",0);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        ComplainIfDebug();

        // Benchmark triangular eigensolver
        EigBenchmark(testCorrectness, print, n, testMatrix);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
