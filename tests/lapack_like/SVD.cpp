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
void TestSequentialSVD
( Int m, Int n, Int rank,
  SVDApproach approach,
  SingularValueToleranceType tolType,
  double tol,
  bool time,
  bool progress,
  bool wantU,
  bool wantV,
  bool useQR,
  bool penalizeDerivative,
  Int divideCutoff,
  bool print )
{
    typedef Base<F> Real;
    Timer timer;

    Output("Sequential test with ",TypeName<F>());

    Matrix<F> A;
    {
        Matrix<F> X, Y; 
        Uniform( X, m, rank );
        Uniform( Y, rank, n );
        Gemm( NORMAL, NORMAL, F(1), X, Y, A );
    }
    if( print )
        Print( A, "ASeq" );

    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.useQR = useQR;
    ctrl.bidiagSVDCtrl.wantU = wantU; 
    ctrl.bidiagSVDCtrl.wantV = wantV;
    ctrl.bidiagSVDCtrl.approach = approach;
    ctrl.bidiagSVDCtrl.tolType = tolType;
    ctrl.bidiagSVDCtrl.tol = tol;
    ctrl.bidiagSVDCtrl.progress = progress;
    ctrl.bidiagSVDCtrl.dcCtrl.cutoff = divideCutoff;
    ctrl.bidiagSVDCtrl.dcCtrl.secularCtrl.penalizeDerivative =
      penalizeDerivative;
    ctrl.bidiagSVDCtrl.dcCtrl.secularCtrl.progress = progress;
    ctrl.time = time;

    Matrix<Real> s;
    Matrix<F> U, V;
    timer.Start();
    auto info = SVD( A, U, s, V, ctrl );
    Output("Sequential SVD: ",timer.Stop());

    const auto& qrInfo = info.bidiagSVDInfo.qrInfo;
    const auto& dcInfo = info.bidiagSVDInfo.dcInfo;
    if( qrInfo.numIterations > 0 ) 
    {
        Output("  num QR Iterations: ",qrInfo.numIterations);
        Output("    numZeroShiftForward: ",
          qrInfo.numZeroShiftForwardIterations);
        Output("    numZeroShiftBackward: ",
          qrInfo.numZeroShiftForwardIterations);
        Output("    numNonzeroShiftForward: ",
          qrInfo.numNonzeroShiftForwardIterations);
        Output("    numNonzeroShiftBackward: ",
          qrInfo.numZeroShiftForwardIterations);
        Output("  numInnerLoops: ",qrInfo.numInnerLoops); 
        Output("    numZeroShiftForward: ",
          qrInfo.numZeroShiftForwardInnerLoops);
        Output("    numZeroShiftBackward: ",
          qrInfo.numZeroShiftForwardInnerLoops);
        Output("    numNonzeroShiftForward: ",
          qrInfo.numNonzeroShiftForwardInnerLoops);
        Output("    numNonzeroShiftBackward: ",
          qrInfo.numZeroShiftForwardInnerLoops);
    }
    else
    {
        const auto& secularInfo = dcInfo.secularInfo;
        Output("  num D&C deflations: ",secularInfo.numDeflations); 
        Output("    small diagonal: ",secularInfo.numSmallDiagonalDeflations);
        Output("    close diagonal: ",secularInfo.numCloseDiagonalDeflations);
        Output("    small update;   ",secularInfo.numSmallUpdateDeflations);
        Output("  num secular iterations: ",secularInfo.numIterations);
        Output("  num secular alternations: ",secularInfo.numAlternations);
        Output("  num secular cubic iter's: ",secularInfo.numCubicIterations);
        Output("  num secular cubic failures: ",secularInfo.numCubicFailures);
    }
    if( print )
    {
        if( wantU )
            Print( U, "USeq" );
        Print( s, "sSeq" );
        if( wantV )
            Print( V, "VSeq" );
    }
    // Check that U and V are unitary
    Matrix<F> E;
    if( wantU )
    {
        Identity( E, U.Width(), U.Width() );
        Herk( LOWER, ADJOINT, Real(-1), U, Real(1), E );
        const Real UOrthErr = HermitianMaxNorm( LOWER, E );
        Output("|| I - U^H U ||_max = ",UOrthErr);
    }
    if( wantV )
    {
        Identity( E, V.Width(), V.Width() );
        Herk( LOWER, ADJOINT, Real(-1), V, Real(1), E );
        const Real VOrthErr = HermitianMaxNorm( LOWER, E );
        Output("|| I - V^H V ||_max = ",VOrthErr);
    }

    // Compute the residual error
    const Real twoNormA = MaxNorm( s );
    const Real maxNormA = MaxNorm( A );
    Output("|| A ||_max   = ",maxNormA);
    Output("|| A ||_2     = ",twoNormA);
    const Int numSingVals = s.Height();
    if( wantU && wantV )
    {
        auto UL = U( ALL, IR(0,numSingVals) );
        auto VL = V( ALL, IR(0,numSingVals) );
        DiagonalScale( RIGHT, NORMAL, s, UL );
        E = A;
        Gemm( NORMAL, ADJOINT, F(-1), UL, VL, F(1), E );
        if( print )
            Print( E, "A - U S V'" );
        const Real maxNormE = MaxNorm( E );
        const Real frobNormE = FrobeniusNorm( E );
        const Real eps = limits::Epsilon<Real>();
        const Real scaledResidual = frobNormE / (Max(m,n)*eps*twoNormA);
        Output("||A - U Sigma V^H||_max = ",maxNormE);
        Output("||A - U Sigma V^H||_F   = ",frobNormE);
        Output
        ("||A - U Sigma V_H||_F / (max(m,n) eps ||A||_2) = ",scaledResidual);
        // TODO(poulson): Provide a rigorous motivation for this bound
        if( scaledResidual > Real(50) )
            LogicError("SVD residual was unacceptably large");
    }
    Output("");
}

template<typename F>
void TestDistributedSVD
( Int m, Int n, Int rank,
  SVDApproach approach,
  SingularValueToleranceType tolType,
  double tol,
  bool time,
  bool progress,
  bool scalapack,
  bool wantU,
  bool wantV,
  bool useQR,
  bool penalizeDerivative,
  Int divideCutoff,
  bool print )
{
    typedef Base<F> Real;
    const int commRank = mpi::Rank();
    Timer timer;

    Grid g( mpi::COMM_WORLD );
    if( commRank == 0 )
        Output("Grid is ",g.Height()," x ",g.Width());
    DistMatrix<F> A(g), X(g), Y(g);
    Uniform( X, m, rank );
    Uniform( Y, rank, n );
    Gemm( NORMAL, NORMAL, F(1), X, Y, A );
    if( print )
        Print( A, "A" );

    // Compute the SVD of A 
    SVDCtrl<Real> ctrl;
    ctrl.bidiagSVDCtrl.useQR = useQR;
    ctrl.bidiagSVDCtrl.wantU = wantU; 
    ctrl.bidiagSVDCtrl.wantV = wantV;
    ctrl.bidiagSVDCtrl.approach = approach;
    ctrl.bidiagSVDCtrl.tolType = tolType;
    ctrl.bidiagSVDCtrl.tol = tol;
    ctrl.bidiagSVDCtrl.progress = progress;
    ctrl.bidiagSVDCtrl.dcCtrl.cutoff = divideCutoff;
    ctrl.bidiagSVDCtrl.dcCtrl.secularCtrl.penalizeDerivative =
      penalizeDerivative;
    ctrl.bidiagSVDCtrl.dcCtrl.secularCtrl.progress = progress;
    ctrl.time = time;

    ctrl.time = time;
    ctrl.useScaLAPACK = scalapack;
    if( commRank == 0 )
        timer.Start();
    DistMatrix<F> U(g), V(g);
    DistMatrix<Real,VR,STAR> s(g);
    auto info = SVD( A, U, s, V, ctrl );
    if( commRank == 0 )
    {
        Output("  SVD time: ",timer.Stop());
        const auto& qrInfo = info.bidiagSVDInfo.qrInfo;
        const auto& dcInfo = info.bidiagSVDInfo.dcInfo;
        if( qrInfo.numIterations > 0 ) 
        {
            Output("  num QR Iterations: ",qrInfo.numIterations);
            Output("    numZeroShiftForward: ",
              qrInfo.numZeroShiftForwardIterations);
            Output("    numZeroShiftBackward: ",
              qrInfo.numZeroShiftForwardIterations);
            Output("    numNonzeroShiftForward: ",
              qrInfo.numNonzeroShiftForwardIterations);
            Output("    numNonzeroShiftBackward: ",
              qrInfo.numZeroShiftForwardIterations);
            Output("  numInnerLoops: ",qrInfo.numInnerLoops); 
            Output("    numZeroShiftForward: ",
              qrInfo.numZeroShiftForwardInnerLoops);
            Output("    numZeroShiftBackward: ",
              qrInfo.numZeroShiftForwardInnerLoops);
            Output("    numNonzeroShiftForward: ",
              qrInfo.numNonzeroShiftForwardInnerLoops);
            Output("    numNonzeroShiftBackward: ",
              qrInfo.numZeroShiftForwardInnerLoops);
        }
        else
        {
            const auto& secularInfo = dcInfo.secularInfo;
            Output("  num D&C deflations: ",secularInfo.numDeflations); 
            Output
            ("    small diagonal: ",secularInfo.numSmallDiagonalDeflations);
            Output
            ("    close diagonal: ",secularInfo.numCloseDiagonalDeflations);
            Output
            ("    small update;   ",secularInfo.numSmallUpdateDeflations);
            Output
            ("  num secular iterations: ",secularInfo.numIterations);
            Output
            ("  num secular alternations: ",secularInfo.numAlternations);
            Output
            ("  num secular cubic iter's: ",secularInfo.numCubicIterations);
            Output
            ("  num secular cubic failures: ",secularInfo.numCubicFailures);
        }
    }
    if( print )
    {
        if( wantU )
            Print( U, "U" );
        Print( s, "s" );
        if( wantV )
            Print( V, "V" );
    }

    // Check that U and V are unitary
    DistMatrix<F> E(g);
    if( wantU )
    {
        Identity( E, U.Width(), U.Width() );
        Herk( LOWER, ADJOINT, Real(-1), U, Real(1), E );
        const Real UOrthErr = HermitianMaxNorm( LOWER, E );
        if( commRank == 0 )
            Output("|| I - U^H U ||_max = ",UOrthErr);
    }
    if( wantV )
    {
        Identity( E, V.Width(), V.Width() );
        Herk( LOWER, ADJOINT, Real(-1), V, Real(1), E );
        const Real VOrthErr = HermitianMaxNorm( LOWER, E );
        if( commRank == 0 )
            Output("|| I - V^H V ||_max = ",VOrthErr);
    }

    // Compute the residual error
    const Real twoNormA = MaxNorm( s );
    const Real maxNormA = MaxNorm( A );
    const Int numSingVals = s.Height();
    if( commRank == 0 )
    {
        Output("|| A ||_max   = ",maxNormA);
        Output("|| A ||_2     = ",twoNormA);
    }
    if( wantU && wantV )
    {
        auto UL = U( ALL, IR(0,numSingVals) );
        auto VL = V( ALL, IR(0,numSingVals) );
        DiagonalScale( RIGHT, NORMAL, s, UL );
        E = A;
        Gemm( NORMAL, ADJOINT, F(-1), UL, VL, F(1), E );
        if( print )
            Print( E, "A - U S V'" );
        const Real maxNormE = MaxNorm( E );
        const Real frobNormE = FrobeniusNorm( E );
        const Real eps = limits::Epsilon<Real>();
        const Real scaledResidual = frobNormE / (Max(m,n)*eps*twoNormA);
        if( commRank == 0 )
        {
            Output("||A - U Sigma V^H||_max = ",maxNormE);
            Output("||A - U Sigma V^H||_F   = ",frobNormE);
            Output
            ("||A - U Sigma V_H||_F / (max(m,n) eps ||A||_2) = ",
             scaledResidual);
        }
        // TODO(poulson): Provide a rigorous motivation for this bound
        if( scaledResidual > Real(50) )
            LogicError("SVD residual was unacceptably large");
    }
    if( commRank == 0 )
        Output("");
}

template<typename F>
void TestSVD
( Int m, Int n, Int rank,
  SVDApproach approach,
  SingularValueToleranceType tolType,
  double tol,
  bool time,
  bool progress,
  bool scalapack,
  bool testSeq,
  bool testDist,
  bool wantU,
  bool wantV,
  bool useQR,
  bool penalizeDerivative,
  Int divideCutoff,
  bool print )
{
    const int commRank = mpi::Rank();
    if( testSeq && commRank == 0 )
    {
        TestSequentialSVD<F>
        ( m, n, rank, approach, tolType, tol, time, progress, wantU, wantV,
          useQR, penalizeDerivative, divideCutoff, print );
    }
    if( testDist )
    {
        TestDistributedSVD<F> 
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          wantU, wantV, useQR, penalizeDerivative, divideCutoff, print );
    }
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int rank = Input("--rank","rank of matrix",10);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",32);
        const Int approachInt = Input("--approach","SVD approach",0);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",false);
        const Int mb = Input("--mb","block height",32);
        const Int nb = Input("--nb","block width",32);
#else
        const bool scalapack = false;
        const Int mb = 32;
        const Int nb = 32;
#endif
        const bool time = Input("--time","time SVD components?",true);
        const bool progress = Input("--progress","print progress?",false);
        const Int tolTypeInt = Input("--tolTypeInt","tolerance type int",1);
        const double tol = Input("--tol","threshold tol",double(0));

        const bool testSeq = Input("--testSeq","test sequential SVD?",true);
        const bool testDist = Input("--testDist","test distributed SVD?",true);
        const bool wantU = Input("--wantU","compute U?",true);
        const bool wantV = Input("--wantV","compute V?",true);
        const bool useQR = Input("--useQR","force use of QR algorithm?",false);
        const bool penalizeDerivative =
          Input
          ("--penalizeDerivative","penalize secular derivative in D&C?",false);
        const Int divideCutoff = Input("--divideCutoff","D&C cutoff?",60);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const SVDApproach approach = static_cast<SVDApproach>(approachInt);
        const SingularValueToleranceType tolType =
          static_cast<SingularValueToleranceType>(tolTypeInt);
        if( mpi::Rank() == 0 )
        {
            if( tolType == RELATIVE_TO_SELF_SING_VAL_TOL )
              Output("Testing with RELATIVE_TO_SELF_SING_VAL_TOL");
            else if( tolType == RELATIVE_TO_MAX_SING_VAL_TOL )
              Output("Testing with RELATIVE_TO_MAX_SING_VAL_TOL");
            else
              Output("Testing with ABSOLUTE_SING_VAL_TOL");
        }

        SetBlocksize( blocksize );
        SetDefaultBlockHeight( mb );
        SetDefaultBlockWidth( nb );

        TestSVD<float>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
        TestSVD<Complex<float>>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );

        TestSVD<double>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
        TestSVD<Complex<double>>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );

#ifdef EL_HAVE_QD
        TestSVD<DoubleDouble>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
        TestSVD<Complex<DoubleDouble>>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );

        TestSVD<QuadDouble>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
        TestSVD<Complex<QuadDouble>>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
#endif

#ifdef EL_HAVE_QUAD
        TestSVD<Quad>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
        TestSVD<Complex<Quad>>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
#endif

#ifdef EL_HAVE_MPC
        TestSVD<BigFloat>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
        TestSVD<Complex<BigFloat>>
        ( m, n, rank, approach, tolType, tol, time, progress, scalapack,
          testSeq, testDist, wantU, wantV, useQR, penalizeDerivative,
          divideCutoff, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
