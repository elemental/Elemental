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
void PrintSVDResiduals
( UpperOrLower uplo,
  const Matrix<F>& mainDiag,
  const Matrix<F>& offDiag,
  const Matrix<F>& U,
  const Matrix<Base<F>>& s,
  const Matrix<F>& V,
  bool print )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = U.Height();
    const Int n = V.Height();
    const Int minDim = Min( m, n );
    Output("m=",m,", n=",n,", minDim=",minDim);
    if( print )
    {
        Print( U, "U" ); 
        Print( s, "s" );
        Print( V, "V" );
    }

    // Explicitly form A
    Matrix<F> A;
    Zeros( A, m, n );
    SetDiagonal( A, mainDiag, 0 );
    if( uplo == UPPER )
        SetDiagonal( A, offDiag, 1 );
    else
        SetDiagonal( A, offDiag, -1 );
    if( print )
        Print( A, "A" );
    const Real AFrob = FrobeniusNorm( A );
    Output("|| A ||_F = ",AFrob);

    // Check || A - U Sigma V^T ||_F
    // TODO(poulson): Introduce diagonally-scaled general outer product
    auto UMod( U );
    auto VMod( V );
    auto UModMin = UMod( ALL, IR(0,minDim) );
    auto VModMin = VMod( ALL, IR(0,minDim) );
    DiagonalScale( RIGHT, NORMAL, s(IR(0,minDim),ALL), UModMin );
    Gemm( NORMAL, ADJOINT, F(-1), UModMin, VModMin, F(1), A );
    if( print )
        Print( A, "E" );
    const Real residFrob = FrobeniusNorm( A );
    Output("|| A - U Sigma V' ||_F / || A ||_F = ",residFrob/AFrob);
    // TODO(poulson): Provide a rigorous motivation for this bound
    if( (residFrob/AFrob) > Real(50) )
        LogicError("SVD residual was unacceptably large");

    // Check the unitarity of U
    Matrix<F> E; 
    Identity( E, U.Width(), U.Width() );
    Herk( LOWER, ADJOINT, Real(-1), U, Real(1), E );
    const Real UOrthogFrob = HermitianFrobeniusNorm( LOWER, E );
    Output("|| I - U' U ||_F = ",UOrthogFrob);
    // TODO(poulson): Failure condition
    
    // Check the unitarity of V
    Identity( E, V.Width(), V.Width() );
    Herk( LOWER, ADJOINT, Real(-1), V, Real(1), E );
    const Real VOrthogFrob = HermitianFrobeniusNorm( LOWER, E );
    Output("|| I - V' V ||_F = ",UOrthogFrob);
    // TODO(poulson): Failure condition
}

template<typename F>
void TestDivideAndConquer
( Int m,
  UpperOrLower uplo, 
  bool wantU,
  bool wantV,
  Int cutoff,
  Int maxIter,  
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print )
{
    Output("Testing DivideAndConquer(",cutoff,") with ",TypeName<F>());
    typedef Base<F> Real;

    // TODO(poulson): Make this configurable
    const bool penalizeDerivative = false;

    BidiagSVDCtrl<Real> ctrl;
    ctrl.wantU = wantU;
    ctrl.wantV = wantV;
    ctrl.progress = progress;
    ctrl.dcCtrl.exploitStructure = true;
    ctrl.dcCtrl.cutoff = cutoff;
    ctrl.dcCtrl.secularCtrl.maxIterations = maxIter;
    ctrl.dcCtrl.secularCtrl.negativeFix = negativeFix;
    ctrl.dcCtrl.secularCtrl.penalizeDerivative = penalizeDerivative;
    ctrl.dcCtrl.secularCtrl.progress = progress;
    ctrl.dcCtrl.secularCtrl.cubicCtrl.maxIterations = maxCubicIter;
    ctrl.dcCtrl.secularCtrl.cubicCtrl.negativeFix = negativeFix;

    const bool square = true;
    const Int n = ( square ? m : m+1 );
    const Int mainDiagHeight = ( uplo==UPPER ? m : n );
    const Int offDiagHeight = ( uplo==UPPER ? n-1 : m-1 );
    Matrix<F> mainDiag, offDiag;
    Uniform( mainDiag, mainDiagHeight, 1 );
    Uniform( offDiag, offDiagHeight, 1 );
    if( print )
    {
        Print( mainDiag, "mainDiag" );
        Print( offDiag, "offDiag" );
    }

    Timer timer;

    Matrix<Real> s;
    Matrix<F> U, V;
    timer.Start();
    auto info = BidiagSVD( uplo, mainDiag, offDiag, U, s, V, ctrl );
    const auto& dcInfo = info.dcInfo;
    const auto& secularInfo = dcInfo.secularInfo;
    Output("Bidiag D&C: ",timer.Stop()," seconds");
    Output("  num deflations: ",secularInfo.numDeflations);
    Output("    small diagonal: ",secularInfo.numSmallDiagonalDeflations);
    Output("    close diagonal: ",secularInfo.numCloseDiagonalDeflations); 
    Output("    small update:   ",secularInfo.numSmallUpdateDeflations);
    Output("  num secular iterations: ",secularInfo.numIterations);
    Output("  num secular alternations: ",secularInfo.numAlternations);
    Output("  num secular cubic iter's: ",secularInfo.numCubicIterations);
    Output("  num secular cubic failures: ",secularInfo.numCubicFailures);
    if( print )
    {
        Print( U, "U" );
        Print( s, "s" );
        Print( V, "V" );
    }

    if( wantU && wantV )
    {
        // Compute the residuals
        Output("Residuals after D&C:");
        PushIndent();
        PrintSVDResiduals( uplo, mainDiag, offDiag, U, s, V, print );
        PopIndent();
    }

    // Compute the residual with the QR algorithm (with relative-to-max tol)
    timer.Start();
    ctrl.useQR = true;
    BidiagSVD( uplo, mainDiag, offDiag, U, s, V, ctrl );
    Output("QR algorithm: ",timer.Stop()," seconds");
    if( wantU && wantV )
    {
        Output("Residuals with QR:"); 
        PushIndent();
        PrintSVDResiduals( uplo, mainDiag, offDiag, U, s, V, print );
        PopIndent();
    }
    Output("");
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","matrix size",100);
        const bool upper = Input("--upper","upper bidiagonal?",true);
        const Int maxIter = Input("--maxIter","max iterations",400);
        const Int maxCubicIter = Input("--maxCubicIter","max cubic iter's",40);
        const bool clip = Input("--clip","clip negative?",true);
        const Int divideCutoff = Input("--divideCutoff","D&C cutoff",60);
        const bool wantU = Input("--wantU","compute U?",true);
        const bool wantV = Input("--wantV","compute V?",true);
        const bool progress = Input("--progress","print progress?",false);
        const bool print = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();

        const UpperOrLower uplo = ( upper ? UPPER : LOWER );
        const FlipOrClip negativeFix =
          ( clip ? CLIP_NEGATIVES : FLIP_NEGATIVES );

        TestDivideAndConquer<float>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
        TestDivideAndConquer<Complex<float>>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );

        TestDivideAndConquer<double>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
        TestDivideAndConquer<Complex<double>>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );

#ifdef EL_HAVE_QD
        TestDivideAndConquer<DoubleDouble>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
        TestDivideAndConquer<Complex<DoubleDouble>>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );

        TestDivideAndConquer<QuadDouble>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
        TestDivideAndConquer<Complex<QuadDouble>>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
#endif
#ifdef EL_HAVE_QUAD
        TestDivideAndConquer<Quad>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
        TestDivideAndConquer<Complex<Quad>>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
#endif
#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        TestDivideAndConquer<BigFloat>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
        TestDivideAndConquer<Complex<BigFloat>>
        ( n, uplo, wantU, wantV, divideCutoff, maxIter, maxCubicIter,
          negativeFix, progress, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
