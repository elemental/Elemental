/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

extern "C" {

void EL_LAPACK(slasd4)
( const BlasInt* n,
  const BlasInt* i,
  const float* d,
  const float* z,
  float* dMinusShift,
  const float* rho,
  float* sigma,
  float* dPlusShift,
  BlasInt* info );

void EL_LAPACK(dlasd4)
( const BlasInt* n,
  const BlasInt* i,
  const double* d,
  const double* z,
  double* dMinusShift,
  const double* rho,
  double* sigma,
  double* dPlusShift,
  BlasInt* info );

} // extern "C"

template<typename Real>
void TestLAPACK
( const Matrix<Real>& d, const Real& rho, const Matrix<Real>& z );

void TestLAPACK
( const Matrix<float>& d, const float& rho, const Matrix<float>& z )
{
    typedef float Real;
    const Int n = d.Height();
    Timer timer;
    Matrix<Real> wLAPACK(n,1), dPlusShift(n,1), dMinusShift(n,1);
    BlasInt nBLAS = n;
    BlasInt infoLAPACK;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        Real sigmaLAPACK;
        const BlasInt ip1 = i+1;
        EL_LAPACK(slasd4)
        ( &nBLAS, &ip1, d.LockedBuffer(), z.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &sigmaLAPACK,
          dPlusShift.Buffer(), &infoLAPACK );
        wLAPACK(i) = sigmaLAPACK;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK secular singular value time: ",lapackTime," seconds");
}

void TestLAPACK
( const Matrix<double>& d, const double& rho, const Matrix<double>& z )
{
    typedef double Real;
    const Int n = d.Height();
    Timer timer;
    Matrix<Real> wLAPACK(n,1), dPlusShift(n,1), dMinusShift(n,1);
    BlasInt nBLAS = n;
    BlasInt infoLAPACK;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        Real sigmaLAPACK;
        const BlasInt ip1 = i+1;
        EL_LAPACK(dlasd4)
        ( &nBLAS, &ip1, d.LockedBuffer(), z.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &sigmaLAPACK,
          dPlusShift.Buffer(), &infoLAPACK );
        wLAPACK(i) = sigmaLAPACK;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK secular singular value time: ",lapackTime," seconds");
}

template<typename Real>
void GenerateData
( Int n, Matrix<Real>& d, Real& rho, Matrix<Real>& z, bool print )
{
    // Implicitly form a matrix
    //
    //   M = | sqrt(rho)*z(0), sqrt(rho)*z(1), ..., sqrt(rho)*z(n-1) |
    //       |                      d(1),                            |
    //       |                                 .                     |
    //       |                                                d(n-1) |
    //
    // where 0 = d(0) <= d(1) <= d(2) <= ... <= d(n-1).
    //
    Uniform( d, n, 1, Real(2), Real(2) );
    Sort( d );
    d(0) = 0;
    Gaussian( z, n, 1 );
    z *= Real(1) / FrobeniusNorm( z );
    rho = SampleUniform( Real(1), Real(1)/Real(2) );
    if( print )
    {
        Print( d, "d" );
        Output( "rho=", rho );
        Print( z, "z" );
    }
}

template<typename Real>
void PrintSVDResiduals
( UpperOrLower uplo,
  const Matrix<Real>& mainDiag,
  const Matrix<Real>& offDiag,
  const Matrix<Real>& U,
  const Matrix<Real>& s,
  const Matrix<Real>& V,
  bool print )
{
    DEBUG_CSE
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
    Matrix<Real> A;
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
    Gemm( NORMAL, ADJOINT, Real(-1), UModMin, VModMin, Real(1), A );
    if( print )
        Print( A, "E" );
    const Real residFrob = FrobeniusNorm( A );
    Output("|| A - U Sigma V' ||_F / || A ||_F = ",residFrob/AFrob);
    // TODO(poulson): Failure condition

    // Check the unitarity of U
    Matrix<Real> E; 
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

template<typename Real>
void TestDivideAndConquer
( Int m,
  bool wantU,
  bool wantV,
  Int cutoff,
  Int maxIter,  
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print )
{
    Output("Testing DivideAndConquer(",cutoff,") with ",TypeName<Real>());

    BidiagSVDCtrl<Real> ctrl;
    ctrl.wantU = wantU;
    ctrl.wantV = wantV;
    ctrl.progress = progress;
    ctrl.dcCtrl.exploitStructure = true;
    ctrl.dcCtrl.cutoff = cutoff;
    ctrl.dcCtrl.secularCtrl.maxIterations = maxIter;
    ctrl.dcCtrl.secularCtrl.maxCubicIterations = maxCubicIter;
    ctrl.dcCtrl.secularCtrl.negativeFix = negativeFix;
    ctrl.dcCtrl.secularCtrl.progress = progress;

    const bool square = true;
    const Int n = ( square ? m : m+1 );
    Matrix<Real> mainDiag, superDiag;
    Uniform( mainDiag, m, 1 );
    Uniform( superDiag, n-1, 1 );
    if( print )
    {
        Print( mainDiag, "mainDiag" );
        Print( superDiag, "superDiag" );
    }

    Timer timer;

    Matrix<Real> s;
    Matrix<Real> U, V;
    timer.Start();
    auto info = BidiagSVD( UPPER, mainDiag, superDiag, U, s, V, ctrl );
    const auto& dcInfo = info.dcInfo;
    const auto& secularInfo = dcInfo.secularInfo;
    const auto& deflationInfo = dcInfo.deflationInfo;
    Output("Bidiag D&C: ",timer.Stop()," seconds");
    Output("  num deflations: ",deflationInfo.numDeflations);
    Output("    small diagonal: ",deflationInfo.numSmallDiagonalDeflations);
    Output("    close diagonal: ",deflationInfo.numCloseDiagonalDeflations); 
    Output("    small update:   ",deflationInfo.numSmallUpdateDeflations);
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
        PrintSVDResiduals( UPPER, mainDiag, superDiag, U, s, V, print );
        PopIndent();
    }

    // Compute the residual with the QR algorithm (with relative-to-max tol)
    timer.Start();
    ctrl.useQR = true;
    BidiagSVD( UPPER, mainDiag, superDiag, U, s, V, ctrl );
    Output("QR algorithm: ",timer.Stop()," seconds");
    if( wantU && wantV )
    {
        Output("Residuals with QR:"); 
        PushIndent();
        PrintSVDResiduals( UPPER, mainDiag, superDiag, U, s, V, print );
        PopIndent();
    }
}

template<typename Real>
void TestSecularHelper
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull,
  bool wantU,
  bool wantV,
  Int divideCutoff )
{
    const Int n = d.Height();

    Timer timer;

    SecularSingularValueCtrl<Real> ctrl;
    ctrl.maxIterations = maxIter;
    ctrl.maxCubicIterations = maxCubicIter;
    ctrl.negativeFix = negativeFix;
    ctrl.progress = progress;

    Matrix<Real> s(n,1), wSecular(n,1);
    Int measMinIter=1e9, measMaxIter=0, measTotalIter=0,
        measMinCubicIter=1e9, measMaxCubicIter=0, measTotalCubicIter=0,
        measMinCubicFails=1e9, measMaxCubicFails=0, measTotalCubicFails=0;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        auto info = SecularSingularValue( i, d, rho, z, s(i), ctrl );
        wSecular(i) = s(i) * s(i);

        measMinIter = Min( measMinIter, info.numIterations );
        measMaxIter = Max( measMaxIter, info.numIterations );
        measTotalIter += info.numIterations;

        measMinCubicIter = Min( measMinCubicIter, info.numCubicIterations );
        measMaxCubicIter = Max( measMaxCubicIter, info.numCubicIterations );
        measTotalCubicIter += info.numCubicIterations;

        measMinCubicFails = Min( measMinCubicFails, info.numCubicFailures );
        measMaxCubicFails = Max( measMaxCubicFails, info.numCubicFailures );
        measTotalCubicFails += info.numCubicFailures;
    }
    const Real secularTime = timer.Stop();
    Output("Secular: ",secularTime," seconds");
    Output
    ("Iterations [min/max/total]: ",
     measMinIter,"/",measMaxIter,"/",measTotalIter);
    Output
    ("Cubic iter's [min/max/total]: ",
     measMinCubicIter,"/",measMaxCubicIter,"/",measTotalCubicIter);
    Output
    ("Cubic failures [min/max/total]: ",
     measMinCubicFails,"/",measMaxCubicFails,"/",measTotalCubicFails);
    Output("");

    // Now compute the singular values and vectors. We recompute the singular
    // values to avoid interfering with the timing experiment above.
    Matrix<Real> U, V;
    timer.Start();
    SecularSVD( d, rho, z, U, s, V, ctrl );
    const double secularSVDTime = timer.Stop();
    Output("Secular SVD: ",secularSVDTime," seconds");
    if( print )
    {
        Print( U, "U" );
        Print( s, "s" );
        Print( V, "V" );
    }

    // Explicitly form the matrix M
    Matrix<Real> M;
    Zeros( M, n, n );
    for( Int j=0; j<n; ++j )
        M(0,j) = z(j)*Sqrt(rho);
    for( Int j=1; j<n; ++j )
        M(j,j) = d(j);
    const Real MFrob = FrobeniusNorm( M );
    Output("|| M ||_F = ",MFrob);
    if( print )
        Print( M, "M" );

    // Test the Singular Value Decomposition of M
    Matrix<Real> UScaled( U );
    DiagonalScale( RIGHT, NORMAL, s, UScaled );
    Matrix<Real> E( M );
    Gemm( NORMAL, ADJOINT, Real(-1), UScaled, V, Real(1), E );
    const Real EFrob = FrobeniusNorm( E );
    Output("|| M - U Sigma V' ||_F = ",EFrob);

    // Test the orthonormality of U and V
    Identity( E, n, n );
    Gemm( NORMAL, ADJOINT, Real(-1), U, U, Real(1), E );
    const Real UOrthError = FrobeniusNorm( E );
    Output("|| I - U U' ||_F = ",UOrthError);
    Identity( E, n, n );
    Gemm( NORMAL, ADJOINT, Real(-1), V, V, Real(1), E );
    const Real VOrthError = FrobeniusNorm( E );
    Output("|| I - V V' ||_F = ",VOrthError);

    if( testFull )
    {
        Matrix<Real> A, w;
        Matrix<Real> dSquared;
        Hadamard( d, d, dSquared );
        Diagonal( A, dSquared );
        Syrk( LOWER, NORMAL, rho, z, Real(1), A );
        timer.Start();
        HermitianEig( LOWER, A, w );
        const Real fullTime = timer.Stop();
        Output("Full Hermitian: ",fullTime," seconds");
        if( print )
            Print( w, "w" );

        auto wDiff( w );
        wDiff -= wSecular;
        const Real diffNorm = FrobeniusNorm( wDiff );
        Output("|| w - wSecular ||_F = ", diffNorm);
        Output("");
    }

    TestDivideAndConquer<Real>
    ( n, wantU, wantV, divideCutoff, maxIter, maxCubicIter, negativeFix,
      progress, print );
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void TestSecular
( Int n, 
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull,
  bool wantU,
  bool wantV,
  Int divideCutoff,
  bool lapack )
{
    Output("Testing with ",TypeName<Real>());
    Matrix<Real> d, z;
    Real rho;
    GenerateData( n, d, rho, z, print );

    if( lapack )
    {
        TestLAPACK( d, rho, z );
    }
    TestSecularHelper<Real>
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, progress, print,
      testFull, wantU, wantV, divideCutoff );
    Output("");
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void TestSecular
( Int n, 
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull,
  bool wantU,
  bool wantV,
  Int divideCutoff )
{
    Output("Testing with ",TypeName<Real>());
    Matrix<Real> d, z;
    Real rho;
    GenerateData( n, d, rho, z, print );

    TestSecularHelper<Real>
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, progress, print,
      testFull, wantU, wantV, divideCutoff );
    Output("");
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","matrix size",100);
        const Int maxIter = Input("--maxIter","max iterations",400);
        const Int maxCubicIter = Input("--maxCubicIter","max cubic iter's",40);
        const Int flipOrClipInt = Input("--flipOrClip","0: flip, 1: clip",1);
        const Int divideCutoff = Input("--divideCutoff","D&C cutoff",60);
        const bool wantU = Input("--wantU","compute U?",true);
        const bool wantV = Input("--wantV","compute V?",true);
        const bool progress = Input("--progress","print progress?",false);
        const bool testFull = Input("--testFull","test full eigensolver?",true);
        const bool lapack =
          Input("--lapack","test against LAPACK's secular solver?",true);
        const bool print = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();

        FlipOrClip negativeFix = static_cast<FlipOrClip>(flipOrClipInt);

        TestSecular<float>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          wantU, wantV, divideCutoff, lapack );
        TestSecular<double>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          wantU, wantV, divideCutoff, lapack );

#ifdef EL_HAVE_QD
        TestSecular<DoubleDouble>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          wantU, wantV, divideCutoff );
        TestSecular<QuadDouble>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          wantU, wantV, divideCutoff );
#endif
#ifdef EL_HAVE_QUAD
        TestSecular<Quad>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          wantU, wantV, divideCutoff );
#endif
#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        TestSecular<BigFloat>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          wantU, wantV, divideCutoff );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
