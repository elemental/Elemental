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

void EL_LAPACK(slaed4)
( const BlasInt* n,
  const BlasInt* i,
  const float* d,
  const float* z,
  float* dMinusShift,
  const float* rho,
  float* sigma,
  BlasInt* info );

void EL_LAPACK(dlaed4)
( const BlasInt* n,
  const BlasInt* i,
  const double* d,
  const double* z,
  double* dMinusShift,
  const double* rho,
  double* sigma,
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
    Matrix<Real> w(n,1), dMinusShift(n,1);
    BlasInt nBLAS = n;
    BlasInt info;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        Real lambda;
        const BlasInt ip1 = i+1;
        EL_LAPACK(slaed4)
        ( &nBLAS, &ip1, d.LockedBuffer(), z.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &lambda, &info );
        if( info != 0 )
            RuntimeError("LAPACK's slaed4 did not converge");
        w(i) = lambda;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK secular eigenvalue time: ",lapackTime," seconds");
}

void TestLAPACK
( const Matrix<double>& d, const double& rho, const Matrix<double>& z )
{
    typedef double Real;
    const Int n = d.Height();
    Timer timer;
    Matrix<Real> w(n,1), dMinusShift(n,1);
    BlasInt nBLAS = n;
    BlasInt info;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        Real lambda;
        const BlasInt ip1 = i+1;
        EL_LAPACK(dlaed4)
        ( &nBLAS, &ip1, d.LockedBuffer(), z.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &lambda, &info );
        if( info != 0 )
            RuntimeError("LAPACK's dlaed4 did not converge");
        w(i) = lambda;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK secular eigenvalue time: ",lapackTime," seconds");
}

template<typename Real>
void GenerateData
( Int n, Matrix<Real>& d, Real& rho, Matrix<Real>& z, bool print )
{
    // Implicitly form a matrix d + rho z z^T, with
    // d(0) <= d(1) <= d(2) <= ... <= d(n-1) and || z ||_2 = 1.
    //
    Uniform( d, n, 1, Real(0), Real(2) );
    Sort( d );
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
void TestSecularHelper
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool penalizeDerivative,
  bool progress,
  bool print,
  bool testFull )
{
    const Int n = d.Height();

    Timer timer;

    SecularEVDCtrl<Real> ctrl;
    ctrl.maxIterations = maxIter;
    ctrl.negativeFix = negativeFix;
    ctrl.penalizeDerivative = penalizeDerivative;
    ctrl.progress = progress;
    ctrl.cubicCtrl.maxIterations = maxCubicIter;
    ctrl.cubicCtrl.negativeFix = negativeFix;

    Matrix<Real> w(n,1);
    Int measMinIter=1e9, measMaxIter=0, measTotalIter=0,
        measMinCubicIter=1e9, measMaxCubicIter=0, measTotalCubicIter=0,
        measMinCubicFails=1e9, measMaxCubicFails=0, measTotalCubicFails=0;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        auto info = SecularEigenvalue( i, d, rho, z, w(i), ctrl );

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

    // Now compute the eigenvalues and vectors. We recompute the eigenvalues
    // to avoid interfering with the timing experiment above.
    Matrix<Real> Q;
    timer.Start();
    SecularEVD( d, rho, z, w, Q, ctrl );
    const double secularEVDTime = timer.Stop();
    Output("Secular EVD: ",secularEVDTime," seconds");
    if( print )
    {
        Print( w, "w" );
        Print( Q, "Q" );
    }

    // Explicitly form the matrix M = diag(d) + rho z z^T
    Matrix<Real> M;
    Diagonal( M, d );
    Geru( rho, z, z, M );
    const Real MFrob = FrobeniusNorm( M );
    Output("|| M ||_F = ",MFrob);
    if( print )
        Print( M, "M" );

    // Test the EigenValue Decomposition of M
    Matrix<Real> QScaled( Q );
    DiagonalScale( RIGHT, NORMAL, w, QScaled );
    Matrix<Real> E( M );
    Gemm( NORMAL, ADJOINT, Real(-1), QScaled, Q, Real(1), E );
    const Real EFrob = FrobeniusNorm( E );
    Output("|| M - Q diag(w) Q' ||_F = ",EFrob);

    // Test the orthonormality of Q
    Identity( E, n, n );
    Gemm( NORMAL, ADJOINT, Real(-1), Q, Q, Real(1), E );
    const Real QOrthError = FrobeniusNorm( E );
    Output("|| I - Q Q' ||_F = ",QOrthError);

    if( testFull )
    {
        Matrix<Real> A, wFull, QFull;
        Diagonal( A, d );
        Syrk( LOWER, NORMAL, rho, z, Real(1), A );
        timer.Start();
        HermitianEig( LOWER, A, wFull, QFull );
        const Real fullTime = timer.Stop();
        Output("Full Hermitian: ",fullTime," seconds");
        if( print )
            Print( wFull, "wFull" );

        auto wDiff( wFull );
        wDiff -= w;
        const Real diffNorm = FrobeniusNorm( wDiff );
        Output("|| wFull - w ||_F = ", diffNorm);
        Output("");
    }
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void TestSecular
( Int n, 
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool penalizeDerivative,
  bool progress,
  bool print,
  bool testFull,
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
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, penalizeDerivative,
      progress, print, testFull );
    Output("");
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void TestSecular
( Int n, 
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool penalizeDerivative,
  bool progress,
  bool print,
  bool testFull )
{
    Output("Testing with ",TypeName<Real>());
    Matrix<Real> d, z;
    Real rho;
    GenerateData( n, d, rho, z, print );

    TestSecularHelper<Real>
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, penalizeDerivative,
      progress, print, testFull );
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
        const bool clip = Input("--clip","clip negative?",true);
        const bool penalizeDerivative =
          Input("--penalizeDerivative","penalize derivative?",true);
        const bool progress = Input("--progress","print progress?",false);
        const bool testFull = Input("--testFull","test full eigensolver?",true);
        const bool lapack =
          Input("--lapack","test against LAPACK's secular solver?",true);
        const bool print = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();

        const FlipOrClip negativeFix = 
          ( clip ? CLIP_NEGATIVES : FLIP_NEGATIVES );

        TestSecular<float>
        ( n, maxIter, maxCubicIter, negativeFix, penalizeDerivative, progress,
          print, testFull, lapack );
        TestSecular<double>
        ( n, maxIter, maxCubicIter, negativeFix, penalizeDerivative, progress,
          print, testFull, lapack );

#ifdef EL_HAVE_QD
        TestSecular<DoubleDouble>
        ( n, maxIter, maxCubicIter, negativeFix, penalizeDerivative, progress,
          print, testFull );
        TestSecular<QuadDouble>
        ( n, maxIter, maxCubicIter, negativeFix, penalizeDerivative, progress,
          print, testFull );
#endif
#ifdef EL_HAVE_QUAD
        TestSecular<Quad>
        ( n, maxIter, maxCubicIter, negativeFix, penalizeDerivative, progress,
          print, testFull );
#endif
#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        TestSecular<BigFloat>
        ( n, maxIter, maxCubicIter, negativeFix, penalizeDerivative, progress,
          print, testFull );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
