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
( const Matrix<Real>& d, const Real& rho, const Matrix<Real>& u );

void TestLAPACK
( const Matrix<float>& d, const float& rho, const Matrix<float>& u )
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
        ( &nBLAS, &ip1, d.LockedBuffer(), u.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &sigmaLAPACK,
          dPlusShift.Buffer(), &infoLAPACK );
        wLAPACK(i) = sigmaLAPACK;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK: ",lapackTime," seconds");
    Output("");
}

void TestLAPACK
( const Matrix<double>& d, const double& rho, const Matrix<double>& u )
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
        ( &nBLAS, &ip1, d.LockedBuffer(), u.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &sigmaLAPACK,
          dPlusShift.Buffer(), &infoLAPACK );
        wLAPACK(i) = sigmaLAPACK;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK: ",lapackTime," seconds");
    Output("");
}

template<typename Real>
void GenerateData
( Int n, Matrix<Real>& d, Real& rho, Matrix<Real>& u, bool print )
{
    Uniform( d, n, 1, Real(2), Real(2) );
    Sort( d );
    Gaussian( u, n, 1 );
    u *= Real(1) / FrobeniusNorm( u );
    rho = SampleUniform( Real(1), Real(1)/Real(2) );
    if( print )
    {
        Print( d, "d" );
        Output( "rho=", rho );
        Print( u, "u" );
    }
}

template<typename Real>
void TestSecularHelper
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull )
{
    Output("Testing with ",TypeName<Real>());
    const Int n = d.Height();

    Timer timer;

    SecularSingularValueCtrl<Real> ctrl;
    ctrl.maxIterations = maxIter;
    ctrl.maxCubicIterations = maxCubicIter;
    ctrl.negativeFix = negativeFix;
    ctrl.progress = progress;

    Matrix<Real> wSecular(n,1);
    Int measMinIter=1e9, measMaxIter=0, measTotalIter=0,
        measMinCubicIter=1e9, measMaxCubicIter=0, measTotalCubicIter=0,
        measMinCubicFails=1e9, measMaxCubicFails=0, measTotalCubicFails=0;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        auto info = SecularSingularValue( i, d, rho, u, ctrl );
        wSecular(i) = info.singularValue*info.singularValue;

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

    if( testFull )
    {
        Matrix<Real> A, w;
        Matrix<Real> dSquared;
        Hadamard( d, d, dSquared );
        Diagonal( A, dSquared );
        Syrk( LOWER, NORMAL, rho, u, Real(1), A );
        timer.Start();
        auto hermEigInfo = HermitianEig( LOWER, A, w );
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
  bool lapack )
{
    Matrix<Real> d, u;
    Real rho;
    GenerateData( n, d, rho, u, print );

    TestSecularHelper<Real>
    ( d, rho, u, maxIter, maxCubicIter, negativeFix, progress, print,
      testFull );
    if( lapack )
    {
        TestLAPACK( d, rho, u );
    }
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void TestSecular
( Int n, 
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull )
{
    Matrix<Real> d, u;
    Real rho;
    GenerateData( n, d, rho, u, print );

    TestSecularHelper<Real>
    ( d, rho, u, maxIter, maxCubicIter, negativeFix, progress, print,
      testFull );
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
          lapack );
        TestSecular<double>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          lapack );

#ifdef EL_HAVE_QD
        TestSecular<DoubleDouble>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
        TestSecular<QuadDouble>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
#endif
#ifdef EL_HAVE_QUAD
        TestSecular<Quad>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
#endif
#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        TestSecular<BigFloat>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
#endif

    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
