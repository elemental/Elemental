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
    Output("LAPACK: ",lapackTime," seconds");
    Output("");
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
    Output("LAPACK: ",lapackTime," seconds");
    Output("");
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
void TestSecularHelper
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
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

    Matrix<Real> s(n,1), wSecular(n,1);
    Int measMinIter=1e9, measMaxIter=0, measTotalIter=0,
        measMinCubicIter=1e9, measMaxCubicIter=0, measTotalCubicIter=0,
        measMinCubicFails=1e9, measMaxCubicFails=0, measTotalCubicFails=0;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        auto info = SecularSingularValue( i, d, rho, z, ctrl );
        s(i) = info.singularValue;
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

    // Now compute the singular values and vectors. We recompute the singular
    // values to avoid interfering with the timing experiment above.
    Matrix<Real> DMinusShift(n,n), DPlusShift(n,n);
    timer.Start();
    for( Int j=0; j<n; ++j )
    {
        auto dMinusShift = DMinusShift(ALL,IR(j));
        auto dPlusShift = DPlusShift(ALL,IR(j));
        auto info =
          SecularSingularValue( j, d, rho, z, dMinusShift, dPlusShift, ctrl );
        s(j) = info.singularValue;
    }
    const double secularStore = timer.Stop();
    Output("Secular solver with shift storing: ",secularStore," seconds");

    // Compute a vector z which would produce the given singular values to
    // high relative accuracy. Keep in mind that the following absorbs the
    // sqrt(rho) factor into zCorrect.
    Matrix<Real> zCorrect(n,1);
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        // See Eq. (3.6) from Gu and Eisenstat's Technical Report
        // "A Divide-and-Conquer Algorithm for the Bidiagonal SVD"
        // [CITATION].
        Real prod = DPlusShift(i,n-1)*DMinusShift(i,n-1);
        for( Int k=0; k<i; ++k )
        {
            const Real dSqDiff = (d(k)+d(i))*(d(k)-d(i));
            prod *= (DPlusShift(i,k)*DMinusShift(i,k)) / dSqDiff;
        }
        for( Int k=i; k<n-1; ++k )
        {
            const Real dSqDiff = (d(k+1)+d(i))*(d(k+1)-d(i));
            prod *= (DPlusShift(i,k)*DMinusShift(i,k)) / dSqDiff;
        }
        zCorrect(i) = Sgn(z(i),false)*Sqrt(Abs(prod));
    }
    const double correctedVecTime = timer.Stop();
    Output("Corrected vector formation time: ",correctedVecTime," seconds");
    if( print )
    {
        Print( zCorrect, "zCorrect" );
        auto zScaled( z );
        zScaled *= Sqrt(rho);
        Print( zScaled, "sqrt(rho) z" );
    }

    Matrix<Real> U(n,n), V(n,n);
    timer.Start();
    for( Int j=0; j<n; ++j )
    {
        // Compute the j'th left singular vectors via Eq. (3.4).
        auto u = U(ALL,IR(j));
        u(0) = -1;
        for( Int i=1; i<n; ++i )
        {
            u(i) = (d(i)*zCorrect(i)) / (DPlusShift(i,j)*DMinusShift(i,j));
        }
        u *= Real(1) / FrobeniusNorm( u );

        // Compute the j'th right singular vector via Eq. (3.3)
        auto v = V(ALL,IR(j));
        for( Int i=0; i<n; ++i )
        {
            v(i) = zCorrect(i) / (DPlusShift(i,j)*DMinusShift(i,j));
        }
        v *= Real(1) / FrobeniusNorm( v );
    }
    const double secularVecTime = timer.Stop();
    Output("Singular vector formation: ",secularVecTime," seconds");
    if( print )
    {
        Print( U, "U" );
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
    Matrix<Real> d, z;
    Real rho;
    GenerateData( n, d, rho, z, print );

    TestSecularHelper<Real>
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, progress, print,
      testFull );
    if( lapack )
    {
        TestLAPACK( d, rho, z );
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
    Matrix<Real> d, z;
    Real rho;
    GenerateData( n, d, rho, z, print );

    TestSecularHelper<Real>
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, progress, print,
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
