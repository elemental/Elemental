/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

using El::BlasInt;
using El::scomplex;
using El::dcomplex;

extern "C" {

// Machine constants
float EL_LAPACK(slamch)( const char* cmach );
double EL_LAPACK(dlamch)( const char* cmach );

// Safe norms
float  EL_LAPACK(slapy2)( const float * alpha, const float * beta );
double EL_LAPACK(dlapy2)( const double* alpha, const double* beta );
float EL_LAPACK(slapy3)
( const float * alpha, const float * beta, const float * gamma );
double EL_LAPACK(dlapy3)
( const double* alpha, const double* beta, const double* gamma );

// Copy matrices
void EL_LAPACK(slacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const float* A, const BlasInt* lda, float* B, const BlasInt* ldb );
void EL_LAPACK(dlacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const double* A, const BlasInt* lda, double* B, const BlasInt* ldb );
void EL_LAPACK(clacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const scomplex* A, const BlasInt* lda, scomplex* B, const BlasInt* ldb );
void EL_LAPACK(zlacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const dcomplex* A, const BlasInt* lda, dcomplex* B, const BlasInt* ldb );

// Safely compute a Givens rotation
void EL_LAPACK(slartg)
( const float* phi, const float* gamma,
  float* c, float* s, float* rho );
void EL_LAPACK(dlartg)
( const double* phi, const double* gamma,
  double* c, double* s, double* rho );
void EL_LAPACK(clartg)
( const scomplex* phi, const scomplex* gamma,
  float* c, scomplex* s, scomplex* rho );
void EL_LAPACK(zlartg)
( const dcomplex* phi, const dcomplex* gamma,
  double* c, dcomplex* s, dcomplex* rho );

// Symmetric tridiagonal eigensolvers (via MRRR)
void EL_LAPACK(sstevr)
( const char* job, const char* range, const BlasInt* n,
  float* d, float* e, const float* vl, const float* vu, 
  const BlasInt* il, const BlasInt* iu, const float* absTol, BlasInt * m, 
  float* w, float* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  float* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );
void EL_LAPACK(dstevr)
( const char* job, const char* range, const BlasInt* n,
  double* d, double* e, const double* vl, const double* vu, 
  const BlasInt* il, const BlasInt* iu, const double* absTol, BlasInt * m, 
  double* w, double* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  double* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );

// Hermitian eigensolvers (via MRRR)
void EL_LAPACK(ssyevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  float* A, const BlasInt* ldA, const float* vl, const float* vu, 
  const BlasInt* il, const BlasInt* iu, const float* absTol, BlasInt * m, 
  float* w, float* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  float* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );
void EL_LAPACK(dsyevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  double* A, const BlasInt* ldA, const double* vl, const double* vu, 
  const BlasInt* il, const BlasInt* iu, const double* absTol, BlasInt * m, 
  double* w, double* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  double* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );
void EL_LAPACK(cheevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  scomplex* A, const BlasInt* ldA, const float* vl, const float* vu, 
  const BlasInt* il, const BlasInt* iu, const float* absTol, BlasInt* m,
  float* w, scomplex* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  scomplex* work, const BlasInt* workSize, 
  float* rWork, const BlasInt* rWorkSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, BlasInt* info );
void EL_LAPACK(zheevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA, const double* vl, const double* vu, 
  const BlasInt* il, const BlasInt* iu, const double* absTol, BlasInt* m,
  double* w, dcomplex* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  dcomplex* work, const BlasInt* workSize, 
  double* rWork, const BlasInt* rWorkSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, BlasInt* info );

// Bidiagonal DQDS
void EL_LAPACK(slasq1)
( const BlasInt* n, float* d, float* e, float* work, BlasInt* info );
void EL_LAPACK(dlasq1)
( const BlasInt* n, double* d, double* e, double* work, BlasInt* info );

// Bidiagonal QR
void EL_LAPACK(sbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVT, const BlasInt* numRowsU,
  const BlasInt* numColsC, float* d, float* e, 
  float* VTrans, const BlasInt* ldVT,
  float* U, const BlasInt* ldU, float* C, const BlasInt* ldC, 
  float* work, BlasInt* info );
void EL_LAPACK(dbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVT, const BlasInt* numRowsU,
  const BlasInt* numColsC, double* d, double* e,
  double* VTrans, const BlasInt* ldVT, double* U, const BlasInt* ldU,
  double* C, const BlasInt* ldC, double* work, BlasInt* info );
void EL_LAPACK(cbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVH, const BlasInt* numRowsU,
  const BlasInt* numColsC, float* d, float* e,
  scomplex* VH, const BlasInt* ldVH, scomplex* U, const BlasInt* ldU,
  scomplex* C, const BlasInt* ldC, float* work, BlasInt* info );
void EL_LAPACK(zbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVH, const BlasInt* numRowsU,
  const BlasInt* numColsC, double* d, double* e,
  dcomplex* VH, const BlasInt* ldVH, dcomplex* U, const BlasInt* ldU,
  dcomplex* C, const BlasInt* ldC, double* work, BlasInt* info );

// Divide and Conquer SVD
void EL_LAPACK(sgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n, 
  float* A, const BlasInt* ldA,
  float* s, float* U, const BlasInt* ldu, 
  float* VTrans, const BlasInt* ldvt,
  float* work, const BlasInt* workSize, BlasInt* iWork, BlasInt* info );
void EL_LAPACK(dgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n, 
  double* A, const BlasInt* ldA,
  double* s, double* U, const BlasInt* ldu, 
  double* VTrans, const BlasInt* ldvt,
  double* work, const BlasInt* workSize, BlasInt* iWork, BlasInt* info );
void EL_LAPACK(cgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n,
  scomplex* A, const BlasInt* ldA, float* s,
  scomplex* U, const BlasInt* ldu, scomplex* VTrans, const BlasInt* ldvt,
  scomplex* work, const BlasInt* workSize, float* rWork,
  BlasInt* iWork, BlasInt* info );
void EL_LAPACK(zgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA, double* s,
  dcomplex* U, const BlasInt* ldu, dcomplex* VH, const BlasInt* ldva,
  dcomplex* work, const BlasInt* workSize, double* rWork,
  BlasInt* iWork, BlasInt* info );

// QR-algorithm SVD [DQDS when no singular vectors desired]
void EL_LAPACK(sgesvd)
( const char* jobU, const char* jobVT, const BlasInt* m, const BlasInt* n,
  float* A, const BlasInt* ldA,
  float* s, float* U, const BlasInt* ldu, float* VTrans, const BlasInt* ldvt,
  float* work, const BlasInt* workSize, BlasInt* info );
void EL_LAPACK(dgesvd)
( const char* jobU, const char* jobVT, const BlasInt* m, const BlasInt* n,
  double* A, const BlasInt* ldA,
  double* s, double* U, const BlasInt* ldu, double* VTrans, const BlasInt* ldvt,
  double* work, const BlasInt* workSize, BlasInt* info );
void EL_LAPACK(cgesvd)
( const char* jobU, const char* jobVH, const BlasInt* m, const BlasInt* n,
  scomplex* A, const BlasInt* ldA, float* s,
  scomplex* U, const BlasInt* ldu, scomplex* VTrans, const BlasInt* ldvt,
  scomplex* work, const BlasInt* workSize, float* rWork, BlasInt* info );
void EL_LAPACK(zgesvd)
( const char* jobU, const char* jobVH, const BlasInt* m, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA, double* s,
  dcomplex* U, const BlasInt* ldu, dcomplex* VH, const BlasInt* ldva,
  dcomplex* work, const BlasInt* workSize, double* rWork, BlasInt* info );

// Reduction to Hessenberg form
void EL_LAPACK(sgehrd)
( const BlasInt* n, const BlasInt* ilo, const BlasInt* ihi, 
  float* A, const BlasInt* ldA,
  float* tau, float* work, const BlasInt* workSize, BlasInt* info );
void EL_LAPACK(dgehrd)
( const BlasInt* n, const BlasInt* ilo, const BlasInt* ihi, 
  double* A, const BlasInt* ldA,
  double* tau, double* work, const BlasInt* workSize, BlasInt* info );
void EL_LAPACK(cgehrd)
( const BlasInt* n, const BlasInt* ilo, const BlasInt* ihi, 
  scomplex* A, const BlasInt* ldA,
  scomplex* tau, scomplex* work, const BlasInt* workSize, BlasInt* info );
void EL_LAPACK(zgehrd)
( const BlasInt* n, const BlasInt* ilo, const BlasInt* ihi, 
  dcomplex* A, const BlasInt* ldA,
  dcomplex* tau, dcomplex* work, const BlasInt* workSize, BlasInt* info );

// Generates a unitary matrix defined as the product of Householder reflectors
void EL_LAPACK(sorghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  float* A, const BlasInt* ldA,
  const float* tau,
  float* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(dorghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  double* A, const BlasInt* ldA,
  const double* tau,
  double* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(cunghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  scomplex* A, const BlasInt* ldA,
  const scomplex* tau,
  scomplex* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(zunghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  dcomplex* A, const BlasInt* ldA,
  const dcomplex* tau,
  dcomplex* work, const BlasInt* workSize,
  BlasInt* info );

// Hessenberg QR algorithm
void EL_LAPACK(shseqr)
( const char* job, const char* compZ, const BlasInt* n, 
  const BlasInt* ilo, const BlasInt* ihi,
  float* H, const BlasInt* ldH, 
  float* wr, float* wi,
  float* Z, const BlasInt* ldZ, 
  float* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(dhseqr)
( const char* job, const char* compZ, const BlasInt* n, 
  const BlasInt* ilo, const BlasInt* ihi,
  double* H, const BlasInt* ldH, 
  double* wr, double* wi,
  double* Z, const BlasInt* ldZ, 
  double* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(chseqr)
( const char* job, const char* compZ, const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi,
  scomplex* H, const BlasInt* ldH,
  scomplex* w,
  scomplex* Z, const BlasInt* ldZ,
  scomplex* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(zhseqr)
( const char* job, const char* compZ, const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi,
  dcomplex* H, const BlasInt* ldH,
  dcomplex* w,
  dcomplex* Z, const BlasInt* ldZ,
  dcomplex* work, const BlasInt* workSize,
  BlasInt* info );

// Compute eigenpairs of a general matrix using the QR algorithm followed
// by a sequence of careful triangular solves
void EL_LAPACK(sgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n, 
  float* A, const BlasInt* ldA,
  float* wr, float* wi, 
  float* VLPacked, const BlasInt* ldVL,
  float* VRPacked, const BlasInt* ldVR,
  float* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(dgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n, 
  double* A, const BlasInt* ldA,
  double* wr, double* wi, 
  double* VLPacked, const BlasInt* ldVL,
  double* VRPacked, const BlasInt* ldVR,
  double* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(cgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n,
  scomplex* A, const BlasInt* ldA,
  scomplex* w,
  scomplex* VL, const BlasInt* ldVL,
  scomplex* VR, const BlasInt* ldVR,
  scomplex* work, const BlasInt* workSize,
  float* rWork,
  BlasInt* info );
void EL_LAPACK(zgeev)
( const char* jobVL, const char* jobVR, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA,
  dcomplex* w,
  dcomplex* VL, const BlasInt* ldVL,
  dcomplex* VR, const BlasInt* ldVR,
  dcomplex* work, const BlasInt* workSize,
  double* rWork,
  BlasInt* info );

} // extern "C"

namespace El {
namespace lapack {

// Machine constants
// =================

template<>
float MachineEpsilon<float>()
{
    const char cmach = 'E';
    static float eps = EL_LAPACK(slamch)( &cmach ); 
    return eps;
}

template<> 
double MachineEpsilon<double>()
{
    const char cmach = 'E';
    static double eps = EL_LAPACK(dlamch)( &cmach );
    return eps;
}

template<> 
float MachineSafeMin<float>()
{
    const char cmach = 'S';
    static float safeMin = EL_LAPACK(slamch)( &cmach );
    return safeMin;
}

template<> 
double MachineSafeMin<double>()
{
    const char cmach = 'S';
    static double safeMin = EL_LAPACK(dlamch)( &cmach );
    return safeMin;
}

template<> 
float MachineBase<float>()
{
    const char cmach = 'B';
    static float base = EL_LAPACK(slamch)( &cmach );
    return base;
}

template<> 
double MachineBase<double>()
{
    const char cmach = 'B';
    static double base = EL_LAPACK(dlamch)( &cmach );
    return base;
}

template<>
float MachinePrecision<float>()
{
    const char cmach = 'P';
    static float prec = EL_LAPACK(slamch)( &cmach );
    return prec;
}

template<> 
double MachinePrecision<double>()
{
    const char cmach = 'P';
    static double prec = EL_LAPACK(dlamch)( &cmach );
    return prec;
}

template<> 
float MachineUnderflowExponent<float>()
{
    const char cmach = 'M';
    static float underExp = EL_LAPACK(slamch)( &cmach );
    return underExp;
}

template<> 
double MachineUnderflowExponent<double>()
{
    const char cmach = 'M';
    static double underExp = EL_LAPACK(dlamch)( &cmach );
    return underExp;
}

template<>
float MachineUnderflowThreshold<float>()
{
    const char cmach = 'U';
    static float underThresh = EL_LAPACK(slamch)( &cmach );
    return underThresh;
}

template<> 
double MachineUnderflowThreshold<double>()
{
    const char cmach = 'U';
    static double underThresh = EL_LAPACK(dlamch)( &cmach );
    return underThresh;
}

template<> 
float MachineOverflowExponent<float>()
{
    const char cmach = 'L';
    static float overExp = EL_LAPACK(slamch)( &cmach );
    return overExp;
}

template<> 
double MachineOverflowExponent<double>()
{
    const char cmach = 'L';
    static double overExp = EL_LAPACK(dlamch)( &cmach );
    return overExp;
}

template<> 
float MachineOverflowThreshold<float>()
{
    const char cmach = 'O';
    static float overThresh = EL_LAPACK(slamch)( &cmach );
    return overThresh;
}

template<> 
double MachineOverflowThreshold<double>()
{
    const char cmach = 'O';
    static double overThresh = EL_LAPACK(dlamch)( &cmach );
    return overThresh;
}

// Safely compute norms
// ====================

template<typename Real>
Real SafeNorm( const Real& alpha, const Real& beta )
{
    Real scale = 0;
    Real scaledSquare = 1;
    UpdateScaledSquare( alpha, scale, scaledSquare );
    UpdateScaledSquare( beta, scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}
template float SafeNorm( const float& alpha, const float& beta );
#ifdef EL_HAVE_QD
template DoubleDouble SafeNorm
( const DoubleDouble& alpha, const DoubleDouble& beta );
template QuadDouble SafeNorm
( const QuadDouble& alpha, const QuadDouble& beta );
#endif
#ifdef EL_HAVE_QUAD
template Quad SafeNorm
( const Quad& alpha, const Quad& beta );
#endif
#ifdef EL_HAVE_MPC
template BigFloat SafeNorm
( const BigFloat& alpha, const BigFloat& beta );
#endif

double SafeNorm( const double& alpha, const double& beta )
{ return EL_LAPACK(dlapy2)( &alpha, &beta ); }

template<typename Real>
Real SafeNorm
( const Real& alpha,
  const Real& beta,
  const Real& gamma )
{
    Real scale = 0;
    Real scaledSquare = 1;
    UpdateScaledSquare( alpha, scale, scaledSquare );
    UpdateScaledSquare( beta, scale, scaledSquare );
    UpdateScaledSquare( gamma, scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}
template float SafeNorm
( const float& alpha,
  const float& beta,
  const float& gamma );
#ifdef EL_HAVE_QD
template DoubleDouble
SafeNorm
( const DoubleDouble& alpha,
  const DoubleDouble& beta,
  const DoubleDouble& gamma );
template QuadDouble
SafeNorm
( const QuadDouble& alpha,
  const QuadDouble& beta,
  const QuadDouble& gamma );
#endif
#ifdef EL_HAVE_QUAD
template Quad SafeNorm
( const Quad& alpha,
  const Quad& beta,
  const Quad& gamma );
#endif
#ifdef EL_HAVE_MPC
template BigFloat SafeNorm
( const BigFloat& alpha,
  const BigFloat& beta,
  const BigFloat& gamma );
#endif

double SafeNorm
( const double& alpha,
  const double& beta,
  const double& gamma )
{ return EL_LAPACK(dlapy3)( &alpha, &beta, &gamma ); }

template<typename Real>
Real SafeNorm( const Complex<Real>& alpha, const Real& beta )
{ return SafeNorm( alpha.real(), alpha.imag(), beta ); }
template float SafeNorm( const Complex<float>& alpha, const float& beta );
template double SafeNorm( const Complex<double>& alpha, const double& beta );
#ifdef EL_HAVE_QUAD
template Quad SafeNorm( const Complex<Quad>& alpha, const Quad& beta );
#endif

template<typename Real>
Real SafeNorm( const Real& alpha, const Complex<Real>& beta )
{ return SafeNorm( beta, alpha ); }
template float SafeNorm( const float& alpha, const Complex<float>& beta );
template double SafeNorm( const double& alpha, const Complex<double>& beta );
#ifdef EL_HAVE_QUAD
template Quad SafeNorm( const Quad& alpha, const Complex<Quad>& beta );
#endif

// Copy a matrix
// =============
template<typename T>
void Copy
( char uplo, BlasInt m, BlasInt n, const T* A, BlasInt lda, T* B, BlasInt ldb )
{
    if( std::toupper(uplo) == 'L' )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<m; ++i )
                B[i+j*ldb] = A[i+j*lda];
    }
    else if( std::toupper(uplo) == 'U' ) 
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                B[i+j*ldb] = A[i+j*lda];
    }
    else
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B[i+j*ldb] = A[i+j*lda];
    }
}
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Int* A, BlasInt lda, Int* B, BlasInt ldb );
#ifdef EL_HAVE_QD
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const DoubleDouble* A, BlasInt lda, DoubleDouble* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const QuadDouble* A, BlasInt lda, QuadDouble* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_QUAD
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Quad* A, BlasInt lda, Quad* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Complex<Quad>* A, BlasInt lda, Complex<Quad>* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_MPC
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const BigInt* A, BlasInt lda, BigInt* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const BigFloat* A, BlasInt lda, BigFloat* B, BlasInt ldb );
#endif

void Copy
( char uplo, BlasInt m, BlasInt n, 
  const float* A, BlasInt lda, float* B, BlasInt ldb )
{ EL_LAPACK(slacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const double* A, BlasInt lda, double* B, BlasInt ldb )
{ EL_LAPACK(dlacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb )
{ EL_LAPACK(clacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb )
{ EL_LAPACK(zlacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }

// Safely compute Givens rotations (using Demmel and Kahan's algorithm)
// ====================================================================

float Givens
( const float& phi, const float& gamma, float* c, float* s )
{ float rho; EL_LAPACK(slartg)( &phi, &gamma, c, s, &rho ); return rho; }

double Givens
( const double& phi, const double& gamma, double* c, double* s )
{ double rho; EL_LAPACK(dlartg)( &phi, &gamma, c, s, &rho ); return rho; }

scomplex Givens
( const scomplex& phi, const scomplex& gamma, float* c, scomplex* s )
{ scomplex rho; EL_LAPACK(clartg)( &phi, &gamma, c, s, &rho ); return rho; }

dcomplex Givens
( const dcomplex& phi, const dcomplex& gamma, double* c, dcomplex* s )
{ dcomplex rho; EL_LAPACK(zlartg)( &phi, &gamma, c, s, &rho ); return rho; }

template<typename Real>
Real Givens( const Real& phi, const Real& gamma, Real* c, Real* s )
{
    // TODO: Switch to the approach of LAPACK's dlartg instead of the
    //       zrotg-like implementation
    return blas::Givens( phi, gamma, c, s );
}
template<typename Real>
Complex<Real> Givens
( const Complex<Real>& phi,
  const Complex<Real>& gamma,
  Real* c,
  Complex<Real>* s )
{
    // TODO: Switch to the approach of LAPACK's zlartg instead of the
    //       zrotg-like implementation
    return blas::Givens( phi, gamma, c, s );
}
#ifdef EL_HAVE_QD
template DoubleDouble Givens
( const DoubleDouble& phi,
  const DoubleDouble& gamma,
  DoubleDouble* c,
  DoubleDouble* s );
template QuadDouble Givens
( const QuadDouble& phi,
  const QuadDouble& gamma,
  QuadDouble* c,
  QuadDouble* s );
#endif
#ifdef EL_HAVE_QUAD
template Quad Givens
( const Quad& phi, const Quad& gamma,
  Quad* c,
  Quad* s );
template Complex<Quad> Givens
( const Complex<Quad>& phi,
  const Complex<Quad>& gamma,
  Quad* c,
  Complex<Quad>* s );
#endif
#ifdef EL_HAVE_MPC
template BigFloat Givens
( const BigFloat& phi,
  const BigFloat& gamma,
  BigFloat* c,
  BigFloat* s );
#endif

// Compute the EVD of a symmetric tridiagonal matrix
// =================================================

BlasInt SymmetricTridiagEigWrapper
( char job, char range, BlasInt n, float* d, float* e, float vl, float vu,
  BlasInt il, BlasInt iu, float absTol, float* w, float* Z, BlasInt ldZ )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEigWrapper"));
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    float workDummy;
    EL_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<float> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sstevr's failed");
    return m;
}

BlasInt SymmetricTridiagEigWrapper
( char job, char range, BlasInt n, double* d, double* e, double vl, double vu,
  BlasInt il, BlasInt iu, double absTol, double* w, double* Z, BlasInt ldZ )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEigWrapper"));
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    double workDummy;
    EL_LAPACK(dstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<double> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(dstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dstevr's failed");
    return m;
}

// Compute eigenvalues
// -------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'A', n, d, e, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'A', n, d, e, 0, 0, 0, 0, absTol, w, 0, 1 );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float vl, float vu, 
  float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'N', 'V', n, d, e, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double vl, double vu, 
  double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'N', 'V', n, d, e, vl, vu, 0, 0, absTol, w, 0, 1 );
}

// Index range
// ^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, BlasInt il, BlasInt iu, 
  float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, BlasInt il, BlasInt iu, 
  double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}

// Compute eigenpairs
// ------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ, 
  float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'A', n, d, e, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ, 
  double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'A', n, d, e, 0, 0, 0, 0, absTol, w, Z, ldZ );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'V', 'V', n, d, e, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'V', 'V', n, d, e, vl, vu, 0, 0, absTol, w, Z, ldZ );
}

// Index range
// ^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ, 
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}

// Compute the EVD of a Hermitian matrix
// =====================================

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, float* A, BlasInt ldA, 
  float vl, float vu, BlasInt il, BlasInt iu, float absTol, 
  float* w, float* Z, BlasInt ldZ )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    float workDummy;
    EL_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<float> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("ssyevr's failed");
    return m;
}

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, double* A, BlasInt ldA, 
  double vl, double vu, BlasInt il, BlasInt iu, double absTol, 
  double* w, double* Z, BlasInt ldZ )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    double workDummy;
    EL_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, 
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<double> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dsyevr's failed");
    return m;
}

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, scomplex* A, BlasInt ldA, 
  float vl, float vu, BlasInt il, BlasInt iu, float absTol, 
  float* w, scomplex* Z, BlasInt ldZ )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, rWorkSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    float rWorkDummy;
    scomplex workDummy;
    EL_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &rWorkDummy, &rWorkSize,
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy.real();
    rWorkSize = rWorkDummy;
    iWorkSize = iWorkDummy;
    vector<scomplex> work(workSize);
    vector<float> rWork(rWorkSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      rWork.data(), &rWorkSize, iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cheevr's failed");
    return m;
}

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, dcomplex* A, BlasInt ldA, 
  double vl, double vu, BlasInt il, BlasInt iu, double absTol, 
  double* w, dcomplex* Z, BlasInt ldZ )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, rWorkSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    double rWorkDummy;
    dcomplex workDummy;
    EL_LAPACK(zheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &rWorkDummy, &rWorkSize,
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy.real();
    rWorkSize = rWorkDummy;
    iWorkSize = iWorkDummy;
    vector<dcomplex> work(workSize);
    vector<double> rWork(rWorkSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(zheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      rWork.data(), &rWorkSize, iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zheevr's failed");
    return m;
}

// Compute the eigenvalues
// -----------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, 
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w, 
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w, 
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w, 
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, 
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w, 
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w, 
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w, 
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}

// Compute the eigenpairs
// ----------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, float* Z, BlasInt ldZ, 
  float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA, double* w, double* Z, BlasInt ldZ,
  double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA, float* w, scomplex* Z, BlasInt ldZ,
  float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA, double* w, dcomplex* Z, BlasInt ldZ,
  double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt HermitianEig
( char uplo, BlasInt n, 
  float* A, BlasInt ldA, float* w, float* Z, BlasInt ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA, double* w, double* Z, BlasInt ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA, float* w, scomplex* Z, BlasInt ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA, double* w, dcomplex* Z, BlasInt ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, 
  float* A, BlasInt ldA, float* w, float* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA, double* w, double* Z, BlasInt ldZ, 
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA, float* w, scomplex* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA, double* w, dcomplex* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_ONLY(CSE cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}

// Bidiagonal DQDS for singular values
// ===================================

void BidiagDQDS( BlasInt n, float* d, float* e )
{
    DEBUG_ONLY(CSE cse("lapack::BidiagDQDS"))
    BlasInt info;
    vector<float> work( 4*n );
    EL_LAPACK(slasq1)( &n, d, e, work.data(), &info );
    if( info != 0 )
    {
        ostringstream msg;
        if( info < 0 )
            msg << "Argument " << -info << " had an illegal value";
        else if( info == 1 )
            msg << "A split was marked in a positive value in E";
        else if( info == 2 )
            msg << "Current block of Z not bidiagonalized after 30*k its";
        else if( info == 3 )
            msg << "Termination criterion of outer while loop not met";
        RuntimeError( msg.str() );
    } 
}

void BidiagDQDS( BlasInt n, double* d, double* e )
{
    DEBUG_ONLY(CSE cse("lapack::BidiagDQDS"))
    BlasInt info;
    vector<double> work( 4*n );
    EL_LAPACK(dlasq1)( &n, d, e, work.data(), &info );
    if( info != 0 )
    {
        ostringstream msg;
        if( info < 0 )
            msg << "Argument " << -info << " had an illegal value";
        else if( info == 1 )
            msg << "A split was marked in a positive value in E";
        else if( info == 2 )
            msg << "Current block of Z not bidiagonalized after 30*k its";
        else if( info == 3 )
            msg << "Termination criterion of outer while loop not met";
        RuntimeError( msg.str() );
    } 
}

// Bidiagonal QR algorithm for SVD
// ===============================

void BidiagQRAlg
( char uplo, BlasInt n, BlasInt numColsVT, BlasInt numRowsU,
  float* d, float* e, float* VTrans, BlasInt ldVT, float* U, BlasInt ldU )
{
    DEBUG_ONLY(CSE cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    BlasInt info;
    float* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<float> work( 4*n );
    EL_LAPACK(sbdsqr)
    ( &uplo, &n, &numColsVT, &numRowsU, &numColsC,
      d, e,
      VTrans, &ldVT, U, &ldU, C, &ldC,
      work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, BlasInt n, BlasInt numColsVT, BlasInt numRowsU, 
  double* d, double* e, double* VTrans, BlasInt ldVT, double* U, BlasInt ldU )
{
    DEBUG_ONLY(CSE cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    BlasInt info;
    double* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<double> work( 4*n );
    EL_LAPACK(dbdsqr)
    ( &uplo, &n, &numColsVT, &numRowsU, &numColsC,
      d, e, 
      VTrans, &ldVT, U, &ldU, C, &ldC,
      work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, BlasInt n, BlasInt numColsVH, BlasInt numRowsU, 
  float* d, float* e,
  scomplex* VH, BlasInt ldVH,
  scomplex* U,  BlasInt ldU )
{
    DEBUG_ONLY(CSE cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    BlasInt info;
    scomplex* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<float> realWork;
    if( numColsVH == 0 && numRowsU == 0 && numColsC == 0 )
        realWork.resize( 2*n );
    else
        realWork.resize( Max(1,4*n-4) );
    EL_LAPACK(cbdsqr)
    ( &uplo, &n, &numColsVH, &numRowsU, &numColsC,
      d, e,
      VH, &ldVH, U, &ldU, C, &ldC,
      realWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, BlasInt n, BlasInt numColsVH, BlasInt numRowsU, 
  double* d, double* e,
  dcomplex* VH, BlasInt ldVH,
  dcomplex* U,  BlasInt ldU )
{
    DEBUG_ONLY(CSE cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    BlasInt info;
    dcomplex* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<double> realWork;
    if( numColsVH == 0 && numRowsU == 0 && numColsC == 0 )
        realWork.resize( 2*n );
    else
        realWork.resize( Max(1,4*n-4) );
    EL_LAPACK(zbdsqr)
    ( &uplo, &n, &numColsVH, &numRowsU, &numColsC,
      d, e,
      VH, &ldVH, U, &ldU, C, &ldC,
      realWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zbdsqr had ",info," elements of e not converge");
}

// Divide and Conquer SVD
// ======================

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  float* A, BlasInt ldA, 
  float* s,
  float* U, BlasInt ldu,
  float* VTrans, BlasInt ldvt,
  bool thin )
{
    DEBUG_ONLY(CSE cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    float workDummy;
    const BlasInt k = Min(m,n);
    vector<BlasInt> iWork(8*k);

    EL_LAPACK(sgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      &workDummy, &workSize,
      iWork.data(),
      &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(sgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      work.data(), &workSize,
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesdd's updating process failed");
}

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  double* A, BlasInt ldA, 
  double* s,
  double* U, BlasInt ldu,
  double* VTrans, BlasInt ldvt,
  bool thin )
{
    DEBUG_ONLY(CSE cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    double workDummy;
    const BlasInt k = Min(m,n);
    vector<BlasInt> iWork(8*k);

    EL_LAPACK(dgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      &workDummy, &workSize,
      iWork.data(),
      &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      work.data(), &workSize,
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesdd's updating process failed");
}

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  scomplex* A, BlasInt ldA, 
  float* s,
  scomplex* U, BlasInt ldu,
  scomplex* VH, BlasInt ldva,
  bool thin )
{
    DEBUG_ONLY(CSE cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    const BlasInt k = Min(m,n);
    const BlasInt K = Max(m,n);
    const BlasInt rWorkSize = k*Max(5*k+7,2*K+2*k+1);
    vector<float> rWork(rWorkSize);
    vector<BlasInt> iWork(8*k);

    scomplex workDummy;
    EL_LAPACK(cgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      &workDummy, &workSize,
      rWork.data(),
      iWork.data(),
      &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(cgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      work.data(), &workSize,
      rWork.data(),
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesdd's updating process failed");
}

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  dcomplex* A, BlasInt ldA, 
  double* s,
  dcomplex* U, BlasInt ldu,
  dcomplex* VH, BlasInt ldva,
  bool thin )
{
    DEBUG_ONLY(CSE cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    dcomplex workDummy;
    const BlasInt k = Min(m,n);
    const BlasInt K = Max(m,n);
    const BlasInt rWorkSize = k*Max(5*k+7,2*K+2*k+1);
    vector<double> rWork(rWorkSize);
    vector<BlasInt> iWork(8*k);

    EL_LAPACK(zgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      &workDummy, &workSize,
      rWork.data(),
      iWork.data(),
      &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      work.data(), &workSize,
      rWork.data(),
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesdd's updating process failed");
}

// QR-algorithm SVD
// ================

void QRSVD
( BlasInt m, BlasInt n,
  float* A, BlasInt ldA, 
  float* s,
  float* U, BlasInt ldu,
  float* VTrans, BlasInt ldvt,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_ONLY(CSE cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVT= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    float workDummy;

    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      &workDummy, &workSize,
      &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void QRSVD
( BlasInt m, BlasInt n,
  double* A, BlasInt ldA, 
  double* s,
  double* U, BlasInt ldu,
  double* VTrans, BlasInt ldvt,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_ONLY(CSE cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVT= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    double workDummy;

    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      &workDummy, &workSize,
      &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void QRSVD
( BlasInt m, BlasInt n,
  scomplex* A, BlasInt ldA, 
  float* s,
  scomplex* U, BlasInt ldu,
  scomplex* VH, BlasInt ldva,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_ONLY(CSE cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVH= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    const BlasInt k = Min(m,n);
    vector<float> rWork(5*k);

    scomplex workDummy;
    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      &workDummy, &workSize,
      rWork.data(),
      &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      work.data(), &workSize,
      rWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void QRSVD
( BlasInt m, BlasInt n,
  dcomplex* A, BlasInt ldA, 
  double* s,
  dcomplex* U, BlasInt ldu,
  dcomplex* VH, BlasInt ldva,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_ONLY(CSE cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVH= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    dcomplex workDummy;
    const BlasInt k = Min(m,n);
    vector<double> rWork(5*k);

    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      &workDummy, &workSize,
      rWork.data(),
      &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      work.data(), &workSize,
      rWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Compute singular values (with DQDS)
// ===================================

void SVD( BlasInt m, BlasInt n, float* A, BlasInt ldA, float* s )
{
    DEBUG_ONLY(CSE cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVT='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    float workDummy;

    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void SVD( BlasInt m, BlasInt n, double* A, BlasInt ldA, double* s )
{
    DEBUG_ONLY(CSE cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVT='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    double workDummy;

    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void SVD( BlasInt m, BlasInt n, scomplex* A, BlasInt ldA, float* s )
{
    DEBUG_ONLY(CSE cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVH='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    scomplex workDummy;
    const BlasInt k = Min(m,n);
    vector<float> rWork(5*k);

    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void SVD( BlasInt m, BlasInt n, dcomplex* A, BlasInt ldA, double* s )
{
    DEBUG_ONLY(CSE cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVH='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    dcomplex workDummy;
    const BlasInt k = Min(m,n);
    vector<double> rWork(5*k);

    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Compute the Schur decomposition of an upper Hessenberg matrix
// =============================================================

void HessenbergSchur
( BlasInt n, float* H, BlasInt ldH, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1, workSize=-1, info;
    float workDummy;
    vector<float> wr( n ), wi( n );
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<float>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n, double* H, BlasInt ldH, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1, workSize=-1, info;
    double workDummy;
    vector<double> wr( n ), wi( n );
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");
    
    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<double>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n, scomplex* H, BlasInt ldH, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1, workSize=-1, info;
    scomplex workDummy;
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void HessenbergSchur
( BlasInt n, dcomplex* H, BlasInt ldH, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1, workSize=-1, info;
    dcomplex workDummy;
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

void HessenbergSchur
( BlasInt n, float* H, BlasInt ldH, scomplex* w, float* Q, BlasInt ldQ, 
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    BlasInt ilo=1, ihi=n;
    BlasInt workSize=-1, info;
    float workDummy;
    vector<float> wr( n ), wi( n );
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<float>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n, double* H, BlasInt ldH, dcomplex* w, double* Q, BlasInt ldQ, 
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    BlasInt ilo=1, ihi=n;
    BlasInt workSize=-1, info;
    double workDummy;
    vector<double> wr( n ), wi( n );
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");
    
    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<double>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n, scomplex* H, BlasInt ldH, scomplex* w, scomplex* Q, BlasInt ldQ,
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    BlasInt ilo=1, ihi=n;
    BlasInt workSize=-1, info;
    scomplex workDummy;
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void HessenbergSchur
( BlasInt n, dcomplex* H, BlasInt ldH, dcomplex* w, dcomplex* Q, BlasInt ldQ,
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    BlasInt ilo=1, ihi=n;
    BlasInt workSize=-1, info;
    dcomplex workDummy;
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

// Compute eigenvalues/pairs of an upper Hessenberg matrix
// =======================================================

void HessenbergEig( BlasInt n, float* H, BlasInt ldH, scomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( BlasInt n, double* H, BlasInt ldH, dcomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( BlasInt n, scomplex* H, BlasInt ldH, scomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( BlasInt n, dcomplex* H, BlasInt ldH, dcomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

// TODO: Compute eigenpairs

// Compute the Schur decomposition of a square matrix
// ==================================================

void Schur( BlasInt n, float* A, BlasInt ldA, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    vector<float> tau( n );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    vector<float> wr( n ), wi( n );
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<float> work( workSize );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<float>(wr[i],wi[i]);
}

void Schur( BlasInt n, double* A, BlasInt ldA, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    vector<double> tau( n );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    vector<double> wr( n ), wi( n );
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<double> work( workSize );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<double>(wr[i],wi[i]);
}

void Schur
( BlasInt n, scomplex* A, BlasInt ldA, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    vector<scomplex> tau( n );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<scomplex> work( workSize );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Schur
( BlasInt n, dcomplex* A, BlasInt ldA, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    vector<dcomplex> tau( n );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<dcomplex> work( workSize );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

void Schur
( BlasInt n, float* A, BlasInt ldA, scomplex* w, float* Q, BlasInt ldQ, 
  bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    vector<float> tau( n );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    vector<float> wr( n ), wi( n );
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<float> work( workSize );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<float>(wr[i],wi[i]);
}

void Schur
( BlasInt n, double* A, BlasInt ldA, dcomplex* w, double* Q, BlasInt ldQ, 
  bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    vector<double> tau( n );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    vector<double> wr( n ), wi( n );
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<double> work( workSize );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = El::Complex<double>(wr[i],wi[i]);
}

void Schur
( BlasInt n, scomplex* A, BlasInt ldA, scomplex* w, scomplex* Q, BlasInt ldQ, 
  bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    vector<scomplex> tau( n );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<scomplex> work( workSize );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, work.data(), &workSize, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Schur
( BlasInt n, dcomplex* A, BlasInt ldA, dcomplex* w, dcomplex* Q, BlasInt ldQ, 
  bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    vector<dcomplex> tau( n );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<dcomplex> work( workSize );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, work.data(), &workSize, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

// Compute the eigenvalues/pairs of a square matrix
// ================================================

// Eigenvalues only
// ----------------

void Eig( BlasInt n, float* A, BlasInt ldA, scomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

void Eig( BlasInt n, double* A, BlasInt ldA, dcomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

void Eig( BlasInt n, scomplex* A, BlasInt ldA, scomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

void Eig( BlasInt n, dcomplex* A, BlasInt ldA, dcomplex* w )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

// Eigenpairs
// ----------
// NOTE: When the matrices are real, an BlasInterface is also provided which returns
//       a packing of the eigenvectors which exploits the fact that, if the
//       eigenvalue is real, so is the corresponding eigenvector, otherwise
//       the eigenvalue's complex conjugate is also an eigenvalue, and the 
//       corresponding eigenvector is also the conjugate. Thus, an n x n
//       real matrix can be used to represent the eigenvectors if
//           x(j  ) = X(:,j) + X(:,j+1)*1i,
//           x(j+1) = X(:,j) - X(:,j+1)*1i
//       when the j'th and j+1'th eigenvalues are complex conjugates.

void Eig
( BlasInt n, float* A, BlasInt ldA, scomplex* w, float* XPacked, BlasInt ldX )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    vector<float> wReal(n), wImag(n);
    BlasInt workSize=-1, info;
    float workDummy;
    EL_LAPACK(sgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim, 
      XPacked, &ldX, &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<float> work( workSize );
    EL_LAPACK(sgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim, 
      XPacked, &ldX, work.data(), &workSize, &info );

    // Post-process the eigenvalues
    for( Int j=0; j<n; ++j )
        w[j] = Complex<float>(wReal[j],wImag[j]);
}

void Eig
( BlasInt n, double* A, BlasInt ldA, dcomplex* w, double* XPacked, BlasInt ldX )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    vector<double> wReal(n), wImag(n);
    BlasInt workSize=-1, info;
    double workDummy;
    EL_LAPACK(dgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim,
      XPacked, &ldX, &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<double> work( workSize );
    EL_LAPACK(dgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim,
      XPacked, &ldX, work.data(), &workSize, &info );

    // Post-process the eigenvalues
    for( Int j=0; j<n; ++j )
        w[j] = Complex<double>(wReal[j],wImag[j]);
}

void Eig
( BlasInt n, float* A, BlasInt ldA, scomplex* w, scomplex* X, BlasInt ldX )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    float* XPacked = (float*)X;    
    Eig( n, A, ldA, w, XPacked, ldX );
    // Unpack the eigenvectors
    vector<scomplex> z(n);
    Int j=n-1;
    while( j >= 0 )
    {
        const bool inPair = ( w[j].imag() != float(0) );
        if( inPair )
        {
            for( Int i=0; i<n; ++i )
                z[i] = XPacked[i+(j-1)*ldX] + XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i ) 
            {
                X[i+(j-1)*ldX] =      z[i];
                X[i+ j   *ldX] = Conj(z[i]);
            }
            j -= 2;
        }
        else
        {
            for( Int i=0; i<n; ++i ) 
                z[i] = XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i )
                X[i+j*ldX] = z[i];
            j -= 1; 
        }
    }
}

void Eig
( BlasInt n, double* A, BlasInt ldA, dcomplex* w, dcomplex* X, BlasInt ldX )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    double* XPacked = (double*)X;    
    Eig( n, A, ldA, w, XPacked, ldX );
    // Unpack the eigenvectors
    vector<scomplex> z(n);
    Int j=n-1;
    while( j >= 0 )
    {
        const bool inPair = ( w[j].imag() != double(0) );
        if( inPair )
        {
            for( Int i=0; i<n; ++i )
                z[i] = XPacked[i+(j-1)*ldX] + XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i ) 
            {
                X[i+(j-1)*ldX] =      z[i];
                X[i+ j   *ldX] = Conj(z[i]);
            }
            j -= 2;
        }
        else
        {
            for( Int i=0; i<n; ++i ) 
                z[i] = XPacked[i+j*ldX];
            for( Int i=0; i<n; ++i )
                X[i+j*ldX] = z[i];
            j -= 1; 
        }
    }
}

void Eig
( BlasInt n, scomplex* A, BlasInt ldA, scomplex* w, scomplex* X, BlasInt ldX )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    vector<float> rWork( 2*n );
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    BlasInt workSize=-1, info;
    scomplex workDummy;
    EL_LAPACK(cgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<scomplex> work( workSize );
    EL_LAPACK(cgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX, 
      work.data(), &workSize, rWork.data(), &info );
}

void Eig
( BlasInt n, dcomplex* A, BlasInt ldA, dcomplex* w, dcomplex* X, BlasInt ldX )
{
    DEBUG_ONLY(CSE cse("lapack::Eig"))
    vector<double> rWork( 2*n );
    const char jobVL='N', jobVR='V';
    const BlasInt fakeLDim = 1;

    BlasInt workSize=-1, info;
    dcomplex workDummy;
    EL_LAPACK(zgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX,
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<dcomplex> work( workSize );
    EL_LAPACK(zgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX,
      work.data(), &workSize, rWork.data(), &info );
}

// TODO: Return the left eigenvectors?

} // namespace lapack
} // namespace El
