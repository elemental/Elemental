/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

using elem::scomplex;
using elem::dcomplex;

extern "C" {

// Machine constants
float ELEM_LAPACK(slamch)( const char* cmach );
double ELEM_LAPACK(dlamch)( const char* cmach );

// Safe norms
float  ELEM_LAPACK(slapy2)( const float * alpha, const float * beta );
double ELEM_LAPACK(dlapy2)( const double* alpha, const double* beta );
float ELEM_LAPACK(slapy3)
( const float * alpha, const float * beta, const float * gamma );
double ELEM_LAPACK(dlapy3)
( const double* alpha, const double* beta, const double* gamma );

// Safely compute a Givens rotation
void ELEM_LAPACK(slartg)
( const float* phi, const float* gamma,
  float* c, float* s, float* rho );
void ELEM_LAPACK(dlartg)
( const double* phi, const double* gamma,
  double* c, double* s, double* rho );
void ELEM_LAPACK(clartg)
( const scomplex* phi, const scomplex* gamma,
  float* c, scomplex* s, scomplex* rho );
void ELEM_LAPACK(zlartg)
( const dcomplex* phi, const dcomplex* gamma,
  double* c, dcomplex* s, dcomplex* rho );

// Symmetric tridiagonal eigensolvers (via MRRR)
void ELEM_LAPACK(sstevr)
( const char* job, const char* range, const int* n,
  float* d, float* e, const float* vl, const float* vu, 
  const int* il, const int* iu, const float* absTol, int * m, 
  float* w, float* Z, const int* ldZ, int* isuppZ, 
  float* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );
void ELEM_LAPACK(dstevr)
( const char* job, const char* range, const int* n,
  double* d, double* e, const double* vl, const double* vu, 
  const int* il, const int* iu, const double* absTol, int * m, 
  double* w, double* Z, const int* ldZ, int* isuppZ, 
  double* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );

// Hermitian eigensolvers (via MRRR)
void ELEM_LAPACK(ssyevr)
( const char* job, const char* range, const char* uplo, const int* n,
  float* A, const int* ldA, const float* vl, const float* vu, 
  const int* il, const int* iu, const float* absTol, int * m, 
  float* w, float* Z, const int* ldZ, int* isuppZ, 
  float* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );
void ELEM_LAPACK(dsyevr)
( const char* job, const char* range, const char* uplo, const int* n,
  double* A, const int* ldA, const double* vl, const double* vu, 
  const int* il, const int* iu, const double* absTol, int * m, 
  double* w, double* Z, const int* ldZ, int* isuppZ, 
  double* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );
void ELEM_LAPACK(cheevr)
( const char* job, const char* range, const char* uplo, const int* n,
  scomplex* A, const int* ldA, const float* vl, const float* vu, 
  const int* il, const int* iu, const float* absTol, int* m,
  float* w, scomplex* Z, const int* ldZ, int* isuppZ, 
  scomplex* work, const int* workSize, float* rWork, const int* rWorkSize, 
  int* iWork, const int* iWorkSize, int* info );
void ELEM_LAPACK(zheevr)
( const char* job, const char* range, const char* uplo, const int* n,
  dcomplex* A, const int* ldA, const double* vl, const double* vu, 
  const int* il, const int* iu, const double* absTol, int* m,
  double* w, dcomplex* Z, const int* ldZ, int* isuppZ, 
  dcomplex* work, const int* workSize, double* rWork, const int* rWorkSize, 
  int* iWork, const int* iWorkSize, int* info );

// Bidiagonal DQDS
void ELEM_LAPACK(slasq1)
( const int* n, float* d, float* e, float* work, int* info );
void ELEM_LAPACK(dlasq1)
( const int* n, double* d, double* e, double* work, int* info );

// Bidiagonal QR
void ELEM_LAPACK(sbdsqr)
( const char* uplo, const int* n, const int* numColsVT, const int* numRowsU,
  const int* numColsC, float* d, float* e, float* VTrans, const int* ldVT,
  float* U, const int* ldU, float* C, const int* ldC, float* work, int* info );
void ELEM_LAPACK(dbdsqr)
( const char* uplo, const int* n, const int* numColsVT, const int* numRowsU,
  const int* numColsC, double* d, double* e,
  double* VTrans, const int* ldVT, double* U, const int* ldU,
  double* C, const int* ldC, double* work, int* info );
void ELEM_LAPACK(cbdsqr)
( const char* uplo, const int* n, const int* numColsVH, const int* numRowsU,
  const int* numColsC, float* d, float* e,
  scomplex* VH, const int* ldVH, scomplex* U, const int* ldU,
  scomplex* C, const int* ldC, float* work, int* info );
void ELEM_LAPACK(zbdsqr)
( const char* uplo, const int* n, const int* numColsVH, const int* numRowsU,
  const int* numColsC, double* d, double* e,
  dcomplex* VH, const int* ldVH, dcomplex* U, const int* ldU,
  dcomplex* C, const int* ldC, double* work, int* info );

// Divide and Conquer SVD
void ELEM_LAPACK(sgesdd)
( const char* jobz, const int* m, const int* n, float* A, const int* ldA,
  float* s, float* U, const int* ldu, float* VTrans, const int* ldvt,
  float* work, const int* workSize, int* iWork, int* info );
void ELEM_LAPACK(dgesdd)
( const char* jobz, const int* m, const int* n, double* A, const int* ldA,
  double* s, double* U, const int* ldu, double* VTrans, const int* ldvt,
  double* work, const int* workSize, int* iWork, int* info );
void ELEM_LAPACK(cgesdd)
( const char* jobz, const int* m, const int* n,
  scomplex* A, const int* ldA, float* s,
  scomplex* U, const int* ldu, scomplex* VTrans, const int* ldvt,
  scomplex* work, const int* workSize, float* rWork,
  int* iWork, int* info );
void ELEM_LAPACK(zgesdd)
( const char* jobz, const int* m, const int* n,
  dcomplex* A, const int* ldA, double* s,
  dcomplex* U, const int* ldu, dcomplex* VH, const int* ldva,
  dcomplex* work, const int* workSize, double* rWork,
  int* iWork, int* info );

// QR-algorithm SVD [DQDS when no singular vectors desired]
void ELEM_LAPACK(sgesvd)
( const char* jobU, const char* jobVT, const int* m, const int* n,
  float* A, const int* ldA,
  float* s, float* U, const int* ldu, float* VTrans, const int* ldvt,
  float* work, const int* workSize, int* info );
void ELEM_LAPACK(dgesvd)
( const char* jobU, const char* jobVT, const int* m, const int* n,
  double* A, const int* ldA,
  double* s, double* U, const int* ldu, double* VTrans, const int* ldvt,
  double* work, const int* workSize, int* info );
void ELEM_LAPACK(cgesvd)
( const char* jobU, const char* jobVH, const int* m, const int* n,
  scomplex* A, const int* ldA, float* s,
  scomplex* U, const int* ldu, scomplex* VTrans, const int* ldvt,
  scomplex* work, const int* workSize, float* rWork, int* info );
void ELEM_LAPACK(zgesvd)
( const char* jobU, const char* jobVH, const int* m, const int* n,
  dcomplex* A, const int* ldA, double* s,
  dcomplex* U, const int* ldu, dcomplex* VH, const int* ldva,
  dcomplex* work, const int* workSize, double* rWork, int* info );

// Reduction to Hessenberg form
void ELEM_LAPACK(sgehrd)
( const int* n, const int* ilo, const int* ihi, float* A, const int* ldA,
  float* tau, float* work, const int* workSize, int* info );
void ELEM_LAPACK(dgehrd)
( const int* n, const int* ilo, const int* ihi, double* A, const int* ldA,
  double* tau, double* work, const int* workSize, int* info );
void ELEM_LAPACK(cgehrd)
( const int* n, const int* ilo, const int* ihi, scomplex* A, const int* ldA,
  scomplex* tau, scomplex* work, const int* workSize, int* info );
void ELEM_LAPACK(zgehrd)
( const int* n, const int* ilo, const int* ihi, dcomplex* A, const int* ldA,
  dcomplex* tau, dcomplex* work, const int* workSize, int* info );

// Generates a unitary matrix defined as the product of Householder reflectors
void ELEM_LAPACK(sorghr)
( const int* n, const int* ilo, const int* ihi, float* A, const int* ldA,
  const float* tau, float* work, const int* workSize, int* info );
void ELEM_LAPACK(dorghr)
( const int* n, const int* ilo, const int* ihi, double* A, const int* ldA,
  const double* tau, double* work, const int* workSize, int* info );
void ELEM_LAPACK(cunghr)
( const int* n, const int* ilo, const int* ihi, scomplex* A, const int* ldA,
  const scomplex* tau, scomplex* work, const int* workSize, int* info );
void ELEM_LAPACK(zunghr)
( const int* n, const int* ilo, const int* ihi, dcomplex* A, const int* ldA,
  const dcomplex* tau, dcomplex* work, const int* workSize, int* info );

// Hessenberg QR algorithm
void ELEM_LAPACK(shseqr)
( const char* job, const char* compZ, const int* n, 
  const int* ilo, const int* ihi, float* H, const int* ldH, 
  float* wr, float* wi, float* Z, const int* ldZ, 
  float* work, const int* workSize, int* info );
void ELEM_LAPACK(dhseqr)
( const char* job, const char* compZ, const int* n, 
  const int* ilo, const int* ihi, double* H, const int* ldH, 
  double* wr, double* wi, double* Z, const int* ldZ, 
  double* work, const int* workSize, int* info );
void ELEM_LAPACK(chseqr)
( const char* job, const char* compZ, const int* n,
  const int* ilo, const int* ihi, scomplex* H, const int* ldH,
  scomplex* w, scomplex* Z, const int* ldZ,
  scomplex* work, const int* workSize, int* info );
void ELEM_LAPACK(zhseqr)
( const char* job, const char* compZ, const int* n,
  const int* ilo, const int* ihi, dcomplex* H, const int* ldH,
  dcomplex* w, dcomplex* Z, const int* ldZ,
  dcomplex* work, const int* workSize, int* info );

// Compute eigenpairs of a general matrix using the QR algorithm followed
// by a sequence of careful triangular solves
void ELEM_LAPACK(sgeev)
( const char* jobVL, const char* jobVR, const int* n, 
  float* A, const int* ldA, float* wr, float* wi, 
  float* VLPacked, const int* ldVL, float* VRPacked, const int* ldVR,
  float* work, const int* workSize, int* info );
void ELEM_LAPACK(dgeev)
( const char* jobVL, const char* jobVR, const int* n, 
  double* A, const int* ldA, double* wr, double* wi, 
  double* VLPacked, const int* ldVL, double* VRPacked, const int* ldVR,
  double* work, const int* workSize, int* info );
void ELEM_LAPACK(cgeev)
( const char* jobVL, const char* jobVR, const int* n,
  scomplex* A, const int* ldA, scomplex* w,
  scomplex* VL, const int* ldVL, scomplex* VR, const int* ldVR,
  scomplex* work, const int* workSize, float* rWork, int* info );
void ELEM_LAPACK(zgeev)
( const char* jobVL, const char* jobVR, const int* n,
  dcomplex* A, const int* ldA, dcomplex* w,
  dcomplex* VL, const int* ldVL, dcomplex* VR, const int* ldVR,
  dcomplex* work, const int* workSize, double* rWork, int* info );

} // extern "C"

namespace elem {
namespace lapack {

// Machine constants
// =================

template<>
float MachineEpsilon<float>()
{
    const char cmach = 'E';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachineEpsilon<double>()
{
    const char cmach = 'E';
    return ELEM_LAPACK(dlamch)( &cmach );
}

template<> 
float MachineSafeMin<float>()
{
    const char cmach = 'S';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachineSafeMin<double>()
{
    const char cmach = 'S';
    return ELEM_LAPACK(dlamch)( &cmach );
}

template<> 
float MachineBase<float>()
{
    const char cmach = 'B';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachineBase<double>()
{
    const char cmach = 'B';
    return ELEM_LAPACK(dlamch)( &cmach );
}

template<>
float MachinePrecision<float>()
{
    const char cmach = 'P';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachinePrecision<double>()
{
    const char cmach = 'P';
    return ELEM_LAPACK(dlamch)( &cmach );
}

template<> 
float MachineUnderflowExponent<float>()
{
    const char cmach = 'M';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachineUnderflowExponent<double>()
{
    const char cmach = 'M';
    return ELEM_LAPACK(dlamch)( &cmach );
}

template<>
float MachineUnderflowThreshold<float>()
{
    const char cmach = 'U';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachineUnderflowThreshold<double>()
{
    const char cmach = 'U';
    return ELEM_LAPACK(dlamch)( &cmach );
}

template<> 
float MachineOverflowExponent<float>()
{
    const char cmach = 'L';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachineOverflowExponent<double>()
{
    const char cmach = 'L';
    return ELEM_LAPACK(dlamch)( &cmach );
}

template<> 
float MachineOverflowThreshold<float>()
{
    const char cmach = 'O';
    return ELEM_LAPACK(slamch)( &cmach );
}

template<> 
double MachineOverflowThreshold<double>()
{
    const char cmach = 'O';
    return ELEM_LAPACK(dlamch)( &cmach );
}

// Safely compute norms
// ====================

float SafeNorm( float alpha, float beta )
{ return ELEM_LAPACK(slapy2)( &alpha, &beta ); }

double SafeNorm( double alpha, double beta )
{ return ELEM_LAPACK(dlapy2)( &alpha, &beta ); }

float SafeNorm( float alpha, float beta, float gamma )
{ return ELEM_LAPACK(slapy3)( &alpha, &beta, &gamma ); }

double SafeNorm( double alpha, double beta, double gamma )
{ return ELEM_LAPACK(dlapy3)( &alpha, &beta, &gamma ); }

float SafeNorm( Complex<float> alpha, float beta )
{ return SafeNorm( alpha.real(), alpha.imag(), beta ); }

double SafeNorm( Complex<double> alpha, double beta )
{ return SafeNorm( alpha.real(), alpha.imag(), beta ); }

float SafeNorm( float alpha, Complex<float> beta )
{ return SafeNorm( beta, alpha ); }

double SafeNorm( double alpha, Complex<double> beta )
{ return SafeNorm( beta, alpha ); }

// Safely compute Givens rotations (using Demmel and Kahan's algorithm)
// ====================================================================

float Givens( float phi, float gamma, float* c, float* s )
{ float rho; ELEM_LAPACK(slartg)( &phi, &gamma, c, s, &rho ); return rho; }

double Givens( double phi, double gamma, double* c, double* s )
{ double rho; ELEM_LAPACK(dlartg)( &phi, &gamma, c, s, &rho ); return rho; }

scomplex Givens( scomplex phi, scomplex gamma, float* c, scomplex* s )
{ scomplex rho; ELEM_LAPACK(clartg)( &phi, &gamma, c, s, &rho ); return rho; }

dcomplex Givens( dcomplex phi, dcomplex gamma, double* c, dcomplex* s )
{ dcomplex rho; ELEM_LAPACK(zlartg)( &phi, &gamma, c, s, &rho ); return rho; }

// Compute the EVD of a symmetric tridiagonal matrix
// =================================================

int SymmetricTridiagEigWrapper
( char job, char range, int n, float* d, float* e, float vl, float vu,
  int il, int iu, float absTol, float* w, float* Z, int ldZ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEigWrapper"));
    if( n == 0 )
        return 0;

    std::vector<int> isuppZ( 2*n );

    int workSize=-1, iWorkSize=-1, m, info;
    int iWorkDummy;
    float workDummy;
    ELEM_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    std::vector<float> work(workSize);
    std::vector<int> iWork(iWorkSize);
    ELEM_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sstevr's failed");
    return m;
}

int SymmetricTridiagEigWrapper
( char job, char range, int n, double* d, double* e, double vl, double vu,
  int il, int iu, double absTol, double* w, double* Z, int ldZ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEigWrapper"));
    if( n == 0 )
        return 0;

    std::vector<int> isuppZ( 2*n );

    int workSize=-1, iWorkSize=-1, m, info;
    int iWorkDummy;
    double workDummy;
    ELEM_LAPACK(dstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    std::vector<double> work(workSize);
    std::vector<int> iWork(iWorkSize);
    ELEM_LAPACK(dstevr)
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
( int n, float* d, float* e, float* w, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'A', n, d, e, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void SymmetricTridiagEig
( int n, double* d, double* e, double* w, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'A', n, d, e, 0, 0, 0, 0, absTol, w, 0, 1 );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int SymmetricTridiagEig
( int n, float* d, float* e, float* w, float vl, float vu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'N', 'V', n, d, e, vl, vu, 0, 0, absTol, w, 0, 1 );
}
int SymmetricTridiagEig
( int n, double* d, double* e, double* w, double vl, double vu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'N', 'V', n, d, e, vl, vu, 0, 0, absTol, w, 0, 1 );
}

// Index range
// ^^^^^^^^^^^^^
void SymmetricTridiagEig
( int n, float* d, float* e, float* w, int il, int iu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void SymmetricTridiagEig
( int n, double* d, double* e, double* w, int il, int iu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'N', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}

// Compute eigenpairs
// ------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void SymmetricTridiagEig
( int n, float* d, float* e, float* w, float* Z, int ldZ, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'A', n, d, e, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void SymmetricTridiagEig
( int n, double* d, double* e, double* w, double* Z, int ldZ, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'A', n, d, e, 0, 0, 0, 0, absTol, w, Z, ldZ );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int SymmetricTridiagEig
( int n, float* d, float* e, float* w, float* Z, int ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'V', 'V', n, d, e, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
int SymmetricTridiagEig
( int n, double* d, double* e, double* w, double* Z, int ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    return SymmetricTridiagEigWrapper
    ( 'V', 'V', n, d, e, vl, vu, 0, 0, absTol, w, Z, ldZ );
}

// Index range
// ^^^^^^^^^^^^^
void SymmetricTridiagEig
( int n, float* d, float* e, float* w, float* Z, int ldZ, 
  int il, int iu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void SymmetricTridiagEig
( int n, double* d, double* e, double* w, double* Z, int ldZ,
  int il, int iu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    SymmetricTridiagEigWrapper
    ( 'V', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}

// Compute the EVD of a Hermitian matrix
// =====================================

int HermitianEigWrapper
( char job, char range, char uplo, int n, float* A, int ldA, 
  float vl, float vu, int il, int iu, float absTol, 
  float* w, float* Z, int ldZ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppZ( 2*n );

    int workSize=-1, iWorkSize=-1, m, info;
    int iWorkDummy;
    float workDummy;
    ELEM_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    std::vector<float> work(workSize);
    std::vector<int> iWork(iWorkSize);
    ELEM_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("ssyevr's failed");
    return m;
}

int HermitianEigWrapper
( char job, char range, char uplo, int n, double* A, int ldA, 
  double vl, double vu, int il, int iu, double absTol, 
  double* w, double* Z, int ldZ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppZ( 2*n );

    int workSize=-1, iWorkSize=-1, m, info;
    int iWorkDummy;
    double workDummy;
    ELEM_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, 
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    std::vector<double> work(workSize);
    std::vector<int> iWork(iWorkSize);
    ELEM_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dsyevr's failed");
    return m;
}

int HermitianEigWrapper
( char job, char range, char uplo, int n, scomplex* A, int ldA, 
  float vl, float vu, int il, int iu, float absTol, 
  float* w, scomplex* Z, int ldZ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppZ( 2*n );

    int workSize=-1, rWorkSize=-1, iWorkSize=-1, m, info;
    int iWorkDummy;
    float rWorkDummy;
    scomplex workDummy;
    ELEM_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &rWorkDummy, &rWorkSize,
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy.real();
    rWorkSize = rWorkDummy;
    iWorkSize = iWorkDummy;
    std::vector<scomplex> work(workSize);
    std::vector<float> rWork(rWorkSize);
    std::vector<int> iWork(iWorkSize);
    ELEM_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      rWork.data(), &rWorkSize, iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cheevr's failed");
    return m;
}

int HermitianEigWrapper
( char job, char range, char uplo, int n, dcomplex* A, int ldA, 
  double vl, double vu, int il, int iu, double absTol, 
  double* w, dcomplex* Z, int ldZ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEigWrapper"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppZ( 2*n );

    int workSize=-1, rWorkSize=-1, iWorkSize=-1, m, info;
    int iWorkDummy;
    double rWorkDummy;
    dcomplex workDummy;
    ELEM_LAPACK(zheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &rWorkDummy, &rWorkSize,
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy.real();
    rWorkSize = rWorkDummy;
    iWorkSize = iWorkDummy;
    std::vector<dcomplex> work(workSize);
    std::vector<double> rWork(rWorkSize);
    std::vector<int> iWork(iWorkSize);
    ELEM_LAPACK(zheevr)
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
( char uplo, int n, float* A, int ldA, float* w, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int HermitianEig
( char uplo, int n, float* A, int ldA, float* w, 
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
int HermitianEig
( char uplo, int n, double* A, int ldA, double* w, 
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
int HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, 
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
int HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, 
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, int ldA, float* w, 
  int il, int iu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w, 
  int il, int iu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, 
  int il, int iu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, 
  int il, int iu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}

// Compute the eigenpairs
// ----------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, int ldA, float* w, float* Z, int ldZ, 
  float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double* Z, int ldZ,
  double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, scomplex* Z, int ldZ,
  float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, dcomplex* Z, int ldZ,
  double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int HermitianEig
( char uplo, int n, float* A, int ldA, float* w, float* Z, int ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
int HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double* Z, int ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
int HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, scomplex* Z, int ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
int HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, dcomplex* Z, int ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, int ldA, float* w, float* Z, int ldZ,
  int il, int iu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double* Z, int ldZ, 
  int il, int iu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, scomplex* Z, int ldZ,
  int il, int iu, float absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, dcomplex* Z, int ldZ,
  int il, int iu, double absTol )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}

// Bidiagonal DQDS for singular values
// ===================================

void BidiagDQDS( int n, float* d, float* e )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagDQDS"))
    int info;
    std::vector<float> work( 4*n );
    ELEM_LAPACK(slasq1)( &n, d, e, work.data(), &info );
    if( info != 0 )
    {
        std::ostringstream msg;
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

void BidiagDQDS( int n, double* d, double* e )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagDQDS"))
    int info;
    std::vector<double> work( 4*n );
    ELEM_LAPACK(dlasq1)( &n, d, e, work.data(), &info );
    if( info != 0 )
    {
        std::ostringstream msg;
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
( char uplo, int n, int numColsVT, int numRowsU,
  float* d, float* e, float* VTrans, int ldVT, float* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    float* C=0;
    const int numColsC=0, ldC=1;
    std::vector<float> work( 4*n );
    ELEM_LAPACK(sbdsqr)
    ( &uplo, &n, &numColsVT, &numRowsU, &numColsC, d, e, VTrans, &ldVT,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, int n, int numColsVT, int numRowsU, 
  double* d, double* e, double* VTrans, int ldVT, double* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    double* C=0;
    const int numColsC=0, ldC=1;
    std::vector<double> work( 4*n );
    ELEM_LAPACK(dbdsqr)
    ( &uplo, &n, &numColsVT, &numRowsU, &numColsC, d, e, VTrans, &ldVT,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, int n, int numColsVH, int numRowsU, 
  float* d, float* e, scomplex* VH, int ldVH, scomplex* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    scomplex* C=0;
    const int numColsC=0, ldC=1;
    std::vector<float> work( 4*n );
    ELEM_LAPACK(cbdsqr)
    ( &uplo, &n, &numColsVH, &numRowsU, &numColsC, d, e, VH, &ldVH,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, int n, int numColsVH, int numRowsU, 
  double* d, double* e, dcomplex* VH, int ldVH, dcomplex* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    dcomplex* C=0;
    const int numColsC=0, ldC=1;
    std::vector<double> work( 4*n );
    ELEM_LAPACK(zbdsqr)
    ( &uplo, &n, &numColsVH, &numRowsU, &numColsC, d, e, VH, &ldVH,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zbdsqr had ",info," elements of e not converge");
}

// Divide and Conquer SVD
// ======================

void DivideAndConquerSVD
( int m, int n, float* A, int ldA, 
  float* s, float* U, int ldu, float* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int workSize=-1, info;
    float workDummy;
    const int k = std::min(m,n);
    std::vector<int> iWork(8*k);

    ELEM_LAPACK(sgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, &workDummy, &workSize,
      iWork.data(), &info );

    workSize = workDummy;
    std::vector<float> work(workSize);
    ELEM_LAPACK(sgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, work.data(), &workSize,
      iWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesdd's updating process failed");
}

void DivideAndConquerSVD
( int m, int n, double* A, int ldA, 
  double* s, double* U, int ldu, double* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int workSize=-1, info;
    double workDummy;
    const int k = std::min(m,n);
    std::vector<int> iWork(8*k);

    ELEM_LAPACK(dgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, &workDummy, &workSize,
      iWork.data(), &info );

    workSize = workDummy;
    std::vector<double> work(workSize);
    ELEM_LAPACK(dgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, work.data(), &workSize,
      iWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesdd's updating process failed");
}

void DivideAndConquerSVD
( int m, int n, scomplex* A, int ldA, 
  float* s, scomplex* U, int ldu, scomplex* VH, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int workSize=-1, info;
    const int k = std::min(m,n);
    const int K = std::max(m,n);
    const int rWorkSize = k*std::max(5*k+7,2*K+2*k+1);
    std::vector<float> rWork(rWorkSize);
    std::vector<int> iWork(8*k);

    scomplex workDummy;
    ELEM_LAPACK(cgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, &workDummy, &workSize,
      rWork.data(), iWork.data(), &info );

    workSize = workDummy.real();
    std::vector<scomplex> work(workSize);
    ELEM_LAPACK(cgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, work.data(), &workSize,
      rWork.data(), iWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesdd's updating process failed");
}

void DivideAndConquerSVD
( int m, int n, dcomplex* A, int ldA, 
  double* s, dcomplex* U, int ldu, dcomplex* VH, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int workSize=-1, info;
    dcomplex workDummy;
    const int k = std::min(m,n);
    const int K = std::max(m,n);
    const int rWorkSize = k*std::max(5*k+7,2*K+2*k+1);
    std::vector<double> rWork(rWorkSize);
    std::vector<int> iWork(8*k);

    ELEM_LAPACK(zgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, &workDummy, &workSize,
      rWork.data(), iWork.data(), &info );

    workSize = workDummy.real();
    std::vector<dcomplex> work(workSize);
    ELEM_LAPACK(zgesdd)
    ( &jobz, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, work.data(), &workSize,
      rWork.data(), iWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesdd's updating process failed");
}

// QR-algorithm SVD
// ================

void QRSVD
( int m, int n, float* A, int ldA, 
  float* s, float* U, int ldu, float* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='S', jobVT='S';
    int workSize=-1, info;
    float workDummy;

    ELEM_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<float> work(workSize);
    ELEM_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void QRSVD
( int m, int n, double* A, int ldA, 
  double* s, double* U, int ldu, double* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='S', jobVT='S';
    int workSize=-1, info;
    double workDummy;

    ELEM_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<double> work(workSize);
    ELEM_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, U, &ldu, VTrans, &ldvt, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void QRSVD
( int m, int n, scomplex* A, int ldA, 
  float* s, scomplex* U, int ldu, scomplex* VH, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='S', jobVH='S';
    int workSize=-1, info;
    const int k = std::min(m,n);
    std::vector<float> rWork(5*k);

    scomplex workDummy;
    ELEM_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    std::vector<scomplex> work(workSize);
    ELEM_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void QRSVD
( int m, int n, dcomplex* A, int ldA, 
  double* s, dcomplex* U, int ldu, dcomplex* VH, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='S', jobVH='S';
    int workSize=-1, info;
    dcomplex workDummy;
    const int k = std::min(m,n);
    std::vector<double> rWork(5*k);

    ELEM_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    std::vector<dcomplex> work(workSize);
    ELEM_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, U, &ldu, VH, &ldva, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Compute singular values (with DQDS)
// ===================================

void SVD( int m, int n, float* A, int ldA, float* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVT='N';
    int fakeLDim=1, workSize=-1, info;
    float workDummy;

    ELEM_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<float> work(workSize);
    ELEM_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void SVD( int m, int n, double* A, int ldA, double* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVT='N';
    int fakeLDim=1, workSize=-1, info;
    double workDummy;

    ELEM_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<double> work(workSize);
    ELEM_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void SVD( int m, int n, scomplex* A, int ldA, float* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVH='N';
    int fakeLDim=1, workSize=-1, info;
    scomplex workDummy;
    const int k = std::min(m,n);
    std::vector<float> rWork(5*k);

    ELEM_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    std::vector<scomplex> work(workSize);
    ELEM_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void SVD( int m, int n, dcomplex* A, int ldA, double* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVH='N';
    int fakeLDim=1, workSize=-1, info;
    dcomplex workDummy;
    const int k = std::min(m,n);
    std::vector<double> rWork(5*k);

    ELEM_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    std::vector<dcomplex> work(workSize);
    ELEM_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Compute the Schur decomposition of an upper Hessenberg matrix
// =============================================================

void HessenbergSchur( int n, float* H, int ldH, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, workSize=-1, info;
    float workDummy;
    std::vector<float> wr( n ), wi( n );
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<float> work(workSize);
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
}

void HessenbergSchur
( int n, double* H, int ldH, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, workSize=-1, info;
    double workDummy;
    std::vector<double> wr( n ), wi( n );
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<double> work(workSize);
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");
    
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
}

void HessenbergSchur
( int n, scomplex* H, int ldH, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, workSize=-1, info;
    scomplex workDummy;
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    std::vector<scomplex> work(workSize);
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void HessenbergSchur
( int n, dcomplex* H, int ldH, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, workSize=-1, info;
    dcomplex workDummy;
    ELEM_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    std::vector<dcomplex> work(workSize);
    ELEM_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

void HessenbergSchur
( int n, float* H, int ldH, scomplex* w, float* Q, int ldQ, 
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    int ilo=1, ihi=n;
    int workSize=-1, info;
    float workDummy;
    std::vector<float> wr( n ), wi( n );
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<float> work(workSize);
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
}

void HessenbergSchur
( int n, double* H, int ldH, dcomplex* w, double* Q, int ldQ, 
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    int ilo=1, ihi=n;
    int workSize=-1, info;
    double workDummy;
    std::vector<double> wr( n ), wi( n );
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<double> work(workSize);
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");
    
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
}

void HessenbergSchur
( int n, scomplex* H, int ldH, scomplex* w, scomplex* Q, int ldQ,
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    int ilo=1, ihi=n;
    int workSize=-1, info;
    scomplex workDummy;
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    std::vector<scomplex> work(workSize);
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void HessenbergSchur
( int n, dcomplex* H, int ldH, dcomplex* w, dcomplex* Q, int ldQ,
  bool fullTriangle, bool multiplyQ )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;

    const char job=(fullTriangle?'S':'E'), compZ=(multiplyQ?'V':'I');
    int ilo=1, ihi=n;
    int workSize=-1, info;
    dcomplex workDummy;
    ELEM_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      &workDummy, &workSize, &info );

    workSize = workDummy.real();
    std::vector<dcomplex> work(workSize);
    ELEM_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

// Compute eigenvalues/pairs of an upper Hessenberg matrix
// =======================================================

void HessenbergEig( int n, float* H, int ldH, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( int n, double* H, int ldH, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( int n, scomplex* H, int ldH, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( int n, dcomplex* H, int ldH, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    HessenbergSchur( n, H, ldH, w, false );
}

// TODO: Compute eigenpairs

// Compute the Schur decomposition of a square matrix
// ==================================================

void Schur( int n, float* A, int ldA, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    std::vector<float> tau( n );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    int fakeLDim=1, negOne=-1;
    std::vector<float> wr( n ), wi( n );
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy), workSize );

    // Reduce to Hessenberg form
    std::vector<float> work( workSize );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
}

void Schur( int n, double* A, int ldA, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    std::vector<double> tau( n );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    int fakeLDim=1, negOne=-1;
    std::vector<double> wr( n ), wi( n );
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy), workSize );

    // Reduce to Hessenberg form
    std::vector<double> work( workSize );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
}

void Schur( int n, scomplex* A, int ldA, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    std::vector<scomplex> tau( n );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    int fakeLDim=1, negOne=-1;
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    std::vector<scomplex> work( workSize );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Schur( int n, dcomplex* A, int ldA, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    std::vector<dcomplex> tau( n );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    int fakeLDim=1, negOne=-1;
    ELEM_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    std::vector<dcomplex> work( workSize );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

void Schur
( int n, float* A, int ldA, scomplex* w, float* Q, int ldQ, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    std::vector<float> tau( n );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    std::vector<float> wr( n ), wi( n );
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy), workSize );

    // Reduce to Hessenberg form
    std::vector<float> work( workSize );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
}

void Schur
( int n, double* A, int ldA, dcomplex* w, double* Q, int ldQ, 
  bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    std::vector<double> tau( n );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    std::vector<double> wr( n ), wi( n );
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy), workSize );

    // Reduce to Hessenberg form
    std::vector<double> work( workSize );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
}

void Schur
( int n, scomplex* A, int ldA, scomplex* w, scomplex* Q, int ldQ, 
  bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    std::vector<scomplex> tau( n );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy.real()), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    std::vector<scomplex> work( workSize );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, work.data(), &workSize, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Schur
( int n, dcomplex* A, int ldA, dcomplex* w, dcomplex* Q, int ldQ, 
  bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    std::vector<dcomplex> tau( n );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy.real()), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    ELEM_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = std::max( int(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    std::vector<dcomplex> work( workSize );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(zhseqr)
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

void Eig( int n, float* A, int ldA, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

void Eig( int n, double* A, int ldA, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

void Eig( int n, scomplex* A, int ldA, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

void Eig( int n, dcomplex* A, int ldA, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    Schur( n, A, ldA, w, false );
}

// Eigenpairs
// ----------
// NOTE: When the matrices are real, an interface is also provided which returns
//       a packing of the eigenvectors which exploits the fact that, if the
//       eigenvalue is real, so is the corresponding eigenvector, otherwise
//       the eigenvalue's complex conjugate is also an eigenvalue, and the 
//       corresponding eigenvector is also the conjugate. Thus, an n x n
//       real matrix can be used to represent the eigenvectors if
//           x(j  ) = X(:,j) + X(:,j+1)*1i,
//           x(j+1) = X(:,j) - X(:,j+1)*1i
//       when the j'th and j+1'th eigenvalues are complex conjugates.

void Eig( int n, float* A, int ldA, scomplex* w, float* XPacked, int ldX )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    const char jobVL='N', jobVR='V';
    const int fakeLDim = 1;

    std::vector<float> wReal(n), wImag(n);
    int workSize=-1, info;
    float workDummy;
    ELEM_LAPACK(sgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim, 
      XPacked, &ldX, &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<float> work( workSize );
    ELEM_LAPACK(sgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim, 
      XPacked, &ldX, work.data(), &workSize, &info );

    // Post-process the eigenvalues
    for( Int j=0; j<n; ++j )
        w[j] = Complex<float>(wReal[j],wImag[j]);
}

void Eig( int n, double* A, int ldA, dcomplex* w, double* XPacked, int ldX )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    const char jobVL='N', jobVR='V';
    const int fakeLDim = 1;

    std::vector<double> wReal(n), wImag(n);
    int workSize=-1, info;
    double workDummy;
    ELEM_LAPACK(dgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim,
      XPacked, &ldX, &workDummy, &workSize, &info );

    workSize = workDummy;
    std::vector<double> work( workSize );
    ELEM_LAPACK(dgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, wReal.data(), wImag.data(), 0, &fakeLDim,
      XPacked, &ldX, work.data(), &workSize, &info );

    // Post-process the eigenvalues
    for( Int j=0; j<n; ++j )
        w[j] = Complex<double>(wReal[j],wImag[j]);
}

void Eig( int n, float* A, int ldA, scomplex* w, scomplex* X, int ldX )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    float* XPacked = (float*)X;    
    Eig( n, A, ldA, w, XPacked, ldX );
    // Unpack the eigenvectors
    std::vector<scomplex> z(n);
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

void Eig( int n, double* A, int ldA, dcomplex* w, dcomplex* X, int ldX )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    double* XPacked = (double*)X;    
    Eig( n, A, ldA, w, XPacked, ldX );
    // Unpack the eigenvectors
    std::vector<scomplex> z(n);
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

void Eig( int n, scomplex* A, int ldA, scomplex* w, scomplex* X, int ldX )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    std::vector<float> rWork( 2*n );
    const char jobVL='N', jobVR='V';
    const int fakeLDim = 1;

    int workSize=-1, info;
    scomplex workDummy;
    ELEM_LAPACK(cgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    std::vector<scomplex> work( workSize );
    ELEM_LAPACK(cgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX, 
      work.data(), &workSize, rWork.data(), &info );
}

void Eig( int n, dcomplex* A, int ldA, dcomplex* w, dcomplex* X, int ldX )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    std::vector<double> rWork( 2*n );
    const char jobVL='N', jobVR='V';
    const int fakeLDim = 1;

    int workSize=-1, info;
    dcomplex workDummy;
    ELEM_LAPACK(zgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX,
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    std::vector<dcomplex> work( workSize );
    ELEM_LAPACK(zgeev)
    ( &jobVL, &jobVR, &n, A, &ldA, w, 0, &fakeLDim, X, &ldX,
      work.data(), &workSize, rWork.data(), &info );
}

// TODO: Return the left eigenvectors?

} // namespace lapack
} // namespace elem
