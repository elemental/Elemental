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
  const int* il, const int* iu, const float* abstol, int * m, 
  float* w, float* Z, const int* ldz, int* isuppz, 
  float* work, const int* lwork, int* iwork, const int* liwork, int* info );
void ELEM_LAPACK(dstevr)
( const char* job, const char* range, const int* n,
  double* d, double* e, const double* vl, const double* vu, 
  const int* il, const int* iu, const double* abstol, int * m, 
  double* w, double* Z, const int* ldz, int* isuppz, 
  double* work, const int* lwork, int* iwork, const int* liwork, int* info );

// Hermitian eigensolvers (via MRRR)
void ELEM_LAPACK(ssyevr)
( const char* job, const char* range, const char* uplo, const int* n,
  float* A, const int* lda, const float* vl, const float* vu, 
  const int* il, const int* iu, const float* abstol, int * m, 
  float* w, float* Z, const int* ldz, int* isuppz, 
  float* work, const int* lwork, int* iwork, const int* liwork, int* info );
void ELEM_LAPACK(dsyevr)
( const char* job, const char* range, const char* uplo, const int* n,
  double* A, const int* lda, const double* vl, const double* vu, 
  const int* il, const int* iu, const double* abstol, int * m, 
  double* w, double* Z, const int* ldz, int* isuppz, 
  double* work, const int* lwork, int* iwork, const int* liwork, int* info );
void ELEM_LAPACK(cheevr)
( const char* job, const char* range, const char* uplo, const int* n,
  scomplex* A, const int* lda, const float* vl, const float* vu, 
  const int* il, const int* iu, const float* abstol, int* m,
  float* w, scomplex* Z, const int* ldz, int* isuppz, 
  scomplex* work, const int* lwork, float* rwork, const int* lrwork, 
  int* iwork, const int* liwork, int* info );
void ELEM_LAPACK(zheevr)
( const char* job, const char* range, const char* uplo, const int* n,
  dcomplex* A, const int* lda, const double* vl, const double* vu, 
  const int* il, const int* iu, const double* abstol, int* m,
  double* w, dcomplex* Z, const int* ldz, int* isuppz, 
  dcomplex* work, const int* lwork, double* rwork, const int* lrwork, 
  int* iwork, const int* liwork, int* info );

// Bidiagonal DQDS
void ELEM_LAPACK(slasq1)
( const int* n, float* d, float* e, float* work, int* info );
void ELEM_LAPACK(dlasq1)
( const int* n, double* d, double* e, double* work, int* info );

// Bidiagonal QR
void ELEM_LAPACK(sbdsqr)
( const char* uplo, const int* n, const int* numColsVTrans, const int* numRowsU,
  const int* numColsC, float* d, float* e, float* VTrans, const int* ldVTrans,
  float* U, const int* ldU, float* C, const int* ldC, float* work, int* info );
void ELEM_LAPACK(dbdsqr)
( const char* uplo, const int* n, const int* numColsVTrans, const int* numRowsU,
  const int* numColsC, double* d, double* e,
  double* VTrans, const int* ldVTrans, double* U, const int* ldU,
  double* C, const int* ldC, double* work, int* info );
void ELEM_LAPACK(cbdsqr)
( const char* uplo, const int* n, const int* numColsVAdj, const int* numRowsU,
  const int* numColsC, float* d, float* e,
  scomplex* VAdj, const int* ldVAdj, scomplex* U, const int* ldU,
  scomplex* C, const int* ldC, float* work, int* info );
void ELEM_LAPACK(zbdsqr)
( const char* uplo, const int* n, const int* numColsVAdj, const int* numRowsU,
  const int* numColsC, double* d, double* e,
  dcomplex* VAdj, const int* ldVAdj, dcomplex* U, const int* ldU,
  dcomplex* C, const int* ldC, double* work, int* info );

// Divide and Conquer SVD
void ELEM_LAPACK(sgesdd)
( const char* jobz, const int* m, const int* n, float* A, const int* lda,
  float* s, float* U, const int* ldu, float* VTrans, const int* ldvt,
  float* work, const int* lwork, int* iwork, int* info );
void ELEM_LAPACK(dgesdd)
( const char* jobz, const int* m, const int* n, double* A, const int* lda,
  double* s, double* U, const int* ldu, double* VTrans, const int* ldvt,
  double* work, const int* lwork, int* iwork, int* info );
void ELEM_LAPACK(cgesdd)
( const char* jobz, const int* m, const int* n,
  scomplex* A, const int* lda, float* s,
  scomplex* U, const int* ldu, scomplex* VTrans, const int* ldvt,
  scomplex* work, const int* lwork, float* rwork,
  int* iwork, int* info );
void ELEM_LAPACK(zgesdd)
( const char* jobz, const int* m, const int* n,
  dcomplex* A, const int* lda, double* s,
  dcomplex* U, const int* ldu, dcomplex* VAdj, const int* ldva,
  dcomplex* work, const int* lwork, double* rwork,
  int* iwork, int* info );

// QR-algorithm SVD [DQDS when no singular vectors desired]
void ELEM_LAPACK(sgesvd)
( const char* jobu, const char* jobvt, const int* m, const int* n,
  float* A, const int* lda,
  float* s, float* U, const int* ldu, float* VTrans, const int* ldvt,
  float* work, const int* lwork, int* info );
void ELEM_LAPACK(dgesvd)
( const char* jobu, const char* jobvt, const int* m, const int* n,
  double* A, const int* lda,
  double* s, double* U, const int* ldu, double* VTrans, const int* ldvt,
  double* work, const int* lwork, int* info );
void ELEM_LAPACK(cgesvd)
( const char* jobu, const char* jobva, const int* m, const int* n,
  scomplex* A, const int* lda, float* s,
  scomplex* U, const int* ldu, scomplex* VTrans, const int* ldvt,
  scomplex* work, const int* lwork, float* rwork, int* info );
void ELEM_LAPACK(zgesvd)
( const char* jobu, const char* jobva, const int* m, const int* n,
  dcomplex* A, const int* lda, double* s,
  dcomplex* U, const int* ldu, dcomplex* VAdj, const int* ldva,
  dcomplex* work, const int* lwork, double* rwork, int* info );

// Reduction to Hessenberg form
void ELEM_LAPACK(sgehrd)
( const int* n, const int* ilo, const int* ihi, float* A, const int* lda,
  float* tau, float* work, const int* lwork, int* info );
void ELEM_LAPACK(dgehrd)
( const int* n, const int* ilo, const int* ihi, double* A, const int* lda,
  double* tau, double* work, const int* lwork, int* info );
void ELEM_LAPACK(cgehrd)
( const int* n, const int* ilo, const int* ihi, scomplex* A, const int* lda,
  scomplex* tau, scomplex* work, const int* lwork, int* info );
void ELEM_LAPACK(zgehrd)
( const int* n, const int* ilo, const int* ihi, dcomplex* A, const int* lda,
  dcomplex* tau, dcomplex* work, const int* lwork, int* info );

// Generates a unitary matrix defined as the product of Householder reflectors
void ELEM_LAPACK(sorghr)
( const int* n, const int* ilo, const int* ihi, float* A, const int* lda,
  const float* tau, float* work, const int* lwork, int* info );
void ELEM_LAPACK(dorghr)
( const int* n, const int* ilo, const int* ihi, double* A, const int* lda,
  const double* tau, double* work, const int* lwork, int* info );
void ELEM_LAPACK(cunghr)
( const int* n, const int* ilo, const int* ihi, scomplex* A, const int* lda,
  const scomplex* tau, scomplex* work, const int* lwork, int* info );
void ELEM_LAPACK(zunghr)
( const int* n, const int* ilo, const int* ihi, dcomplex* A, const int* lda,
  const dcomplex* tau, dcomplex* work, const int* lwork, int* info );

// Hessenberg QR algorithm
void ELEM_LAPACK(shseqr)
( const char* job, const char* compz, const int* n, 
  const int* ilo, const int* ihi, float* H, const int* ldh, 
  float* wr, float* wi, float* Z, const int* ldz, 
  float* work, const int* lwork, int* info );
void ELEM_LAPACK(dhseqr)
( const char* job, const char* compz, const int* n, 
  const int* ilo, const int* ihi, double* H, const int* ldh, 
  double* wr, double* wi, double* Z, const int* ldz, 
  double* work, const int* lwork, int* info );
void ELEM_LAPACK(chseqr)
( const char* job, const char* compz, const int* n,
  const int* ilo, const int* ihi, scomplex* H, const int* ldh,
  scomplex* w, scomplex* Z, const int* ldz,
  scomplex* work, const int* lwork, int* info );
void ELEM_LAPACK(zhseqr)
( const char* job, const char* compz, const int* n,
  const int* ilo, const int* ihi, dcomplex* H, const int* ldh,
  dcomplex* w, dcomplex* Z, const int* ldz,
  dcomplex* work, const int* lwork, int* info );

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

int SymmetricTridiagEig
( char job, char range, int n, float* d, float* e, float vl, float vu,
  int il, int iu, float abstol, float* w, float* Z, int ldz )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    if( n == 0 )
        return 0;

    std::vector<int> isuppz( 2*n );

    int lwork=-1, liwork=-1, m, info;
    int dummyIwork;
    float dummyWork;
    ELEM_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), &dummyWork, &lwork, &dummyIwork, &liwork,
      &info );

    lwork = dummyWork;
    liwork = dummyIwork;
    std::vector<float> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sstevr's failed");
    return m;
}

int SymmetricTridiagEig
( char job, char range, int n, double* d, double* e, double vl, double vu,
  int il, int iu, double abstol, double* w, double* Z, int ldz )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SymmetricTridiagEig"));
    if( n == 0 )
        return 0;

    std::vector<int> isuppz( 2*n );

    int lwork=-1, liwork=-1, m, info;
    int dummyIwork;
    double dummyWork;
    ELEM_LAPACK(dstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), &dummyWork, &lwork, &dummyIwork, &liwork,
      &info );

    lwork = dummyWork;
    liwork = dummyIwork;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_LAPACK(dstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dstevr's failed");
    return m;
}

// Compute the EVD of a Hermitian matrix
// =====================================

int HermitianEig
( char job, char range, char uplo, int n, float* A, int lda, 
  float vl, float vu, int il, int iu, float abstol, 
  float* w, float* Z, int ldz )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppz( 2*n );

    int lwork=-1, liwork=-1, m, info;
    int dummyIwork;
    float dummyWork;
    ELEM_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), &dummyWork, &lwork, &dummyIwork, &liwork, 
      &info );

    lwork = dummyWork;
    liwork = dummyIwork;
    std::vector<float> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("ssyevr's failed");
    return m;
}

int HermitianEig
( char job, char range, char uplo, int n, double* A, int lda, 
  double vl, double vu, int il, int iu, double abstol, 
  double* w, double* Z, int ldz )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppz( 2*n );

    int lwork=-1, liwork=-1, m, info;
    int dummyIwork;
    double dummyWork;
    ELEM_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), &dummyWork, &lwork, &dummyIwork, &liwork, 
      &info );

    lwork = dummyWork;
    liwork = dummyIwork;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    ELEM_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dsyevr's failed");
    return m;
}

int HermitianEig
( char job, char range, char uplo, int n, scomplex* A, int lda, 
  float vl, float vu, int il, int iu, float abstol, 
  float* w, scomplex* Z, int ldz )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppz( 2*n );

    int lwork=-1, lrwork=-1, liwork=-1, m, info;
    int dummyIwork;
    float dummyRwork;
    scomplex dummyWork;
    ELEM_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), &dummyWork, &lwork, &dummyRwork, &lrwork,
      &dummyIwork, &liwork, &info );

    lwork = dummyWork.real();
    lrwork = dummyRwork;
    liwork = dummyIwork;
    std::vector<scomplex> work(lwork);
    std::vector<float> rwork(lrwork);
    std::vector<int> iwork(liwork);
    ELEM_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), work.data(), &lwork, rwork.data(), &lrwork,
      iwork.data(), &liwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cheevr's failed");
    return m;
}

int HermitianEig
( char job, char range, char uplo, int n, dcomplex* A, int lda, 
  double vl, double vu, int il, int iu, double abstol, 
  double* w, dcomplex* Z, int ldz )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HermitianEig"))
    if( n == 0 )
        return 0;

    std::vector<int> isuppz( 2*n );

    int lwork=-1, lrwork=-1, liwork=-1, m, info;
    int dummyIwork;
    double dummyRwork;
    dcomplex dummyWork;
    ELEM_LAPACK(zheevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), &dummyWork, &lwork, &dummyRwork, &lrwork,
      &dummyIwork, &liwork, &info );

    lwork = dummyWork.real();
    lrwork = dummyRwork;
    liwork = dummyIwork;
    std::vector<dcomplex> work(lwork);
    std::vector<double> rwork(lrwork);
    std::vector<int> iwork(liwork);
    ELEM_LAPACK(zheevr)
    ( &job, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu, &abstol, &m,
      w, Z, &ldz, isuppz.data(), work.data(), &lwork, rwork.data(), &lrwork,
      iwork.data(), &liwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zheevr's failed");
    return m;
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
( char uplo, int n, int numColsVTrans, int numRowsU,
  float* d, float* e, float* VTrans, int ldVTrans, float* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    float* C=0;
    const int numColsC=0, ldC=1;
    std::vector<float> work( 4*n );
    ELEM_LAPACK(sbdsqr)
    ( &uplo, &n, &numColsVTrans, &numRowsU, &numColsC, d, e, VTrans, &ldVTrans,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, int n, int numColsVTrans, int numRowsU, 
  double* d, double* e, double* VTrans, int ldVTrans, double* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    double* C=0;
    const int numColsC=0, ldC=1;
    std::vector<double> work( 4*n );
    ELEM_LAPACK(dbdsqr)
    ( &uplo, &n, &numColsVTrans, &numRowsU, &numColsC, d, e, VTrans, &ldVTrans,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, int n, int numColsVAdj, int numRowsU, 
  float* d, float* e, scomplex* VAdj, int ldVAdj, scomplex* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    scomplex* C=0;
    const int numColsC=0, ldC=1;
    std::vector<float> work( 4*n );
    ELEM_LAPACK(cbdsqr)
    ( &uplo, &n, &numColsVAdj, &numRowsU, &numColsC, d, e, VAdj, &ldVAdj,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cbdsqr had ",info," elements of e not converge");
}

void BidiagQRAlg
( char uplo, int n, int numColsVAdj, int numRowsU, 
  double* d, double* e, dcomplex* VAdj, int ldVAdj, dcomplex* U, int ldU )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::BidiagQRAlg"))
    if( n==0 )
        return;

    int info;
    dcomplex* C=0;
    const int numColsC=0, ldC=1;
    std::vector<double> work( 4*n );
    ELEM_LAPACK(zbdsqr)
    ( &uplo, &n, &numColsVAdj, &numRowsU, &numColsC, d, e, VAdj, &ldVAdj,
      U, &ldU, C, &ldC, work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zbdsqr had ",info," elements of e not converge");
}

// Divide and Conquer SVD
// ======================

void DivideAndConquerSVD
( int m, int n, float* A, int lda, 
  float* s, float* U, int ldu, float* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int lwork=-1, info;
    float dummyWork;
    const int k = std::min(m,n);
    std::vector<int> iwork(8*k);

    ELEM_LAPACK(sgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, &dummyWork, &lwork,
      iwork.data(), &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    ELEM_LAPACK(sgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, work.data(), &lwork,
      iwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesdd's updating process failed");
}

void DivideAndConquerSVD
( int m, int n, double* A, int lda, 
  double* s, double* U, int ldu, double* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int lwork=-1, info;
    double dummyWork;
    const int k = std::min(m,n);
    std::vector<int> iwork(8*k);

    ELEM_LAPACK(dgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, &dummyWork, &lwork,
      iwork.data(), &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    ELEM_LAPACK(dgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, work.data(), &lwork,
      iwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesdd's updating process failed");
}

void DivideAndConquerSVD
( int m, int n, scomplex* A, int lda, 
  float* s, scomplex* U, int ldu, scomplex* VAdj, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int lwork=-1, info;
    const int k = std::min(m,n);
    const int K = std::max(m,n);
    const int lrwork = k*std::max(5*k+7,2*K+2*k+1);
    std::vector<float> rwork(lrwork);
    std::vector<int> iwork(8*k);

    scomplex dummyWork;
    ELEM_LAPACK(cgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, &dummyWork, &lwork,
      rwork.data(), iwork.data(), &info );

    lwork = dummyWork.real();
    std::vector<scomplex> work(lwork);
    ELEM_LAPACK(cgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, work.data(), &lwork,
      rwork.data(), iwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesdd's updating process failed");
}

void DivideAndConquerSVD
( int m, int n, dcomplex* A, int lda, 
  double* s, dcomplex* U, int ldu, dcomplex* VAdj, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::DivideAndConquerSVD"))
    if( m==0 || n==0 )
        return;

    const char jobz='S';
    int lwork=-1, info;
    dcomplex dummyWork;
    const int k = std::min(m,n);
    const int K = std::max(m,n);
    const int lrwork = k*std::max(5*k+7,2*K+2*k+1);
    std::vector<double> rwork(lrwork);
    std::vector<int> iwork(8*k);

    ELEM_LAPACK(zgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, &dummyWork, &lwork,
      rwork.data(), iwork.data(), &info );

    lwork = dummyWork.real();
    std::vector<dcomplex> work(lwork);
    ELEM_LAPACK(zgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, work.data(), &lwork,
      rwork.data(), iwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesdd's updating process failed");
}

// QR-algorithm SVD
// ================

void QRSVD
( int m, int n, float* A, int lda, 
  float* s, float* U, int ldu, float* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='S', jobvt='S';
    int lwork=-1, info;
    float dummyWork;

    ELEM_LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    ELEM_LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void QRSVD
( int m, int n, double* A, int lda, 
  double* s, double* U, int ldu, double* VTrans, int ldvt )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='S', jobvt='S';
    int lwork=-1, info;
    double dummyWork;

    ELEM_LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    ELEM_LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void QRSVD
( int m, int n, scomplex* A, int lda, 
  float* s, scomplex* U, int ldu, scomplex* VAdj, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='S', jobva='S';
    int lwork=-1, info;
    const int k = std::min(m,n);
    std::vector<float> rwork(5*k);

    scomplex dummyWork;
    ELEM_LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      &dummyWork, &lwork, rwork.data(), &info );

    lwork = dummyWork.real();
    std::vector<scomplex> work(lwork);
    ELEM_LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      work.data(), &lwork, rwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void QRSVD
( int m, int n, dcomplex* A, int lda, 
  double* s, dcomplex* U, int ldu, dcomplex* VAdj, int ldva )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::QRSVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='S', jobva='S';
    int lwork=-1, info;
    dcomplex dummyWork;
    const int k = std::min(m,n);
    std::vector<double> rwork(5*k);

    ELEM_LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      &dummyWork, &lwork, rwork.data(), &info );

    lwork = dummyWork.real();
    std::vector<dcomplex> work(lwork);
    ELEM_LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      work.data(), &lwork, rwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Compute singular values (with DQDS)
// ===================================

void SVD( int m, int n, float* A, int lda, float* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='N', jobvt='N';
    int fakeLDim=1, lwork=-1, info;
    float dummyWork;

    ELEM_LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    ELEM_LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void SVD( int m, int n, double* A, int lda, double* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='N', jobvt='N';
    int fakeLDim=1, lwork=-1, info;
    double dummyWork;

    ELEM_LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    ELEM_LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void SVD( int m, int n, scomplex* A, int lda, float* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='N', jobva='N';
    int fakeLDim=1, lwork=-1, info;
    scomplex dummyWork;
    const int k = std::min(m,n);
    std::vector<float> rwork(5*k);

    ELEM_LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, rwork.data(), &info );

    lwork = dummyWork.real();
    std::vector<scomplex> work(lwork);
    ELEM_LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &lwork, rwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void SVD( int m, int n, dcomplex* A, int lda, double* s )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::SVD"))
    if( m==0 || n==0 )
        return;

    const char jobu='N', jobva='N';
    int fakeLDim=1, lwork=-1, info;
    dcomplex dummyWork;
    const int k = std::min(m,n);
    std::vector<double> rwork(5*k);

    ELEM_LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, rwork.data(), &info );

    lwork = dummyWork.real();
    std::vector<dcomplex> work(lwork);
    ELEM_LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &lwork, rwork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Compute the Schur decomposition of an upper Hessenberg matrix
// =============================================================

void HessenbergEig( int n, float* H, int ldh, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    if( n == 0 )
        return;

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    float dummyWork;
    std::vector<float> wr( n ), wi( n );
    ELEM_LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, wr.data(), wi.data(), 0, &fakeLDim,
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    ELEM_LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
}

void HessenbergEig( int n, double* H, int ldh, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    if( n == 0 )
        return;

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    double dummyWork;
    std::vector<double> wr( n ), wi( n );
    ELEM_LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, wr.data(), wi.data(), 0, &fakeLDim,
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    ELEM_LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");
    
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
}

void HessenbergEig( int n, scomplex* H, int ldh, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    if( n == 0 )
        return;

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    scomplex dummyWork;
    ELEM_LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork.real();
    std::vector<scomplex> work(lwork);
    ELEM_LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void HessenbergEig( int n, dcomplex* H, int ldh, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergEig"))
    if( n == 0 )
        return;

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    dcomplex dummyWork;
    ELEM_LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork.real();
    std::vector<dcomplex> work(lwork);
    ELEM_LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

void HessenbergSchur( int n, float* H, int ldh, float* Q, int ldq, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;
    LogicError("This routine not yet written");
}

void HessenbergSchur
( int n, double* H, int ldh, double* Q, int ldq, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;
    LogicError("This routine not yet written");
}

void HessenbergSchur
( int n, scomplex* H, int ldh, scomplex* Q, int ldq, scomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;
    LogicError("This routine not yet written");
}

void HessenbergSchur
( int n, dcomplex* H, int ldh, dcomplex* Q, int ldq, dcomplex* w )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::HessenbergSchur"))
    if( n == 0 )
        return;
    LogicError("This routine not yet written");
}

// Compute the Schur decomposition of a square matrix
// ==================================================

void Eig( int n, float* A, int lda, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    float dummyWork;
    std::vector<float> tau( n );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='N';
    int fakeLDim=1, negOne=-1;
    std::vector<float> wr( n ), wi( n );
    ELEM_LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), 0, &fakeLDim,
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork), lwork );

    // Reduce to Hessenberg form
    std::vector<float> work( lwork );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
}

void Eig( int n, double* A, int lda, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    double dummyWork;
    std::vector<double> tau( n );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='N';
    int fakeLDim=1, negOne=-1;
    std::vector<double> wr( n ), wi( n );
    ELEM_LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), 0, &fakeLDim,
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork), lwork );

    // Reduce to Hessenberg form
    std::vector<double> work( lwork );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
}

void Eig( int n, scomplex* A, int lda, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    scomplex dummyWork;
    std::vector<scomplex> tau( n );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='N';
    int fakeLDim=1, negOne=-1;
    ELEM_LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, 0, &fakeLDim, 
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork.real()), lwork );

    // Reduce to Hessenberg form
    std::vector<scomplex> work( lwork );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, 0, &fakeLDim, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Eig( int n, dcomplex* A, int lda, dcomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Eig"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    dcomplex dummyWork;
    std::vector<dcomplex> tau( n );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='N';
    int fakeLDim=1, negOne=-1;
    ELEM_LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, 0, &fakeLDim, 
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork.real()), lwork );

    // Reduce to Hessenberg form
    std::vector<dcomplex> work( lwork );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    ELEM_LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, 0, &fakeLDim, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

void Schur
( int n, float* A, int lda, float* Q, int ldq, scomplex* w, bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    float dummyWork;
    std::vector<float> tau( n );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork;

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork), lwork );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='V';
    std::vector<float> wr( n ), wi( n );
    ELEM_LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), Q, &ldq, 
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork), lwork );

    // Reduce to Hessenberg form
    std::vector<float> work( lwork );
    ELEM_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldq], &A[j*lda], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), Q, &ldq, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
}

void Schur
( int n, double* A, int lda, double* Q, int ldq, dcomplex* w, 
  bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    double dummyWork;
    std::vector<double> tau( n );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork;

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork), lwork );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='V';
    std::vector<double> wr( n ), wi( n );
    ELEM_LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), Q, &ldq, 
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork), lwork );

    // Reduce to Hessenberg form
    std::vector<double> work( lwork );
    ELEM_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldq], &A[j*lda], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, wr.data(), wi.data(), Q, &ldq, 
      work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
}

void Schur
( int n, scomplex* A, int lda, scomplex* Q, int ldq, scomplex* w, 
  bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    scomplex dummyWork;
    std::vector<scomplex> tau( n );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork.real();

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork.real()), lwork );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='V';
    ELEM_LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, Q, &ldq, 
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork.real()), lwork );

    // Reduce to Hessenberg form
    std::vector<scomplex> work( lwork );
    ELEM_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldq], &A[j*lda], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, Q, &ldq, work.data(), &lwork, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Schur
( int n, dcomplex* A, int lda, dcomplex* Q, int ldq, dcomplex* w, 
  bool fullTriangle )
{
    DEBUG_ONLY(CallStackEntry cse("lapack::Schur"))
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    int ilo=1, ihi=n, lwork=-1, info;
    dcomplex dummyWork;
    std::vector<dcomplex> tau( n );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), &dummyWork, &lwork, &info );
    lwork = dummyWork.real();

    // Query the explicit Q formation workspace
    int negOne=-1; 
    ELEM_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork.real()), lwork );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compz='V';
    ELEM_LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, Q, &ldq, 
      &dummyWork, &negOne, &info );
    lwork = std::max( int(dummyWork.real()), lwork );

    // Reduce to Hessenberg form
    std::vector<dcomplex> work( lwork );
    ELEM_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &lda, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( int j=0; j<n; ++j )
        MemCopy( &Q[j*ldq], &A[j*lda], n );

    // Form the orthogonal matrix in place
    ELEM_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldq, tau.data(), work.data(), &lwork, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    ELEM_LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, A, &lda, w, Q, &ldq, work.data(), &lwork, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

} // namespace lapack
} // namespace elem
