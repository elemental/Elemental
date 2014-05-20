/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

using El::scomplex;
using El::dcomplex;

extern "C" {

// Level 1 BLAS
// ============
void EL_BLAS(saxpy)
( const int* n, const float* alpha, const float* x, const int* incx,
                                          float* y, const int* incy );
void EL_BLAS(daxpy)
( const int* n, const double* alpha, const double* x, const int* incx,
                                           double* y, const int* incy );
void EL_BLAS(caxpy)
( const int* n, const scomplex* alpha,
  const scomplex* x, const int* incx,
        scomplex* y, const int* incy );
void EL_BLAS(zaxpy)
( const int* n, const dcomplex* alpha,
  const dcomplex* x, const int* incx,
        dcomplex* y, const int* incy );

void EL_BLAS(scopy)
( const int* n, const float* x, const int* incx,
                      float* y, const int* incy );
void EL_BLAS(dcopy)
( const int* n, const double* x, const int* incx,
                      double* y, const int* incy );
void EL_BLAS(ccopy)
( const int* n, const scomplex* x, const int* incx,
                      scomplex* y, const int* incy );
void EL_BLAS(zcopy)
( const int* n, const dcomplex* x, const int* incx,
                      dcomplex* y, const int* incy );

float EL_BLAS(sdot)
( const int* n, const float* x, const int* incx,
                const float* y, const int* incy );
double EL_BLAS(ddot)
( const int* n, const double* x, const int* incx,
                const double* y, const int* incy );
// To avoid the compatibility issue, we simply handroll our own complex dots
float  EL_BLAS(snrm2) ( const int* n, const float   * x, const int* incx );
double EL_BLAS(dnrm2) ( const int* n, const double  * x, const int* incx );
float  EL_BLAS(scnrm2)( const int* n, const scomplex* x, const int* incx );
double EL_BLAS(dznrm2)( const int* n, const dcomplex* x, const int* incx );

// Apply a Givens rotation to a pair of vectors
void EL_BLAS(srot)
( const int* n, float* x, const int* incx, float* y, const int* incy,
  const float* c, const float* s );
void EL_BLAS(drot)
( const int* n, double* x, const int* incx, double* y, const int* incy,
  const double* c, const double* s );
void EL_BLAS(crot)
( const int* n, scomplex* x, const int* incx, scomplex* y, const int* incy,
  const float* c, const scomplex* s );
void EL_BLAS(zrot)
( const int* n, dcomplex* x, const int* incx, dcomplex* y, const int* incy,
  const double* c, const dcomplex* s );

// Quickly compute a Givens rotation
void EL_BLAS(srotg)
( float* alpha, float* beta, float* c, float* s );
void EL_BLAS(drotg)
( double* alpha, double* beta, double* c, double* s );
void EL_BLAS(crotg)
( scomplex* alpha, scomplex* beta, float* c, scomplex* s );
void EL_BLAS(zrotg)
( dcomplex* alpha, dcomplex* beta, double* c, dcomplex* s );

// Scale a vector
void EL_BLAS(sscal)
( const int* n, const float   * alpha, float   * x, const int* incx );
void EL_BLAS(dscal)
( const int* n, const double  * alpha, double  * x, const int* incx );
void EL_BLAS(cscal)
( const int* n, const scomplex* alpha, scomplex* x, const int* incx );
void EL_BLAS(zscal)
( const int* n, const dcomplex* alpha, dcomplex* x, const int* incx );

float  EL_BLAS(sasum) ( const int* n, const float   * x, const int* incx );
double EL_BLAS(dasum) ( const int* n, const double  * x, const int* incx );
float  EL_BLAS(scasum)( const int* n, const scomplex* x, const int* incx );
double EL_BLAS(dzasum)( const int* n, const dcomplex* x, const int* incx );
float  EL_LAPACK(scsum1)( const int* n, const scomplex* x, const int* incx );
double EL_LAPACK(dzsum1)( const int* n, const dcomplex* x, const int* incx );

void EL_BLAS(sswap)
( const int* n, float   * x, const int* incx, float   * y, const int* incy );
void EL_BLAS(dswap)
( const int* n, double  * x, const int* incx, double  * y, const int* incy );
void EL_BLAS(cswap)
( const int* n, scomplex* x, const int* incx, scomplex* y, const int* incy );
void EL_BLAS(zswap)
( const int* n, dcomplex* x, const int* incx, dcomplex* y, const int* incy );

// Level 2 BLAS
// ============
void EL_BLAS(sgemv)
( const char* trans, const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                      const float* x, const int* incx,
  const float* beta,        float* y, const int* incy );
void EL_BLAS(dgemv)
( const char* trans, const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                       const double* x, const int* incx,
  const double* beta,        double* y, const int* incy );
void EL_BLAS(cgemv)
( const char* trans, const int* m, const int* n,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* x, const int* incx,
  const scomplex* beta,
        scomplex* y, const int* incy );
void EL_BLAS(zgemv)
( const char* trans, const int* m, const int* n,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* x, const int* incx,
  const dcomplex* beta,
        dcomplex* y, const int* incy );

void EL_BLAS(sger)
( const int* m, const int* n,
  const float* alpha, const float* x, const int* incx,
                      const float* y, const int* incy,
                            float* A, const int* lda  );
void EL_BLAS(dger)
( const int* m, const int* n,
  const double* alpha, const double* x, const int* incx,
                       const double* y, const int* incy,
                             double* A, const int* lda  );
void EL_BLAS(cgerc)
( const int* m, const int* n,
  const scomplex* alpha,
  const scomplex* x, const int* incx,
  const scomplex* y, const int* incy,
        scomplex* A, const int* lda  );
void EL_BLAS(zgerc)
( const int* m, const int* n,
  const dcomplex* alpha,
  const dcomplex* x, const int* incx,
  const dcomplex* y, const int* incy,
        dcomplex* A, const int* lda  );

void EL_BLAS(cgeru)
( const int* m, const int* n,
  const scomplex* alpha,
  const scomplex* x, const int* incx,
  const scomplex* y, const int* incy,
        scomplex* A, const int* lda  );
void EL_BLAS(zgeru)
( const int* m, const int* n,
  const dcomplex* alpha,
  const dcomplex* x, const int* incx,
  const dcomplex* y, const int* incy,
        dcomplex* A, const int* lda  );

void EL_BLAS(chemv)
( const char* uplo, const int* m,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* x, const int* incx,
  const scomplex* beta,
        scomplex* y, const int* incy );
void EL_BLAS(zhemv)
( const char* uplo, const int* m,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* x, const int* incx,
  const dcomplex* beta,
        dcomplex* y, const int* incy );

void EL_BLAS(cher)
( const char* uplo, const int* m,
  const scomplex* alpha,
  const scomplex* x, const int* incx,
        scomplex* A, const int* lda  );
void EL_BLAS(zher)
( const char* uplo, const int* m,
  const dcomplex* alpha,
  const dcomplex* x, const int* incx,
        dcomplex* A, const int* lda  );

void EL_BLAS(cher2)
( const char* uplo, const int* m,
  const scomplex* alpha,
  const scomplex* x, const int* incx,
  const scomplex* y, const int* incy,
        scomplex* A, const int* lda  );
void EL_BLAS(zher2)
( const char* uplo, const int* m,
  const dcomplex* alpha,
  const dcomplex* x, const int* incx,
  const dcomplex* y, const int* incy,
        dcomplex* A, const int* lda  );

void EL_BLAS(ssymv)
( const char* uplo, const int* m,
  const float* alpha, const float* A, const int* lda,
                      const float* x, const int* incx,
  const float* beta,        float* y, const int* incy );
void EL_BLAS(dsymv)
( const char* uplo, const int* m,
  const double* alpha, const double* A, const int* lda,
                       const double* x, const int* incx,
  const double* beta,        double* y, const int* incy );
// 'csymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csymv)
( const char* uplo, const int* m,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* x, const int* incx,
  const scomplex* beta,
        scomplex* y, const int* incy );
// 'zsymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsymv)
( const char* uplo, const int* m,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* x, const int* incx,
  const dcomplex* beta,
        dcomplex* y, const int* incy );

void EL_BLAS(ssyr)
( const char* uplo, const int* m,
  const float* alpha, const float* x, const int* incx,
                            float* A, const int* lda  );
void EL_BLAS(dsyr)
( const char* uplo, const int* m,
  const double* alpha, const double* x, const int* incx,
                             double* A, const int* lda  );
// 'csyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csyr)
( const char* uplo, const int* m,
  const scomplex* alpha,
  const scomplex* x, const int* incx,
        scomplex* A, const int* lda  );
// 'zsyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsyr)
( const char* uplo, const int* m,
  const dcomplex* alpha,
  const dcomplex* x, const int* incx,
        dcomplex* A, const int* lda  );

void EL_BLAS(ssyr2)
( const char* uplo, const int* m,
  const float* alpha, const float* x, const int* incx,
                      const float* y, const int* incy,
                            float* A, const int* lda  );
void EL_BLAS(dsyr2)
( const char* uplo, const int* m,
  const double* alpha, const double* x, const int* incx,
                       const double* y, const int* incy,
                             double* A, const int* lda  );

void EL_BLAS(strmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const float* A, const int* lda, float* x, const int* incx );
void EL_BLAS(dtrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const double* A, const int* lda, double* x, const int* incx );
void EL_BLAS(ctrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const scomplex* A, const int* lda,
        scomplex* x, const int* incx );
void EL_BLAS(ztrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const dcomplex* A, const int* lda,
        dcomplex* x, const int* incx );

void EL_BLAS(strsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const float* A, const int* lda, float* x, const int* incx );
void EL_BLAS(dtrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const double* A, const int* lda, double* x, const int* incx );
void EL_BLAS(ctrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const scomplex* A, const int* lda,
        scomplex* x, const int* incx );
void EL_BLAS(ztrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const dcomplex* A, const int* lda,
        dcomplex* x, const int* incx );

// Level 3 BLAS 
// ============
void EL_BLAS(sgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );
void EL_BLAS(dgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );
void EL_BLAS(cgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* B, const int* ldb,
  const scomplex* beta,
        scomplex* C, const int* ldc );
void EL_BLAS(zgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* B, const int* ldb,
  const dcomplex* beta,
        dcomplex* C, const int* ldc );

void EL_BLAS(chemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* B, const int* ldb,
  const scomplex* beta,
        scomplex* C, const int* ldc );
void EL_BLAS(zhemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* B, const int* ldb,
  const dcomplex* beta,
        dcomplex* C, const int* ldc );

void EL_BLAS(cher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* B, const int* ldb,
  const scomplex* beta,
        scomplex* C, const int* ldc );
void EL_BLAS(zher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* B, const int* ldb,
  const dcomplex* beta,
        dcomplex* C, const int* ldc );

void EL_BLAS(cherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* beta,
        scomplex* C, const int* ldc );
void EL_BLAS(zherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* beta,
        dcomplex* C, const int* ldc );

void EL_BLAS(ssymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );
void EL_BLAS(dsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );
void EL_BLAS(csymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* B, const int* ldb,
  const scomplex* beta,
        scomplex* C, const int* ldc );
void EL_BLAS(zsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* B, const int* ldb,
  const dcomplex* beta,
        dcomplex* C, const int* ldc );

void EL_BLAS(ssyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );
void EL_BLAS(dsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );
void EL_BLAS(csyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* B, const int* ldb,
  const scomplex* beta,
        scomplex* C, const int* ldc );
void EL_BLAS(zsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* B, const int* ldb,
  const dcomplex* beta,
        dcomplex* C, const int* ldc );

void EL_BLAS(ssyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
  const float* beta,        float* C, const int* ldc );
void EL_BLAS(dsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
  const double* beta,        double* C, const int* ldc );
void EL_BLAS(csyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
  const scomplex* beta,
        scomplex* C, const int* ldc );
void EL_BLAS(zsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
  const dcomplex* beta,
        dcomplex* C, const int* ldc );

void EL_BLAS(strmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                            float* B, const int* ldb );
void EL_BLAS(dtrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                             double* B, const int* ldb );
void EL_BLAS(ctrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
        scomplex* B, const int* ldb );
void EL_BLAS(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
        dcomplex* B, const int* ldb );

void EL_BLAS(strsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                            float* B, const int* ldb );
void EL_BLAS(dtrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                             double* B, const int* ldb );
void EL_BLAS(ctrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const scomplex* alpha,
  const scomplex* A, const int* lda,
        scomplex* B, const int* ldb );
void EL_BLAS(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const dcomplex* alpha,
  const dcomplex* A, const int* lda,
        dcomplex* B, const int* ldb );

} // extern "C"

namespace El {
namespace blas {

// Level 1 BLAS
// ============
void Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy )
{ EL_BLAS(saxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy )
{ EL_BLAS(daxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy )
{ EL_BLAS(caxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy )
{ EL_BLAS(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }

void Copy( int n, const float* x, int incx, float* y, int incy )
{ EL_BLAS(scopy)( &n, x, &incx, y, &incy ); }
void Copy( int n, const double* x, int incx, double* y, int incy )
{ EL_BLAS(dcopy)( &n, x, &incx, y, &incy ); }
void Copy( int n, const scomplex* x, int incx, scomplex* y, int incy )
{ EL_BLAS(ccopy)( &n, x, &incx, y, &incy ); }
void Copy( int n, const dcomplex* x, int incx, dcomplex* y, int incy )
{ EL_BLAS(zcopy)( &n, x, &incx, y, &incy ); }

float Dot( int n, const float* x, int incx, const float* y, int incy )
{ return EL_BLAS(sdot)( &n, x, &incx, y, &incy ); }
double Dot( int n, const double* x, int incx, const double* y, int incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }
scomplex Dot( int n, const scomplex* x, int incx, const scomplex* y, int incy )
{ 
    scomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += Conj(x[i*incx])*y[i*incy];
    return alpha;
}
dcomplex Dot( int n, const dcomplex* x, int incx, const dcomplex* y, int incy )
{
    dcomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += Conj(x[i*incx])*y[i*incy];
    return alpha;
}

float Dotc( int n, const float* x, int incx, const float* y, int incy )
{ return EL_BLAS(sdot)( &n, x, &incx, y, &incy ); }
double Dotc( int n, const double* x, int incx, const double* y, int incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }
scomplex Dotc( int n, const scomplex* x, int incx, const scomplex* y, int incy )
{ 
    scomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += Conj(x[i*incx])*y[i*incy];
    return alpha;
}
dcomplex Dotc( int n, const dcomplex* x, int incx, const dcomplex* y, int incy )
{ 
    dcomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += Conj(x[i*incx])*y[i*incy];
    return alpha;
}

float Dotu( int n, const float* x, int incx, const float* y, int incy )
{ return EL_BLAS(sdot)( &n, x, &incx, y, &incy ); }
double Dotu( int n, const double* x, int incx, const double* y, int incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }
scomplex Dotu( int n, const scomplex* x, int incx, const scomplex* y, int incy )
{
    scomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += x[i*incx]*y[i*incy];
    return alpha;
}
dcomplex Dotu( int n, const dcomplex* x, int incx, const dcomplex* y, int incy )
{
    dcomplex alpha = 0;
    for( int i=0; i<n; ++i ) 
        alpha += x[i*incx]*y[i*incy];
    return alpha;
}

float Nrm2( int n, const float* x, int incx )
{ return EL_BLAS(snrm2)( &n, x, &incx ); }
double Nrm2( int n, const double* x, int incx )
{ return EL_BLAS(dnrm2)( &n, x, &incx ); }
float Nrm2( int n, const scomplex* x, int incx )
{ return EL_BLAS(scnrm2)( &n, x, &incx ); }
double Nrm2( int n, const dcomplex* x, int incx )
{ return EL_BLAS(dznrm2)( &n, x, &incx ); }

float Givens
( float alpha, float beta, float* c, float* s )
{ EL_BLAS(srotg)( &alpha, &beta, c, s ); return alpha; }
double Givens
( double alpha, double beta, double* c, double* s )
{ EL_BLAS(drotg)( &alpha, &beta, c, s ); return alpha; }
scomplex Givens
( scomplex alpha, scomplex beta, float* c, scomplex* s )
{ EL_BLAS(crotg)( &alpha, &beta, c, s ); return alpha; }
dcomplex Givens
( dcomplex alpha, dcomplex beta, double* c, dcomplex* s )
{ EL_BLAS(zrotg)( &alpha, &beta, c, s ); return alpha; }

void Rot
( int n, float* x, int incx, float* y, int incy, float c, float s )
{ EL_BLAS(srot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( int n, double* x, int incx, double* y, int incy, double c, double s )
{ EL_BLAS(drot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( int n, scomplex* x, int incx, scomplex* y, int incy, float c, scomplex s )
{ EL_BLAS(crot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( int n, dcomplex* x, int incx, dcomplex* y, int incy, double c, dcomplex s )
{ EL_BLAS(zrot)( &n, x, &incx, y, &incy, &c, &s ); }

void Scal( int n, float alpha, float* x, int incx )
{ EL_BLAS(sscal)( &n, &alpha, x, &incx ); }
void Scal( int n, double alpha, double* x, int incx )
{ EL_BLAS(dscal)( &n, &alpha, x, &incx ); }
void Scal( int n, scomplex alpha, scomplex* x, int incx )
{ EL_BLAS(cscal)( &n, &alpha, x, &incx ); }
void Scal( int n, dcomplex alpha, dcomplex* x, int incx )
{ EL_BLAS(zscal)( &n, &alpha, x, &incx ); }

// NOTE: 'nrm1' is not the official name but is consistent with 'nrm2'
float Nrm1( int n, const float* x, int incx )
{ return EL_BLAS(sasum)( &n, x, &incx ); }
double Nrm1( int n, const double* x, int incx )
{ return EL_BLAS(dasum)( &n, x, &incx ); }
float Nrm1( int n, const scomplex* x, int incx )
{ return EL_LAPACK(scsum1)( &n, x, &incx ); }
double Nrm1( int n, const dcomplex* x, int incx )
{ return EL_LAPACK(dzsum1)( &n, x, &incx ); }

void Swap( int n, float* x, int incx, float* y, int incy )
{ EL_BLAS(sswap)( &n, x, &incx, y, &incy ); }
void Swap( int n, double* x, int incx, double* y, int incy )
{ EL_BLAS(dswap)( &n, x, &incx, y, &incy ); }
void Swap( int n, scomplex* x, int incx, scomplex* y, int incy )
{ EL_BLAS(cswap)( &n, x, &incx, y, &incy ); }
void Swap( int n, dcomplex* x, int incx, dcomplex* y, int incy )
{ EL_BLAS(zswap)( &n, x, &incx, y, &incy ); }

// Level 2 BLAS
// ============
void Gemv
( char trans, int m, int n,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, int m, int n,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ EL_BLAS(cgemv)
  ( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ EL_BLAS(zgemv)
  ( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Ger
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda  )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ EL_BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ EL_BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ EL_BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ EL_BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ EL_BLAS(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ EL_BLAS(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Hemv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ EL_BLAS(chemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ EL_BLAS(zhemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{ EL_BLAS(cher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{ EL_BLAS(zher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ EL_BLAS(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ EL_BLAS(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Symv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    EL_LAPACK(csymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Symv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    EL_LAPACK(zsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Syr
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda  )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Syr
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Syr
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    EL_LAPACK(csyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

void Syr
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    EL_LAPACK(zsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

void Syr2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Syr2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Syr2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{
    // csyr2 doesn't exist, so we route through csyr2k. However, csyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const int k = 1;
    const scomplex beta = 1.;
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

void Syr2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{
    // zsyr2 doesn't exist, so we route through zsyr2k. However, zsyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const int k = 1;
    const dcomplex beta = 1.;
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

void Trmv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ EL_BLAS(strmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ EL_BLAS(dtrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ EL_BLAS(ctrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ EL_BLAS(ztrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ EL_BLAS(strsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ EL_BLAS(dtrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ EL_BLAS(ctrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ EL_BLAS(ztrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

// Level 3 BLAS
// ============
void Gemm
( char transA, char transB, int m, int n, int k, 
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    EL_BLAS(sgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB,
  int m, int n, int k, 
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    EL_BLAS(dgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB, int m, int n, int k, 
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    EL_BLAS(cgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB, int m, int n, int k, 
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    EL_BLAS(zgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    EL_BLAS(chemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    EL_BLAS(zhemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(ssyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dsyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    EL_BLAS(cher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    EL_BLAS(zher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(ssyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dsyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ EL_BLAS(cherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Herk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ EL_BLAS(zherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Symm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    EL_BLAS(csymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    EL_BLAS(zsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    EL_BLAS(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    EL_BLAS(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syrk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{ EL_BLAS(ssyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{ EL_BLAS(dsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ EL_BLAS(csyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ EL_BLAS(zsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    EL_BLAS(strmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    EL_BLAS(dtrmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    EL_BLAS(ctrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    EL_BLAS(ztrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(strsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dtrsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    EL_BLAS(ctrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    EL_BLAS(ztrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

} // namespace blas
} // namespace El
