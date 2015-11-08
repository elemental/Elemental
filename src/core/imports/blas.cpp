/*
   Copyright (c) 2009-2015, Jack Poulson
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

// Level 1 BLAS
// ============
void EL_BLAS(saxpy)
( const BlasInt* n, const float* alpha, const float* x, const BlasInt* incx,
                                              float* y, const BlasInt* incy );
void EL_BLAS(daxpy)
( const BlasInt* n, const double* alpha, const double* x, const BlasInt* incx,
                                               double* y, const BlasInt* incy );
void EL_BLAS(caxpy)
( const BlasInt* n, const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zaxpy)
( const BlasInt* n, const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(scopy)
( const BlasInt* n, const float* x, const BlasInt* incx,
                          float* y, const BlasInt* incy );
void EL_BLAS(dcopy)
( const BlasInt* n, const double* x, const BlasInt* incx,
                          double* y, const BlasInt* incy );
void EL_BLAS(ccopy)
( const BlasInt* n, const scomplex* x, const BlasInt* incx,
                          scomplex* y, const BlasInt* incy );
void EL_BLAS(zcopy)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx,
                          dcomplex* y, const BlasInt* incy );

float EL_BLAS(sdot)
( const BlasInt* n, const float* x, const BlasInt* incx,
                    const float* y, const BlasInt* incy );
double EL_BLAS(ddot)
( const BlasInt* n, const double* x, const BlasInt* incx,
                    const double* y, const BlasInt* incy );
// To avoid the compatibility issue, we simply handroll our own complex dots
float  EL_BLAS(snrm2) 
( const BlasInt* n, const float   * x, const BlasInt* incx );
double EL_BLAS(dnrm2) 
( const BlasInt* n, const double  * x, const BlasInt* incx );
float  EL_BLAS(scnrm2)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
double EL_BLAS(dznrm2)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );

BlasInt EL_LAPACK(isamax)
( const BlasInt* n, const float* x, const BlasInt* incx );
BlasInt EL_LAPACK(idamax)
( const BlasInt* n, const double* x, const BlasInt* incx );
BlasInt EL_LAPACK(icamax)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
BlasInt EL_LAPACK(izamax)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );

// Apply a Givens rotation to a pair of vectors
void EL_BLAS(srot)
( const BlasInt* n, float* x, const BlasInt* incx, 
                    float* y, const BlasInt* incy,
  const float* c, const float* s );
void EL_BLAS(drot)
( const BlasInt* n, double* x, const BlasInt* incx, 
                    double* y, const BlasInt* incy,
  const double* c, const double* s );
void EL_BLAS(crot)
( const BlasInt* n, scomplex* x, const BlasInt* incx, 
                    scomplex* y, const BlasInt* incy,
  const float* c, const scomplex* s );
void EL_BLAS(zrot)
( const BlasInt* n, dcomplex* x, const BlasInt* incx, 
                    dcomplex* y, const BlasInt* incy,
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
( const BlasInt* n, const float   * alpha, float   * x, const BlasInt* incx );
void EL_BLAS(dscal)
( const BlasInt* n, const double  * alpha, double  * x, const BlasInt* incx );
void EL_BLAS(cscal)
( const BlasInt* n, const scomplex* alpha, scomplex* x, const BlasInt* incx );
void EL_BLAS(zscal)
( const BlasInt* n, const dcomplex* alpha, dcomplex* x, const BlasInt* incx );

float  EL_BLAS(sasum) 
( const BlasInt* n, const float   * x, const BlasInt* incx );
double EL_BLAS(dasum) 
( const BlasInt* n, const double  * x, const BlasInt* incx );
float  EL_BLAS(scasum)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
double EL_BLAS(dzasum)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );
float  EL_LAPACK(scsum1)
( const BlasInt* n, const scomplex* x, const BlasInt* incx );
double EL_LAPACK(dzsum1)
( const BlasInt* n, const dcomplex* x, const BlasInt* incx );

void EL_BLAS(sswap)
( const BlasInt* n, float   * x, const BlasInt* incx, 
                    float   * y, const BlasInt* incy );
void EL_BLAS(dswap)
( const BlasInt* n, double  * x, const BlasInt* incx, 
                    double  * y, const BlasInt* incy );
void EL_BLAS(cswap)
( const BlasInt* n, scomplex* x, const BlasInt* incx, 
                    scomplex* y, const BlasInt* incy );
void EL_BLAS(zswap)
( const BlasInt* n, dcomplex* x, const BlasInt* incx, 
                    dcomplex* y, const BlasInt* incy );

// Level 2 BLAS
// ============
void EL_BLAS(sgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* lda,
                      const float* x, const BlasInt* incx,
  const float* beta,        float* y, const BlasInt* incy );
void EL_BLAS(dgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* lda,
                       const double* x, const BlasInt* incx,
  const double* beta,        double* y, const BlasInt* incy );
void EL_BLAS(cgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(sger)
( const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* x, const BlasInt* incx,
                      const float* y, const BlasInt* incy,
                            float* A, const BlasInt* lda  );
void EL_BLAS(dger)
( const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* x, const BlasInt* incx,
                       const double* y, const BlasInt* incy,
                             double* A, const BlasInt* lda  );
void EL_BLAS(cgerc)
( const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* lda  );
void EL_BLAS(zgerc)
( const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* lda  );

void EL_BLAS(cgeru)
( const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* lda  );
void EL_BLAS(zgeru)
( const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* lda  );

void EL_BLAS(chemv)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zhemv)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(cher)
( const char* uplo, const BlasInt* m,
  const float* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* A, const BlasInt* lda  );
void EL_BLAS(zher)
( const char* uplo, const BlasInt* m,
  const double* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* A, const BlasInt* lda  );

void EL_BLAS(cher2)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* lda  );
void EL_BLAS(zher2)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* lda  );

void EL_BLAS(ssymv)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* A, const BlasInt* lda,
                      const float* x, const BlasInt* incx,
  const float* beta,        float* y, const BlasInt* incy );
void EL_BLAS(dsymv)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* A, const BlasInt* lda,
                       const double* x, const BlasInt* incx,
  const double* beta,        double* y, const BlasInt* incy );
// 'csymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csymv)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
// 'zsymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsymv)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(ssyr)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* x, const BlasInt* incx,
                            float* A, const BlasInt* lda  );
void EL_BLAS(dsyr)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* x, const BlasInt* incx,
                             double* A, const BlasInt* lda  );
// 'csyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csyr)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* A, const BlasInt* lda  );
// 'zsyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsyr)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* A, const BlasInt* lda  );

void EL_BLAS(ssyr2)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* x, const BlasInt* incx,
                      const float* y, const BlasInt* incy,
                            float* A, const BlasInt* lda  );
void EL_BLAS(dsyr2)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* x, const BlasInt* incx,
                       const double* y, const BlasInt* incy,
                             double* A, const BlasInt* lda  );

void EL_BLAS(strmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const float* A, const BlasInt* lda, float* x, const BlasInt* incx );
void EL_BLAS(dtrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const double* A, const BlasInt* lda, double* x, const BlasInt* incx );
void EL_BLAS(ctrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const scomplex* A, const BlasInt* lda,
        scomplex* x, const BlasInt* incx );
void EL_BLAS(ztrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const dcomplex* A, const BlasInt* lda,
        dcomplex* x, const BlasInt* incx );

void EL_BLAS(strsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const float* A, const BlasInt* lda, float* x, const BlasInt* incx );
void EL_BLAS(dtrsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const double* A, const BlasInt* lda, double* x, const BlasInt* incx );
void EL_BLAS(ctrsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const scomplex* A, const BlasInt* lda,
        scomplex* x, const BlasInt* incx );
void EL_BLAS(ztrsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const dcomplex* A, const BlasInt* lda,
        dcomplex* x, const BlasInt* incx );

// Level 3 BLAS 
// ============
void EL_BLAS(sgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* lda,
                      const float* B, const BlasInt* ldb,
  const float* beta,        float* C, const BlasInt* ldc );
void EL_BLAS(dgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* lda,
                       const double* B, const BlasInt* ldb,
  const double* beta,        double* C, const BlasInt* ldc );
void EL_BLAS(cgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* B, const BlasInt* ldb,
  const scomplex* beta,
        scomplex* C, const BlasInt* ldc );
void EL_BLAS(zgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* B, const BlasInt* ldb,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* ldc );

void EL_BLAS(chemm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* B, const BlasInt* ldb,
  const scomplex* beta,
        scomplex* C, const BlasInt* ldc );
void EL_BLAS(zhemm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* B, const BlasInt* ldb,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* ldc );

void EL_BLAS(cher2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* B, const BlasInt* ldb,
  const float* beta,
        scomplex* C, const BlasInt* ldc );
void EL_BLAS(zher2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* B, const BlasInt* ldb,
  const double* beta,
        dcomplex* C, const BlasInt* ldc );

void EL_BLAS(cherk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha,
  const scomplex* A, const BlasInt* lda,
  const float* beta,
        scomplex* C, const BlasInt* ldc );
void EL_BLAS(zherk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha,
  const dcomplex* A, const BlasInt* lda,
  const double* beta,
        dcomplex* C, const BlasInt* ldc );

void EL_BLAS(ssymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* lda,
                      const float* B, const BlasInt* ldb,
  const float* beta,        float* C, const BlasInt* ldc );
void EL_BLAS(dsymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* lda,
                       const double* B, const BlasInt* ldb,
  const double* beta,        double* C, const BlasInt* ldc );
void EL_BLAS(csymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* B, const BlasInt* ldb,
  const scomplex* beta,
        scomplex* C, const BlasInt* ldc );
void EL_BLAS(zsymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* B, const BlasInt* ldb,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* ldc );

void EL_BLAS(ssyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* lda,
                      const float* B, const BlasInt* ldb,
  const float* beta,        float* C, const BlasInt* ldc );
void EL_BLAS(dsyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* lda,
                       const double* B, const BlasInt* ldb,
  const double* beta,        double* C, const BlasInt* ldc );
void EL_BLAS(csyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* B, const BlasInt* ldb,
  const scomplex* beta,
        scomplex* C, const BlasInt* ldc );
void EL_BLAS(zsyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* B, const BlasInt* ldb,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* ldc );

void EL_BLAS(ssyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* lda,
  const float* beta,        float* C, const BlasInt* ldc );
void EL_BLAS(dsyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* lda,
  const double* beta,        double* C, const BlasInt* ldc );
void EL_BLAS(csyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
  const scomplex* beta,
        scomplex* C, const BlasInt* ldc );
void EL_BLAS(zsyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* ldc );

void EL_BLAS(strmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* lda,
                            float* B, const BlasInt* ldb );
void EL_BLAS(dtrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* lda,
                             double* B, const BlasInt* ldb );
void EL_BLAS(ctrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
        scomplex* B, const BlasInt* ldb );
void EL_BLAS(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
        dcomplex* B, const BlasInt* ldb );

void EL_BLAS(strsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* lda,
                            float* B, const BlasInt* ldb );
void EL_BLAS(dtrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* lda,
                             double* B, const BlasInt* ldb );
void EL_BLAS(ctrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* lda,
        scomplex* B, const BlasInt* ldb );
void EL_BLAS(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* lda,
        dcomplex* B, const BlasInt* ldb );

} // extern "C"

namespace El {
namespace blas {

// Level 1 BLAS
// ============
template<typename T>
void Axpy( BlasInt n, T alpha, const T* x, BlasInt incx, T* y, BlasInt incy )
{
    for( BlasInt i=0; i<n; ++i )
        y[i*incy] += alpha*x[i*incx];
}
template void Axpy
( BlasInt n, Int alpha, const Int* x, BlasInt incx, Int* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template void Axpy
( BlasInt n, Quad alpha, 
  const Quad* x, BlasInt incx, Quad* y, BlasInt incy );
template void Axpy
( BlasInt n, Complex<Quad> alpha, 
  const Complex<Quad>* x, BlasInt incx, Complex<Quad>* y, BlasInt incy );
#endif

void Axpy
( BlasInt n, float alpha, const float* x, BlasInt incx, 
                                float* y, BlasInt incy )
{ EL_BLAS(saxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n, double alpha, const double* x, BlasInt incx, 
                                 double* y, BlasInt incy )
{ EL_BLAS(daxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n, scomplex alpha, const scomplex* x, BlasInt incx, 
                                   scomplex* y, BlasInt incy )
{ EL_BLAS(caxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n, dcomplex alpha, const dcomplex* x, BlasInt incx, 
                                   dcomplex* y, BlasInt incy )
{ EL_BLAS(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }

template<typename T>
void Copy( BlasInt n, const T* x, BlasInt incx, T* y, BlasInt incy )
{
    for( BlasInt i=0; i<n; ++i )
        y[i*incy] = x[i*incx];
}
template void Copy
( BlasInt n, const Int* x, BlasInt incx, Int* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template void Copy
( BlasInt n, const Quad* x, BlasInt incx, 
                   Quad* y, BlasInt incy );
template void Copy
( BlasInt n, const Complex<Quad>* x, BlasInt incx, 
                   Complex<Quad>* y, BlasInt incy );
#endif

void Copy
( BlasInt n, const float* x, BlasInt incx, float* y, BlasInt incy )
{ EL_BLAS(scopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n, const double* x, BlasInt incx, double* y, BlasInt incy )
{ EL_BLAS(dcopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n, const scomplex* x, BlasInt incx, scomplex* y, BlasInt incy )
{ EL_BLAS(ccopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n, const dcomplex* x, BlasInt incx, dcomplex* y, BlasInt incy )
{ EL_BLAS(zcopy)( &n, x, &incx, y, &incy ); }

template<typename T>
T Dot( BlasInt n, const T* x, BlasInt incx, const T* y, BlasInt incy )
{
    T alpha = 0;
    for( BlasInt i=0; i<n; ++i )
        alpha += Conj(x[i*incx])*y[i*incy];
    return alpha;
}
template Int Dot
( BlasInt n, const Int* x, BlasInt incx, const Int* y, BlasInt incy );
template float Dot
( BlasInt n, const float* x, BlasInt incx, const float* y, BlasInt incy );
template scomplex Dot
( BlasInt n, const scomplex* x, BlasInt incx, const scomplex* y, BlasInt incy );
template dcomplex Dot
( BlasInt n, const dcomplex* x, BlasInt incx, const dcomplex* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template Quad Dot
( BlasInt n, const Quad* x, BlasInt incx, 
             const Quad* y, BlasInt incy );
template Complex<Quad> Dot
( BlasInt n, const Complex<Quad>* x, BlasInt incx, 
             const Complex<Quad>* y, BlasInt incy );
#endif

// NOTE: I am under the impression that it is generally unsafe to return 
//       anything except a double-precision float to C from Fortran
double Dot
( BlasInt n, const double* x, BlasInt incx, const double* y, BlasInt incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }

template<typename T>
T Dotu( BlasInt n, const T* x, BlasInt incx, const T* y, BlasInt incy )
{
    T alpha = 0;
    for( BlasInt i=0; i<n; ++i )
        alpha += x[i*incx]*y[i*incy];
    return alpha;
}
template Int Dotu
( BlasInt n, const Int* x, BlasInt incx, const Int* y, BlasInt incy );
template float Dotu
( BlasInt n, const float* x, BlasInt incx, const float* y, BlasInt incy );
template scomplex Dotu
( BlasInt n, const scomplex* x, BlasInt incx, const scomplex* y, BlasInt incy );
template dcomplex Dotu
( BlasInt n, const dcomplex* x, BlasInt incx, const dcomplex* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template Quad Dotu
( BlasInt n, const Quad* x, BlasInt incx, 
             const Quad* y, BlasInt incy );
template Complex<Quad> Dotu
( BlasInt n, const Complex<Quad>* x, BlasInt incx, 
             const Complex<Quad>* y, BlasInt incy );
#endif

double Dotu
( BlasInt n, const double* x, BlasInt incx, const double* y, BlasInt incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }

template<typename F>
Base<F> Nrm2( BlasInt n, const F* x, BlasInt incx )
{
    typedef Base<F> Real;
    Real scale = 0; 
    Real scaledSquare = 1;
    for( BlasInt i=0; i<n; ++i )
        UpdateScaledSquare( x[i*incx], scale, scaledSquare );
    return scale*Sqrt(scaledSquare);
}
template float Nrm2( BlasInt n, const float* x, BlasInt incx );
template float Nrm2( BlasInt n, const scomplex* x, BlasInt incx );
#ifdef EL_HAVE_QUAD
template Quad Nrm2( BlasInt n, const Quad* x, BlasInt incx );
template Quad Nrm2( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif

double Nrm2( BlasInt n, const double* x, BlasInt incx )
{ return EL_BLAS(dnrm2)( &n, x, &incx ); }
double Nrm2( BlasInt n, const dcomplex* x, BlasInt incx )
{ return EL_BLAS(dznrm2)( &n, x, &incx ); }

template<typename F>
BlasInt MaxInd( BlasInt n, const F* x, BlasInt incx )
{
    typedef Base<F> Real;
    Real maxAbsVal = -1;
    BlasInt maxAbsInd = -1;
    for( BlasInt i=0; i<n; ++i ) 
    {
        const Real absVal = Abs(x[i*incx]);
        if( absVal > maxAbsVal )
        {
            maxAbsVal = absVal;
            maxAbsInd = i;
        }
    } 
    return maxAbsInd;
}
#ifdef EL_HAVE_QUAD
template BlasInt MaxInd( BlasInt n, const Int* x, BlasInt incx );
template BlasInt MaxInd( BlasInt n, const Quad* x, BlasInt incx );
template BlasInt MaxInd( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif

BlasInt MaxInd( BlasInt n, const float* x, BlasInt incx )
{
    return EL_LAPACK(isamax)( &n, x, &incx ) - 1;
}

BlasInt MaxInd( BlasInt n, const double* x, BlasInt incx )
{
    return EL_LAPACK(idamax)( &n, x, &incx ) - 1;
}

BlasInt MaxInd( BlasInt n, const scomplex* x, BlasInt incx )
{
    return EL_LAPACK(icamax)( &n, x, &incx ) - 1;
}

BlasInt MaxInd( BlasInt n, const dcomplex* x, BlasInt incx )
{
    return EL_LAPACK(izamax)( &n, x, &incx ) - 1;
}

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
( BlasInt n, float* x, BlasInt incx, 
             float* y, BlasInt incy, float c, float s )
{ EL_BLAS(srot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n, double* x, BlasInt incx, 
             double* y, BlasInt incy, double c, double s )
{ EL_BLAS(drot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n, scomplex* x, BlasInt incx, 
             scomplex* y, BlasInt incy, float c, scomplex s )
{ EL_BLAS(crot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n, dcomplex* x, BlasInt incx, 
             dcomplex* y, BlasInt incy, double c, dcomplex s )
{ EL_BLAS(zrot)( &n, x, &incx, y, &incy, &c, &s ); }

template<typename T>
void Scal( BlasInt n, T alpha, T* x, BlasInt incx )
{
    for( BlasInt j=0; j<n; ++j )
        x[j*incx] *= alpha;
}
template<typename T>
void Scal( BlasInt n, T alpha, Complex<T>* x, BlasInt incx )
{
    for( BlasInt j=0; j<n; ++j )
        x[j*incx] *= alpha;
}
template void Scal( BlasInt n, Int alpha, Int* x, BlasInt incx );
#ifdef EL_HAVE_QUAD
template void Scal
( BlasInt n, Quad alpha, Quad* x, BlasInt incx );
template void Scal
( BlasInt n, Complex<Quad> alpha, Complex<Quad>* x, BlasInt incx );
#endif
template void Scal
( BlasInt n, float alpha, Complex<float>* x, BlasInt incx );
template void Scal
( BlasInt n, double alpha, Complex<double>* x, BlasInt incx );
#ifdef EL_HAVE_QUAD
template void Scal
( BlasInt n, Quad alpha, Complex<Quad>* x, BlasInt incx );
#endif

void Scal( BlasInt n, float alpha, float* x, BlasInt incx )
{ EL_BLAS(sscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, double alpha, double* x, BlasInt incx )
{ EL_BLAS(dscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, scomplex alpha, scomplex* x, BlasInt incx )
{ EL_BLAS(cscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, dcomplex alpha, dcomplex* x, BlasInt incx )
{ EL_BLAS(zscal)( &n, &alpha, x, &incx ); }

// NOTE: 'nrm1' is not the official name but is consistent with 'nrm2'
template<typename F>
Base<F> Nrm1( BlasInt n, const F* x, BlasInt incx )
{
    Base<F> sum=0;
    for( BlasInt i=0; i<n; ++i )
        sum += Abs(x[i*incx]);
    return sum;
}
template float Nrm1( BlasInt n, const float* x, BlasInt incx );
template float Nrm1( BlasInt n, const scomplex* x, BlasInt incx );
#ifdef EL_HAVE_QUAD
template Quad Nrm1( BlasInt n, const Quad* x, BlasInt incx );
template Quad Nrm1( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif

double Nrm1( BlasInt n, const double* x, BlasInt incx )
{ return EL_BLAS(dasum)( &n, x, &incx ); }
double Nrm1( BlasInt n, const dcomplex* x, BlasInt incx )
{ return EL_LAPACK(dzsum1)( &n, x, &incx ); }

template<typename T>
void Swap( BlasInt n, T* x, BlasInt incx, T* y, BlasInt incy )
{
    for( BlasInt i=0; i<n; ++i )
    {
        const T temp = x[i*incx];
        x[i*incx] = y[i*incy];
        y[i*incy] = temp;
    }
}
template void Swap( BlasInt n, Int* x, BlasInt incx, Int* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template void Swap( BlasInt n, Quad* x, BlasInt incx, Quad* y, BlasInt incy );
template void Swap
( BlasInt n, Complex<Quad>* x, BlasInt incx, Complex<Quad>* y, BlasInt incy );
#endif

void Swap( BlasInt n, float* x, BlasInt incx, float* y, BlasInt incy )
{ EL_BLAS(sswap)( &n, x, &incx, y, &incy ); }
void Swap( BlasInt n, double* x, BlasInt incx, double* y, BlasInt incy )
{ EL_BLAS(dswap)( &n, x, &incx, y, &incy ); }
void Swap( BlasInt n, scomplex* x, BlasInt incx, scomplex* y, BlasInt incy )
{ EL_BLAS(cswap)( &n, x, &incx, y, &incy ); }
void Swap( BlasInt n, dcomplex* x, BlasInt incx, dcomplex* y, BlasInt incy )
{ EL_BLAS(zswap)( &n, x, &incx, y, &incy ); }

// Level 2 BLAS
// ============
template<typename T>
void Gemv
( char trans, BlasInt m, BlasInt n,
  T alpha, const T* A, BlasInt lda, const T* x, BlasInt incx,
  T beta,        T* y, BlasInt incy )
{
    if( trans == 'N' )
    {
        if( m > 0 && n == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<m; ++i )
                y[i*incy] = 0;
            return;
        }
        Scal( m, beta, y, incy );
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                y[i*incy] += alpha*A[i+j*lda]*x[j*incx];
    }
    else if( trans == 'T' )
    {
        if( n > 0 && m == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<n; ++i )
                y[i*incy] = 0;
            return;
        }
        Scal( n, beta, y, incy );
        for( BlasInt i=0; i<n; ++i )
            for( BlasInt j=0; j<m; ++j )
                y[i*incy] += alpha*A[j+i*lda]*x[j*incx];
    }
    else
    {
        if( n > 0 && m == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<n; ++i )
                y[i*incy] = 0;
            return;
        }
        Scal( n, beta, y, incy );
        for( BlasInt i=0; i<n; ++i )
            for( BlasInt j=0; j<m; ++j )
                y[i*incy] += alpha*Conj(A[j+i*lda])*x[j*incx];
    }
}
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  Int alpha, const Int* A, BlasInt lda, const Int* x, BlasInt incx, 
  Int beta, Int* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  Quad alpha, const Quad* A, BlasInt lda, const Quad* x, BlasInt incx, 
  Quad beta, Quad* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda, 
                       const Complex<Quad>* x, BlasInt incx, 
  Complex<Quad> beta,        Complex<Quad>* y, BlasInt incy );
#endif

void Gemv
( char trans, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, const float* x, BlasInt incx,
  float beta,        float* y, BlasInt incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, const double* x, BlasInt incx,
  double beta,        double* y, BlasInt incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* x, BlasInt incx,
  scomplex beta,        scomplex* y, BlasInt incy )
{ EL_BLAS(cgemv)
  ( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Gemv
( char trans, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* x, BlasInt incx,
  dcomplex beta,        dcomplex* y, BlasInt incy )
{ EL_BLAS(zgemv)
  ( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

template<typename T>
void Ger
( BlasInt m, BlasInt n,
  T alpha, const T* x, BlasInt incx, 
           const T* y, BlasInt incy,
                 T* A, BlasInt lda )
{
    for( BlasInt j=0; j<n; ++j )
        for( BlasInt i=0; i<m; ++i )
            A[i+j*lda] += alpha*x[i*incx]*Conj(y[j*incy]);
}
template void Ger
( BlasInt m, BlasInt n, 
  Int alpha, const Int* x, BlasInt incx, 
             const Int* y, BlasInt incy, 
                   Int* A, BlasInt lda );
#ifdef EL_HAVE_QUAD
template void Ger
( BlasInt m, BlasInt n, 
  Quad alpha, const Quad* x, BlasInt incx, 
              const Quad* y, BlasInt incy, 
                    Quad* A, BlasInt lda );
template void Ger
( BlasInt m, BlasInt n, 
  Complex<Quad> alpha, const Complex<Quad>* x, BlasInt incx, 
                       const Complex<Quad>* y, BlasInt incy, 
                             Complex<Quad>* A, BlasInt lda );
#endif

void Ger
( BlasInt m, BlasInt n,
  float alpha, const float* x, BlasInt incx, 
               const float* y, BlasInt incy,
                     float* A, BlasInt lda )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( BlasInt m, BlasInt n,
  double alpha, const double* x, BlasInt incx, 
                const double* y, BlasInt incy,
                      double* A, BlasInt lda  )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* x, BlasInt incx, 
                  const scomplex* y, BlasInt incy,
                        scomplex* A, BlasInt lda )
{ EL_BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* x, BlasInt incx, 
                  const dcomplex* y, BlasInt incy,
                        dcomplex* A, BlasInt lda )
{ EL_BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

template<typename T>
void Geru
( BlasInt m, BlasInt n,
  T alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt lda )
{
    for( BlasInt j=0; j<n; ++j )
        for( BlasInt i=0; i<m; ++i )
            A[i+j*lda] += alpha*x[i*incx]*y[j*incy];
}
template void Geru
( BlasInt m, BlasInt n, 
  Int alpha,
  const Int* x, BlasInt incx, 
  const Int* y, BlasInt incy, 
        Int* A, BlasInt lda );
#ifdef EL_HAVE_QUAD
template void Geru
( BlasInt m, BlasInt n, 
  Quad alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt lda );
template void Geru
( BlasInt m, BlasInt n, 
  Complex<Quad> alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt lda );
#endif

void Geru
( BlasInt m, BlasInt n,
  float alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt lda )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( BlasInt m, BlasInt n,
  double alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt lda )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( BlasInt m, BlasInt n,
  scomplex alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt lda )
{ EL_BLAS(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( BlasInt m, BlasInt n,
  dcomplex alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt lda )
{ EL_BLAS(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

// TODO: Introduce some sort of blocking
template<typename T>
void Hemv
( char uplo, BlasInt m,
  T alpha, const T* A, BlasInt lda, 
           const T* x, BlasInt incx,
  T beta,        T* y, BlasInt incy )
{
    // y := beta*y
    if( beta == T(0) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] *= beta;
    }

    if( uplo == LOWER )
    {
        // Multiply with the lower triangle
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=j; i<m; ++i )
                y[i*incy] += alpha*A[i+j*lda]*x[j*incx];

        // Multiply with the adjoint of the strictly lower triangle
        for( BlasInt j=0; j<m; ++j ) 
            for( BlasInt i=j+1; i<m; ++i )
                y[j*incy] += alpha*Conj(A[i+j*lda])*x[i*incx];
    }
    else
    {
        // Multiply with the upper triangle
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=0; i<=j; ++i )
                y[i*incy] += alpha*A[i+j*lda]*x[j*incx];

        // Multiply with the adjoint of the strictly upper triangle
        for( BlasInt j=0; j<m; ++j ) 
            for( BlasInt i=0; i<j; ++i )
                y[j*incy] += alpha*Conj(A[i+j*lda])*x[i*incx];
    }
}

template void Hemv
( char uplo, BlasInt m, 
  Int alpha, const Int* A, BlasInt lda, 
             const Int* x, BlasInt incx, 
  Int beta,        Int* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template void Hemv
( char uplo, BlasInt m, 
  Quad alpha, const Quad* A, BlasInt lda, 
              const Quad* x, BlasInt incx, 
  Quad beta,        Quad* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda, 
                       const Complex<Quad>* x, BlasInt incx, 
  Complex<Quad> beta,        Complex<Quad>* y, BlasInt incy );
#endif

void Hemv
( char uplo, BlasInt m,
  float alpha, const float* A, BlasInt lda, 
               const float* x, BlasInt incx,
  float beta,        float* y, BlasInt incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  double alpha, const double* A, BlasInt lda, 
                const double* x, BlasInt incx,
  double beta,        double* y, BlasInt incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* x, BlasInt incx,
  scomplex beta,        scomplex* y, BlasInt incy )
{ EL_BLAS(chemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* x, BlasInt incx,
  dcomplex beta,        dcomplex* y, BlasInt incy )
{ EL_BLAS(zhemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

template<typename T>
void Her
( char uplo, BlasInt m, Base<T> alpha, 
  const T* x, BlasInt incx, T* A, BlasInt lda )
{
    if( uplo == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=j; i<m; ++i )
                A[i+j*lda] += alpha*x[i*incx]*Conj(x[j*incx]);
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=0; i<=j; ++i )
                A[i+j*lda] += alpha*x[i*incx]*Conj(x[j*incx]);
    }
}
template void Her
( char uplo, BlasInt m, Int alpha, const Int* x, BlasInt incx, 
                                         Int* A, BlasInt lda );
#ifdef EL_HAVE_QUAD
template void Her
( char uplo, BlasInt m, 
  Quad alpha, const Quad* x, BlasInt incx, Quad* A, BlasInt lda );
template void Her
( char uplo, BlasInt m, 
  Quad alpha, const Complex<Quad>* x, BlasInt incx, 
                    Complex<Quad>* A, BlasInt lda );
#endif

void Her
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, float* A, BlasInt lda )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, double* A, BlasInt lda )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, BlasInt m,
  float alpha, const scomplex* x, BlasInt incx, scomplex* A, BlasInt lda )
{ EL_BLAS(cher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, BlasInt m,
  double alpha, const dcomplex* x, BlasInt incx, dcomplex* A, BlasInt lda )
{ EL_BLAS(zher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

template<typename T>
void Her2
( char uplo, BlasInt m, 
  T alpha, const T* x, BlasInt incx, const T* y, BlasInt incy,
                 T* A, BlasInt lda )
{
    if( uplo == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=j; i<m; ++i )
                A[i+j*lda] += alpha*      x[i*incx]*Conj(y[j*incy]) + 
                              Conj(alpha)*y[i*incy]*Conj(x[j*incx]);
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=0; i<=j; ++i )
                A[i+j*lda] += alpha*      x[i*incx]*Conj(y[j*incy]) + 
                              Conj(alpha)*y[i*incy]*Conj(x[j*incx]);
    }
}
template void Her2
( char uplo, BlasInt m, 
  Int alpha, const Int* x, BlasInt incx, 
             const Int* y, BlasInt incy, 
                   Int* A, BlasInt lda );
#ifdef EL_HAVE_QUAD
template void Her2
( char uplo, BlasInt m, 
  Quad alpha, const Quad* x, BlasInt incx, 
              const Quad* y, BlasInt incy, 
                    Quad* A, BlasInt lda );
template void Her2
( char uplo, BlasInt m, 
  Complex<Quad> alpha, const Complex<Quad>* x, BlasInt incx, 
                       const Complex<Quad>* y, BlasInt incy, 
                             Complex<Quad>* A, BlasInt lda );
#endif

void Her2
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, 
               const float* y, BlasInt incy,
                     float* A, BlasInt lda )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, 
                const double* y, BlasInt incy,
                      double* A, BlasInt lda )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* x, BlasInt incx, 
                  const scomplex* y, BlasInt incy,
                        scomplex* A, BlasInt lda )
{ EL_BLAS(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* x, BlasInt incx, 
                  const dcomplex* y, BlasInt incy,
                        dcomplex* A, BlasInt lda )
{ EL_BLAS(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

// TODO: Introduce some sort of blocking
template<typename T>
void Symv
( char uplo, BlasInt m,
  T alpha, const T* A, BlasInt lda, 
           const T* x, BlasInt incx,
  T beta,        T* y, BlasInt incy )
{
    // y := beta*y
    if( beta == T(0) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] *= beta;
    }

    if( uplo == LOWER )
    {
        // Multiply with the lower triangle
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=j; i<m; ++i )
                y[i*incy] += alpha*A[i+j*lda]*x[j*incx];

        // Multiply with the transpose of the strictly lower triangle
        for( BlasInt j=0; j<m; ++j ) 
            for( BlasInt i=j+1; i<m; ++i )
                y[j*incy] += alpha*A[i+j*lda]*x[i*incx];
    }
    else
    {
        // Multiply with the upper triangle
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=0; i<=j; ++i )
                y[i*incy] += alpha*A[i+j*lda]*x[j*incx];

        // Multiply with the transpose of the strictly upper triangle
        for( BlasInt j=0; j<m; ++j ) 
            for( BlasInt i=0; i<j; ++i )
                y[j*incy] += alpha*A[i+j*lda]*x[i*incx];
    }
}

template void Symv
( char uplo, BlasInt m, 
  Int alpha, const Int* A, BlasInt lda, 
             const Int* x, BlasInt incx, 
  Int beta,        Int* y, BlasInt incy );
#ifdef EL_HAVE_QUAD
template void Symv
( char uplo, BlasInt m, 
  Quad alpha, const Quad* A, BlasInt lda, 
              const Quad* x, BlasInt incx, 
  Quad beta,        Quad* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda, 
                       const Complex<Quad>* x, BlasInt incx, 
  Complex<Quad> beta,        Complex<Quad>* y, BlasInt incy );
#endif

void Symv
( char uplo, BlasInt m,
  float alpha, const float* A, BlasInt lda, 
               const float* x, BlasInt incx,
  float beta,        float* y, BlasInt incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, BlasInt m,
  double alpha, const double* A, BlasInt lda, 
                const double* x, BlasInt incx,
  double beta,        double* y, BlasInt incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* x, BlasInt incx,
  scomplex beta,        scomplex* y, BlasInt incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    EL_LAPACK(csymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Symv
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* x, BlasInt incx,
  dcomplex beta,        dcomplex* y, BlasInt incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    EL_LAPACK(zsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

template<typename T>
void Syr
( char uplo, BlasInt m, T alpha, const T* x, BlasInt incx, T* A, BlasInt lda )
{
    if( uplo == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=j; i<m; ++i )
                A[i+j*lda] += alpha*x[i*incx]*x[j*incx];
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=0; i<=j; ++i )
                A[i+j*lda] += alpha*x[i*incx]*x[j*incx];
    }
}
template void Syr
( char uplo, BlasInt m, Int alpha, const Int* x, BlasInt incx, 
                                         Int* A, BlasInt lda );
#ifdef EL_HAVE_QUAD
template void Syr
( char uplo, BlasInt m, 
  Quad alpha, const Quad* x, BlasInt incx, Quad* A, BlasInt lda );
template void Syr
( char uplo, BlasInt m, 
  Complex<Quad> alpha, const Complex<Quad>* x, BlasInt incx, 
                             Complex<Quad>* A, BlasInt lda );
#endif

void Syr
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, float* A, BlasInt lda  )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Syr
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, double* A, BlasInt lda )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Syr
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* x, BlasInt incx, scomplex* A, BlasInt lda )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    EL_LAPACK(csyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

void Syr
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* x, BlasInt incx, dcomplex* A, BlasInt lda )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    EL_LAPACK(zsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

template<typename T>
void Syr2
( char uplo, BlasInt m, 
  T alpha, const T* x, BlasInt incx, const T* y, BlasInt incy,
                 T* A, BlasInt lda )
{
    if( uplo == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=j; i<m; ++i )
                A[i+j*lda] += alpha*(x[i*incx]*y[j*incy]+y[i*incy]*x[j*incx]);
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
            for( BlasInt i=0; i<=j; ++i )
                A[i+j*lda] += alpha*(x[i*incx]*y[j*incy]+y[i*incy]*x[j*incx]);
    }
}
template void Syr2
( char uplo, BlasInt m, 
  Int alpha, const Int* x, BlasInt incx, const Int* y, BlasInt incy, 
                   Int* A, BlasInt lda );
#ifdef EL_HAVE_QUAD
template void Syr2
( char uplo, BlasInt m, 
  Quad alpha, const Quad* x, BlasInt incx, 
              const Quad* y, BlasInt incy, 
                    Quad* A, BlasInt lda );
template void Syr2
( char uplo, BlasInt m, 
  Complex<Quad> alpha, const Complex<Quad>* x, BlasInt incx, 
                       const Complex<Quad>* y, BlasInt incy, 
                             Complex<Quad>* A, BlasInt lda );
#endif

void Syr2
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, const float* y, BlasInt incy,
                     float* A, BlasInt lda )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Syr2
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, const double* y, BlasInt incy,
                      double* A, BlasInt lda )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Syr2
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* x, BlasInt incx, 
                  const scomplex* y, BlasInt incy,
                        scomplex* A, BlasInt lda )
{
    // csyr2 doesn't exist, so we route through csyr2k. However, csyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const BlasInt k = 1;
    const scomplex beta = 1.;
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

void Syr2
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* x, BlasInt incx, 
                  const dcomplex* y, BlasInt incy,
                        dcomplex* A, BlasInt lda )
{
    // zsyr2 doesn't exist, so we route through zsyr2k. However, zsyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const BlasInt k = 1;
    const dcomplex beta = 1.;
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

// TODO: Verify correctness and BlasIntroduce blocking
template<typename T>
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const T* A, BlasInt lda, T* x, BlasInt incx )
{
    const bool conj = ( trans == 'C' );
    if( uplo == LOWER )
    {
        if( trans == 'N' )
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                T gamma = x[j*incx]; 
                if( gamma != T(0) )
                {
                    for( BlasInt i=m-1; i>j; --i )
                        x[i*incx] += gamma*A[i+j*lda];
                    if( diag != 'N' )
                        x[j*incx] *= A[j+j*lda];
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<m; ++j )
            {
                T gamma = x[j*incx]; 
                if( conj )
                {
                    if( diag != 'N' )
                        gamma *= Conj(A[j+j*lda]);
                    for( BlasInt i=j+1; i<m; ++i )
                        gamma += Conj(A[i+j*lda])*x[i*incx];
                }
                else
                {
                    if( diag != 'N' )
                        gamma *= A[j+j*lda];
                    for( BlasInt i=j+1; i<m; ++i )
                        gamma += A[i+j*lda]*x[i*incx];
                }
                x[j*incx] = gamma;
            }
        }
    }
    else
    {
        if( trans == 'N' )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                T gamma = x[j*incx]; 
                if( gamma != T(0) )
                {
                    for( BlasInt i=0; i<j; ++i )
                        x[i*incx] += gamma*A[i+j*lda];
                    if( diag != 'N' )
                        x[j*incx] *= A[j+j*lda];
                }
            }
        }
        else
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                T gamma = x[j*incx]; 
                if( conj )
                {
                    if( diag != 'N' )
                        gamma *= Conj(A[j+j*lda]);
                    for( BlasInt i=j-1; i>=0; --i )
                        gamma += Conj(A[i+j*lda])*x[i*incx];
                }
                else
                {
                    if( diag != 'N' )
                        gamma *= A[j+j*lda];
                    for( BlasInt i=j-1; i>=0; --i )
                        gamma += A[i+j*lda]*x[i*incx];
                }
                x[j*incx] = gamma;
            }
        }
    }
}
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Int* A, BlasInt lda, Int* x, BlasInt incx );
#ifdef EL_HAVE_QUAD
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Quad* A, BlasInt lda, Quad* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<Quad>* A, BlasInt lda, Complex<Quad>* x, BlasInt incx );
#endif

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt lda, float* x, BlasInt incx )
{ EL_BLAS(strmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt lda, double* x, BlasInt incx )
{ EL_BLAS(dtrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt lda, scomplex* x, BlasInt incx )
{ EL_BLAS(ctrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt lda, dcomplex* x, BlasInt incx )
{ EL_BLAS(ztrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

template<typename F>
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const F* A, BlasInt lda, F* x, BlasInt incx )
{
    const bool conj = ( trans == 'C' );
    if( uplo == 'L' )
    {
        if( trans == 'N' )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                if( diag == 'N' ) 
                    x[j*incx] /= A[j+j*lda];
                const F gamma = x[j*incx];
                for( BlasInt i=j+1; i<m; ++i )
                    x[i*incx] -= gamma*A[i+j*lda];
            }
        }
        else
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                if( conj )
                {
                    if( diag == 'N' )
                        x[j*incx] /= Conj(A[j+j*lda]);
                    const F gamma = x[j*incx];
                    for( BlasInt i=j-1; i>=0; --i )
                        x[i*incx] -= gamma*Conj(A[j+i*lda]);
                }
                else
                {
                    if( diag == 'N' )
                        x[j*incx] /= A[j+j*lda];
                    const F gamma = x[j*incx];
                    for( BlasInt i=j-1; i>=0; --i )
                        x[i*incx] -= gamma*A[j+i*lda];
                }
            }
        }
    }
    else
    {
        if( trans == 'N' )
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                if( diag == 'N' )
                    x[j*incx] /= A[j+j*lda];
                const F gamma = x[j*incx];
                for( BlasInt i=j-1; i>=0; --i )
                    x[i*incx] -= gamma*A[i+j*lda];
            }
        }
        else
        {
            for( BlasInt j=0; j<m; ++j )
            {
                if( conj )
                {
                    if( diag == 'N' ) 
                        x[j*incx] /= A[j+j*lda];
                    const F gamma = x[j*incx];
                    for( BlasInt i=j+1; i<m; ++i )
                        x[i*incx] -= gamma*A[j+i*lda];
                }
                else
                {
                    if( diag == 'N' ) 
                        x[j*incx] /= Conj(A[j+j*lda]);
                    const F gamma = x[j*incx];
                    for( BlasInt i=j+1; i<m; ++i )
                        x[i*incx] -= gamma*Conj(A[j+i*lda]);
                }
            }
        }
    }
}
#ifdef EL_HAVE_QUAD
template void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const Quad* A, BlasInt lda, Quad* x, BlasInt incx );
template void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<Quad>* A, BlasInt lda, Complex<Quad>* x, BlasInt incx );
#endif

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt lda, float* x, BlasInt incx )
{ EL_BLAS(strsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt lda, double* x, BlasInt incx )
{ EL_BLAS(dtrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt lda, scomplex* x, BlasInt incx )
{ EL_BLAS(ctrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt lda, dcomplex* x, BlasInt incx )
{ EL_BLAS(ztrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

// Level 3 BLAS
// ============
template<typename T>
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  T alpha, const T* A, BlasInt lda, const T* B, BlasInt ldb,
  T beta,        T* C, BlasInt ldc )
{
    if( m > 0 && n > 0 && k == 0 && beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*ldc] = 0;
        return;
    }

    // Scale C
    if( beta != T(1) )
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*ldc] *= beta;

    // Naive implementation
    if( transA == 'N' && transB == 'N' )
    {
        // C := alpha A B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt l=0; l<k; ++l )
            {
                const T gamma = alpha*B[l+j*ldb];
                for( BlasInt i=0; i<m; ++i )
                    C[i+j*ldc] += gamma*A[i+l*lda];
            }
        }
    }
    else if( transA == 'N' )
    {
        if( transB == 'T' )
        {
            // C := alpha A B^T + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt l=0; l<k; ++l )
                {
                    const T gamma = alpha*B[j+l*ldb];
                    for( BlasInt i=0; i<m; ++i )
                        C[i+j*ldc] += gamma*A[i+l*lda];
                }
            }
        }
        else
        {
            // C := alpha A B^H + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt l=0; l<k; ++l )
                {
                    const T gamma = alpha*Conj(B[j+l*ldb]);
                    for( BlasInt i=0; i<m; ++i )
                        C[i+j*ldc] += gamma*A[i+l*lda];
                }
            }
        }
    }
    else if( transB == 'N' )
    {
        if( transA == 'T' )
        {
            // C := alpha A^T B + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    T gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                        gamma += A[l+i*lda]*B[l+j*ldb];
                    C[i+j*ldc] += alpha*gamma;
                }
            }
        }
        else
        {
            // C := alpha A^H B + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    T gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                        gamma += Conj(A[l+i*lda])*B[l+j*ldb];
                    C[i+j*ldc] += alpha*gamma;
                }
            }
        }
    }
    else
    {
        if( transA == 'T' && transB == 'T' )
        {
            // C := alpha A^T B^T + C
            for( BlasInt j=0; j<n; ++j )
                for( BlasInt i=0; i<m; ++i )
                    for( BlasInt l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*B[j+l*ldb];
        }
        else if( transA == 'T' )
        {
            // C := alpha A^T B^H + C
            for( BlasInt j=0; j<n; ++j )
                for( BlasInt i=0; i<m; ++i )
                    for( BlasInt l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*Conj(B[j+l*ldb]);
        }
        else if( transB == 'T' )
        {
            // C := alpha A^H B^T + C
            for( BlasInt j=0; j<n; ++j )
                for( BlasInt i=0; i<m; ++i )
                    for( BlasInt l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*B[j+l*ldb];
        }
        else
        {
            // C := alpha A^H B^H + C
            for( BlasInt j=0; j<n; ++j )
                for( BlasInt i=0; i<m; ++i )
                    for( BlasInt l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*Conj(B[j+l*ldb]);
        }
    }
}
template void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  Int alpha, const Int* A, BlasInt lda, const Int* B, BlasInt ldb,
  Int beta,        Int* C, BlasInt ldc );
#ifdef EL_HAVE_QUAD
template void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  Quad alpha, const Quad* A, BlasInt lda, const Quad* B, BlasInt ldb,
  Quad beta,        Quad* C, BlasInt ldc );
template void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda, 
                       const Complex<Quad>* B, BlasInt ldb,
  Complex<Quad> beta,        Complex<Quad>* C, BlasInt ldc );
#endif

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  float alpha, const float* A, BlasInt lda, const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");
      if( transA == 'N' )
      {
          if( lda < Max(m,1) )
              LogicError("lda was too small: lda=",lda,",m=",m);
      }
      else
      {
          if( lda < Max(k,1) )
              LogicError("lda was too small: lda=",lda,",k=",k);
      }

      if( transB == 'N' )
      {
          if( ldb < Max(k,1) )
              LogicError("ldb was too small: ldb=",ldb,",k=",k);
      }
      else
      {
          if( ldb < Max(n,1) )
              LogicError("ldb was too small: ldb=",ldb,",n=",n);
      }

      if( ldc < Max(m,1) )
          LogicError("ldc was too small: ldc=",ldc,",m=",m);
    )
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    EL_BLAS(sgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");

      if( transA == 'N' )
      {
          if( lda < Max(m,1) )
              LogicError("lda was too small: lda=",lda,",m=",m);
      }
      else
      {
          if( lda < Max(k,1) )
              LogicError("lda was too small: lda=",lda,",k=",k);
      }      

      if( transB == 'N' )
      {
          if( ldb < Max(k,1) )
              LogicError("ldb was too small: ldb=",ldb,",k=",k);
      }
      else
      {
          if( ldb < Max(n,1) )
              LogicError("ldb was too small: ldb=",ldb,",n=",n);
      }

      if( ldc < Max(m,1) )
          LogicError("ldc was too small: ldc=",ldc,",m=",m);
    )
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    EL_BLAS(dgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");

      if( transA == 'N' )
      {
          if( lda < Max(m,1) )
              LogicError("lda was too small: lda=",lda,",m=",m);
      }
      else
      {
          if( lda < Max(k,1) )
              LogicError("lda was too small: lda=",lda,",k=",k);
      }      

      if( transB == 'N' )
      {
          if( ldb < Max(k,1) )
              LogicError("ldb was too small: ldb=",ldb,",k=",k);
      }
      else
      {
          if( ldb < Max(n,1) )
              LogicError("ldb was too small: ldb=",ldb,",n=",n);
      }

      if( ldc < Max(m,1) )
          LogicError("ldc was too small: ldc=",ldc,",m=",m);
    )
    EL_BLAS(cgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");

      if( transA == 'N' )
      {
          if( lda < Max(m,1) )
              LogicError("lda was too small: lda=",lda,",m=",m);
      }
      else
      {
          if( lda < Max(k,1) )
              LogicError("lda was too small: lda=",lda,",k=",k);
      }      

      if( transB == 'N' )
      {
          if( ldb < Max(k,1) )
              LogicError("ldb was too small: ldb=",ldb,",k=",k);
      }
      else
      {
          if( ldb < Max(n,1) )
              LogicError("ldb was too small: ldb=",ldb,",n=",n);
      }

      if( ldc < Max(m,1) )
          LogicError("ldc was too small: ldc=",ldc,",m=",m);
    )
    EL_BLAS(zgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

template<typename T>
void Hemm
( char side, char uplo, BlasInt m, BlasInt n, 
  T alpha, const T* A, BlasInt lda, 
           const T* B, BlasInt ldb,
  T beta,        T* C, BlasInt ldc )
{
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*ldc] = 0;
        return;
    }

    // Scale C
    if( beta != T(1) )
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*ldc] *= beta;

    // Naive implementation
    if( side == 'L' && uplo == 'L' )
    {
        // C := alpha tril(A) B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<=i; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*B[l+j*ldb];

        // C := alpha tril(A,-1)^H B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=i+1; l<m; ++l )
                    C[i+j*ldc] += alpha*Conj(A[l+i*lda])*B[l+j*ldb];
    }
    else if( side == 'L' && uplo == 'U' )
    {
        // C := alpha triu(A) B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=i; l<m; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*B[l+j*ldb];

        // C := alpha triu(A,-)^H B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<i; ++l )
                    C[i+j*ldc] += alpha*Conj(A[l+i*lda])*B[l+j*ldb];
    }
    else if( side == 'R' && uplo == 'L' )
    {
        // C := alpha B tril(A) + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=j; l<n; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*A[l+j*lda];

        // C := alpha B tril(A,-1)^H + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<j; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*Conj(A[j+l*lda]);
    }
    else if( side == 'R' && uplo == 'U' )
    {
        // C := alpha B triu(A) + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<=j; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*A[l+j*lda];

        // C := alpha B triu(A,1)^H + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=j+1; l<n; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*Conj(A[j+l*lda]);
    }
    else
        LogicError("Unsuported Hemm option");
}
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  Int alpha, const Int* A, BlasInt lda, 
             const Int* B, BlasInt ldb,
  Int beta,        Int* C, BlasInt ldc );
#ifdef EL_HAVE_QUAD
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  Quad alpha, const Quad* A, BlasInt lda, 
              const Quad* B, BlasInt ldb,
  Quad beta,        Quad* C, BlasInt ldc );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda, 
                       const Complex<Quad>* B, BlasInt ldb,
  Complex<Quad> beta,        Complex<Quad>* C, BlasInt ldc );
#endif

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc )
{
    EL_BLAS(chemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc )
{
    EL_BLAS(zhemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(ssyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dsyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  float beta,           scomplex* C, BlasInt ldc )
{
    EL_BLAS(cher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  double beta,          dcomplex* C, BlasInt ldc )
{
    EL_BLAS(zher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

template<typename T>
void Herk
( char uplo, BlasInt n, BlasInt k, 
  Base<T> alpha, const T* A, BlasInt lda, 
  Base<T> beta,        T* C, BlasInt ldc )
{
    if( uplo == 'L' )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=j; i<n; ++i )
                for( BlasInt l=0; l<k; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*Conj(A[j+l*lda]);
    }
    else
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<=j; ++i )
                for( BlasInt l=0; l<k; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*Conj(A[j+l*lda]);
    }
}
template void Herk
( char uplo, BlasInt n, BlasInt k, 
  Int alpha, const Int* A, BlasInt lda, 
  Int beta,        Int* C, BlasInt ldc );
#ifdef EL_HAVE_QUAD
template void Herk
( char uplo, BlasInt n, BlasInt k, 
  Quad alpha, const Quad* A, BlasInt lda, 
  Quad beta,        Quad* C, BlasInt ldc );
template void Herk
( char uplo, BlasInt n, BlasInt k,
  Quad alpha, const Complex<Quad>* A, BlasInt lda, 
  Quad beta,        Complex<Quad>* C, BlasInt ldc );
#endif

void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda,
  float beta,        float* C, BlasInt ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(ssyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda,
  double beta,        double* C, BlasInt ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dsyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const scomplex* A, BlasInt lda,
  float beta,        scomplex* C, BlasInt ldc )
{ EL_BLAS(cherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const dcomplex* A, BlasInt lda,
  double beta,        dcomplex* C, BlasInt ldc )
{ EL_BLAS(zherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

template<typename T>
void Symm
( char side, char uplo, BlasInt m, BlasInt n, 
  T alpha, const T* A, BlasInt lda, 
           const T* B, BlasInt ldb,
  T beta,        T* C, BlasInt ldc )
{
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*ldc] = 0;
        return;
    }

    // Scale C
    if( beta != T(1) )
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*ldc] *= beta;

    // Naive implementation
    if( side == 'L' && uplo == 'L' )
    {
        // C := alpha tril(A) B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<=i; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*B[l+j*ldb];

        // C := alpha tril(A,-1)^T B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=i+1; l<m; ++l )
                    C[i+j*ldc] += alpha*A[l+i*lda]*B[l+j*ldb];
    }
    else if( side == 'L' && uplo == 'U' )
    {
        // C := alpha triu(A) B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=i; l<m; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*B[l+j*ldb];

        // C := alpha triu(A,-)^T B + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<i; ++l )
                    C[i+j*ldc] += alpha*A[l+i*lda]*B[l+j*ldb];
    }
    else if( side == 'R' && uplo == 'L' )
    {
        // C := alpha B tril(A) + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=j; l<n; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*A[l+j*lda];

        // C := alpha B tril(A,-1)^T + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<j; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*A[j+l*lda];
    }
    else if( side == 'R' && uplo == 'U' )
    {
        // C := alpha B triu(A) + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=0; l<=j; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*A[l+j*lda];

        // C := alpha B triu(A,1)^T + C
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                for( BlasInt l=j+1; l<n; ++l )
                    C[i+j*ldc] += alpha*B[i+l*ldb]*A[j+l*lda];
    }
    else
        LogicError("Unsuported Hemm option");
}
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  Int alpha, const Int* A, BlasInt lda, 
             const Int* B, BlasInt ldb,
  Int beta,        Int* C, BlasInt ldc );
#ifdef EL_HAVE_QUAD
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  Quad alpha, const Quad* A, BlasInt lda, 
              const Quad* B, BlasInt ldb,
  Quad beta,        Quad* C, BlasInt ldc );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda, 
                       const Complex<Quad>* B, BlasInt ldb,
  Complex<Quad> beta,        Complex<Quad>* C, BlasInt ldc );
#endif

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc )
{
    EL_BLAS(csymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc )
{
    EL_BLAS(zsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc )
{
    EL_BLAS(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc )
{
    EL_BLAS(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  scomplex alpha, const scomplex* A, BlasInt lda, const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc )
{
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc )
{
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

template<typename T>
void Syrk
( char uplo, BlasInt n, BlasInt k, 
  T alpha, const T* A, BlasInt lda, 
  T beta,        T* C, BlasInt ldc )
{
    if( uplo == 'L' )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=j; i<n; ++i )
                for( BlasInt l=0; l<k; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*A[j+l*lda];
    }
    else
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<=j; ++i )
                for( BlasInt l=0; l<k; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*A[j+l*lda];
    }
}
template void Syrk
( char uplo, BlasInt n, BlasInt k, 
  Int alpha, const Int* A, BlasInt lda, 
  Int beta,        Int* C, BlasInt ldc );
#ifdef EL_HAVE_QUAD
template void Syrk
( char uplo, BlasInt n, BlasInt k, 
  Quad alpha, const Quad* A, BlasInt lda, 
  Quad beta,        Quad* C, BlasInt ldc );
template void Syrk
( char uplo, BlasInt n, BlasInt k,
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda, 
  Complex<Quad> beta,        Complex<Quad>* C, BlasInt ldc );
#endif

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda,
  float beta,        float* C, BlasInt ldc )
{ EL_BLAS(ssyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda,
  double beta,        double* C, BlasInt ldc )
{ EL_BLAS(dsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  scomplex alpha, const scomplex* A, BlasInt lda,
  scomplex beta,        scomplex* C, BlasInt ldc )
{ EL_BLAS(csyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  dcomplex alpha, const dcomplex* A, BlasInt lda,
  dcomplex beta,        dcomplex* C, BlasInt ldc )
{ EL_BLAS(zsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, float* B, BlasInt ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    EL_BLAS(strmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, double* B, BlasInt ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    EL_BLAS(dtrmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb )
{
    EL_BLAS(ctrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb )
{
    EL_BLAS(ztrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
}

void Trsm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, float* B, BlasInt ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(strsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, double* B, BlasInt ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dtrsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb )
{
    EL_BLAS(ctrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb )
{
    EL_BLAS(ztrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &lda, B, &ldb );
} 

} // namespace blas
} // namespace El
