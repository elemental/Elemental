/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental.hpp"

extern "C" {

//------------------------------------------------------------------------//
// Level 1 BLAS                                                           //
//------------------------------------------------------------------------//
void BLAS(saxpy)
( const int* n, const float* alpha, const float* x, const int* incx,
                                          float* y, const int* incy );
void BLAS(daxpy)
( const int* n, const double* alpha, const double* x, const int* incx,
                                           double* y, const int* incy );
void BLAS(caxpy)
( const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* x, const int* incx,
        elem::scomplex* y, const int* incy );
void BLAS(zaxpy)
( const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* x, const int* incx,
        elem::dcomplex* y, const int* incy );

void BLAS(scopy)
( const int* n, const float* x, const int* incx,
                      float* y, const int* incy );
void BLAS(dcopy)
( const int* n, const double* x, const int* incx,
                      double* y, const int* incy );
void BLAS(ccopy)
( const int* n, const elem::scomplex* x, const int* incx,
                      elem::scomplex* y, const int* incy );
void BLAS(zcopy)
( const int* n, const elem::dcomplex* x, const int* incx,
                      elem::dcomplex* y, const int* incy );

float BLAS(sdot)
( const int* n, const float* x, const int* incx,
                const float* y, const int* incy );
double BLAS(ddot)
( const int* n, const double* x, const int* incx,
                const double* y, const int* incy );
// To avoid the compatibility issue, we simply handroll our own complex dots
float BLAS(snrm2)
( const int* n, const float* x, const int* incx );
double BLAS(dnrm2)
( const int* n, const double* x, const int* incx );
float BLAS(scnrm2)
( const int* n, const elem::scomplex* x, const int* incx );
double BLAS(dznrm2)
( const int* n, const elem::dcomplex* x, const int* incx );

void BLAS(sscal)
( const int* n, const float* alpha, float* x, const int* incx );
void BLAS(dscal)
( const int* n, const double* alpha, double* x, const int* incx );
void BLAS(cscal)
( const int* n, const elem::scomplex* alpha, elem::scomplex* x,
  const int* incx );
void BLAS(zscal)
( const int* n, const elem::dcomplex* alpha, elem::dcomplex* x,
  const int* incx );

//------------------------------------------------------------------------//
// Level 2 BLAS                                                           //
//------------------------------------------------------------------------//
void BLAS(sgemv)
( const char* trans, const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                      const float* x, const int* incx,
  const float* beta,        float* y, const int* incy );
void BLAS(dgemv)
( const char* trans, const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                       const double* x, const int* incx,
  const double* beta,        double* y, const int* incy );
void BLAS(cgemv)
( const char* trans, const int* m, const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* x, const int* incx,
  const elem::scomplex* beta,
        elem::scomplex* y, const int* incy );
void BLAS(zgemv)
( const char* trans, const int* m, const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* x, const int* incx,
  const elem::dcomplex* beta,
        elem::dcomplex* y, const int* incy );

void BLAS(sger)
( const int* m, const int* n,
  const float* alpha, const float* x, const int* incx,
                      const float* y, const int* incy,
                            float* A, const int* lda  );
void BLAS(dger)
( const int* m, const int* n,
  const double* alpha, const double* x, const int* incx,
                       const double* y, const int* incy,
                             double* A, const int* lda  );
void BLAS(cgerc)
( const int* m, const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* x, const int* incx,
  const elem::scomplex* y, const int* incy,
        elem::scomplex* A, const int* lda  );
void BLAS(zgerc)
( const int* m, const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* x, const int* incx,
  const elem::dcomplex* y, const int* incy,
        elem::dcomplex* A, const int* lda  );

void BLAS(cgeru)
( const int* m, const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* x, const int* incx,
  const elem::scomplex* y, const int* incy,
        elem::scomplex* A, const int* lda  );
void BLAS(zgeru)
( const int* m, const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* x, const int* incx,
  const elem::dcomplex* y, const int* incy,
        elem::dcomplex* A, const int* lda  );

void BLAS(chemv)
( const char* uplo, const int* m,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* x, const int* incx,
  const elem::scomplex* beta,
        elem::scomplex* y, const int* incy );
void BLAS(zhemv)
( const char* uplo, const int* m,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* x, const int* incx,
  const elem::dcomplex* beta,
        elem::dcomplex* y, const int* incy );

void BLAS(cher)
( const char* uplo, const int* m,
  const elem::scomplex* alpha,
  const elem::scomplex* x, const int* incx,
        elem::scomplex* A, const int* lda  );
void BLAS(zher)
( const char* uplo, const int* m,
  const elem::dcomplex* alpha,
  const elem::dcomplex* x, const int* incx,
        elem::dcomplex* A, const int* lda  );

void BLAS(cher2)
( const char* uplo, const int* m,
  const elem::scomplex* alpha,
  const elem::scomplex* x, const int* incx,
  const elem::scomplex* y, const int* incy,
        elem::scomplex* A, const int* lda  );
void BLAS(zher2)
( const char* uplo, const int* m,
  const elem::dcomplex* alpha,
  const elem::dcomplex* x, const int* incx,
  const elem::dcomplex* y, const int* incy,
        elem::dcomplex* A, const int* lda  );

void BLAS(ssymv)
( const char* uplo, const int* m,
  const float* alpha, const float* A, const int* lda,
                      const float* x, const int* incx,
  const float* beta,        float* y, const int* incy );
void BLAS(dsymv)
( const char* uplo, const int* m,
  const double* alpha, const double* A, const int* lda,
                       const double* x, const int* incx,
  const double* beta,        double* y, const int* incy );
// 'csymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void LAPACK(csymv)
( const char* uplo, const int* m,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* x, const int* incx,
  const elem::scomplex* beta,
        elem::scomplex* y, const int* incy );
// 'zsymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void LAPACK(zsymv)
( const char* uplo, const int* m,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* x, const int* incx,
  const elem::dcomplex* beta,
        elem::dcomplex* y, const int* incy );

void BLAS(ssyr)
( const char* uplo, const int* m,
  const float* alpha, const float* x, const int* incx,
                            float* A, const int* lda  );
void BLAS(dsyr)
( const char* uplo, const int* m,
  const double* alpha, const double* x, const int* incx,
                             double* A, const int* lda  );
// 'csyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void LAPACK(csyr)
( const char* uplo, const int* m,
  const elem::scomplex* alpha,
  const elem::scomplex* x, const int* incx,
        elem::scomplex* A, const int* lda  );
// 'zsyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void LAPACK(zsyr)
( const char* uplo, const int* m,
  const elem::dcomplex* alpha,
  const elem::dcomplex* x, const int* incx,
        elem::dcomplex* A, const int* lda  );

void BLAS(ssyr2)
( const char* uplo, const int* m,
  const float* alpha, const float* x, const int* incx,
                      const float* y, const int* incy,
                            float* A, const int* lda  );
void BLAS(dsyr2)
( const char* uplo, const int* m,
  const double* alpha, const double* x, const int* incx,
                       const double* y, const int* incy,
                             double* A, const int* lda  );

void BLAS(strmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const float* A, const int* lda, float* x, const int* incx );
void BLAS(dtrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const double* A, const int* lda, double* x, const int* incx );
void BLAS(ctrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elem::scomplex* A, const int* lda,
        elem::scomplex* x, const int* incx );
void BLAS(ztrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elem::dcomplex* A, const int* lda,
        elem::dcomplex* x, const int* incx );

void BLAS(strsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const float* A, const int* lda, float* x, const int* incx );
void BLAS(dtrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const double* A, const int* lda, double* x, const int* incx );
void BLAS(ctrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elem::scomplex* A, const int* lda,
        elem::scomplex* x, const int* incx );
void BLAS(ztrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elem::dcomplex* A, const int* lda,
        elem::dcomplex* x, const int* incx );

//------------------------------------------------------------------------//
// Level 3 BLAS                                                           //
//------------------------------------------------------------------------//
void BLAS(sgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );
void BLAS(dgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );
void BLAS(cgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* B, const int* ldb,
  const elem::scomplex* beta,
        elem::scomplex* C, const int* ldc );
void BLAS(zgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* B, const int* ldb,
  const elem::dcomplex* beta,
        elem::dcomplex* C, const int* ldc );

void BLAS(chemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* B, const int* ldb,
  const elem::scomplex* beta,
        elem::scomplex* C, const int* ldc );
void BLAS(zhemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* B, const int* ldb,
  const elem::dcomplex* beta,
        elem::dcomplex* C, const int* ldc );

void BLAS(cher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* B, const int* ldb,
  const elem::scomplex* beta,
        elem::scomplex* C, const int* ldc );
void BLAS(zher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* B, const int* ldb,
  const elem::dcomplex* beta,
        elem::dcomplex* C, const int* ldc );

void BLAS(cherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* beta,
        elem::scomplex* C, const int* ldc );
void BLAS(zherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* beta,
        elem::dcomplex* C, const int* ldc );

void BLAS(ssymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );
void BLAS(dsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );
void BLAS(csymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* B, const int* ldb,
  const elem::scomplex* beta,
        elem::scomplex* C, const int* ldc );
void BLAS(zsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* B, const int* ldb,
  const elem::dcomplex* beta,
        elem::dcomplex* C, const int* ldc );

void BLAS(ssyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );
void BLAS(dsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );
void BLAS(csyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* B, const int* ldb,
  const elem::scomplex* beta,
        elem::scomplex* C, const int* ldc );
void BLAS(zsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* B, const int* ldb,
  const elem::dcomplex* beta,
        elem::dcomplex* C, const int* ldc );

void BLAS(ssyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
  const float* beta,        float* C, const int* ldc );
void BLAS(dsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
  const double* beta,        double* C, const int* ldc );
void BLAS(csyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
  const elem::scomplex* beta,
        elem::scomplex* C, const int* ldc );
void BLAS(zsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
  const elem::dcomplex* beta,
        elem::dcomplex* C, const int* ldc );

void BLAS(strmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                            float* B, const int* ldb );
void BLAS(dtrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                             double* B, const int* ldb );
void BLAS(ctrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
        elem::scomplex* B, const int* ldb );
void BLAS(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
        elem::dcomplex* B, const int* ldb );

void BLAS(strsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                            float* B, const int* ldb );
void BLAS(dtrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                             double* B, const int* ldb );
void BLAS(ctrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const elem::scomplex* alpha,
  const elem::scomplex* A, const int* lda,
        elem::scomplex* B, const int* ldb );
void BLAS(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const elem::dcomplex* alpha,
  const elem::dcomplex* A, const int* lda,
        elem::dcomplex* B, const int* ldb );

} // extern "C"

namespace elem {
namespace blas {

//----------------------------------------------------------------------------//
// Level 1 BLAS                                                               //
//----------------------------------------------------------------------------//
void Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy )
{ BLAS(saxpy)( &n, &alpha, x, &incx, y, &incy ); }

void Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy )
{ BLAS(daxpy)( &n, &alpha, x, &incx, y, &incy ); }

void Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy )
{ BLAS(caxpy)( &n, &alpha, x, &incx, y, &incy ); }

void Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy )
{ BLAS(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }

void Copy( int n, const float* x, int incx, float* y, int incy )
{ BLAS(scopy)( &n, x, &incx, y, &incy ); }

void Copy( int n, const double* x, int incx, double* y, int incy )
{ BLAS(dcopy)( &n, x, &incx, y, &incy ); }

void Copy( int n, const scomplex* x, int incx, scomplex* y, int incy )
{ BLAS(ccopy)( &n, x, &incx, y, &incy ); }

void Copy( int n, const dcomplex* x, int incx, dcomplex* y, int incy )
{ BLAS(zcopy)( &n, x, &incx, y, &incy ); }

float Dot( int n, const float* x, int incx, const float* y, int incy )
{ return BLAS(sdot)( &n, x, &incx, y, &incy ); }

double Dot( int n, const double* x, int incx, const double* y, int incy )
{ return BLAS(ddot)( &n, x, &incx, y, &incy ); }

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
{ return BLAS(sdot)( &n, x, &incx, y, &incy ); }

double Dotc( int n, const double* x, int incx, const double* y, int incy )
{ return BLAS(ddot)( &n, x, &incx, y, &incy ); }

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
{ return BLAS(sdot)( &n, x, &incx, y, &incy ); }

double Dotu( int n, const double* x, int incx, const double* y, int incy )
{ return BLAS(ddot)( &n, x, &incx, y, &incy ); }

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
{ return BLAS(snrm2)( &n, x, &incx ); }

double Nrm2( int n, const double* x, int incx )
{ return BLAS(dnrm2)( &n, x, &incx ); }

float Nrm2( int n, const scomplex* x, int incx )
{ return BLAS(scnrm2)( &n, x, &incx ); }

double Nrm2( int n, const dcomplex* x, int incx )
{ return BLAS(dznrm2)( &n, x, &incx ); }

void Scal( int n, float alpha, float* x, int incx )
{ BLAS(sscal)( &n, &alpha, x, &incx ); }

void Scal( int n, double alpha, double* x, int incx )
{ BLAS(dscal)( &n, &alpha, x, &incx ); }

void Scal( int n, scomplex alpha, scomplex* x, int incx )
{ BLAS(cscal)( &n, &alpha, x, &incx ); }

void Scal( int n, dcomplex alpha, dcomplex* x, int incx )
{ BLAS(zscal)( &n, &alpha, x, &incx ); }

//----------------------------------------------------------------------------//
// Level 2 BLAS                                                               //
//----------------------------------------------------------------------------//
void Gemv
( char trans, int m, int n,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, int m, int n,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ BLAS(cgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ BLAS(zgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Ger
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda  )
{ BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Ger
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Gerc
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Geru
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

void Hemv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ BLAS(chemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ BLAS(zhemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda )
{ BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{ BLAS(cher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{ BLAS(zher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Her2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ BLAS(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Her2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ BLAS(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Symv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ BLAS(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ BLAS(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    LAPACK(csymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Symv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    LAPACK(zsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

void Syr
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda  )
{ BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Syr
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

void Syr
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    LAPACK(csyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

void Syr
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    LAPACK(zsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

void Syr2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

void Syr2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

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
    BLAS(csyr2k)
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
    BLAS(zsyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

void Trmv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ BLAS(strmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ BLAS(dtrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ BLAS(ctrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ BLAS(ztrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ BLAS(strsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ BLAS(dtrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ BLAS(ctrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ BLAS(ztrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

//----------------------------------------------------------------------------//
// Level 3 BLAS                                                               //
//----------------------------------------------------------------------------//
void Gemm
( char transA, char transB, int m, int n, int k, 
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    BLAS(sgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
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
    BLAS(dgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB, int m, int n, int k, 
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(cgemm)( &transA, &transB, &m, &n, &k,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Gemm
( char transA, char transB, int m, int n, int k, 
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zgemm)( &transA, &transB, &m, &n, &k,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    BLAS(ssymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    BLAS(dsymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(chemm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Hemm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zhemm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(ssyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(dsyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(cher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(ssyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    BLAS(dsyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

void Herk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ BLAS(cherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Herk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ BLAS(zherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Symm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    BLAS(ssymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    BLAS(dsymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(csymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Symm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zsymm)( &side, &uplo, &m, &n,
                 &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    BLAS(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    BLAS(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    BLAS(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syr2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    BLAS(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

void Syrk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{ BLAS(ssyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{ BLAS(dsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ BLAS(csyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Syrk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ BLAS(zsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    BLAS(strmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    BLAS(dtrmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    BLAS(ctrmm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    BLAS(ztrmm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
}

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(strsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    BLAS(dtrsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    BLAS(ctrsm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    BLAS(ztrsm)( &side, &uplo, &trans, &unit, &m, &n,
                 &alpha, A, &lda, B, &ldb );
} 

} // namespace blas
} // namespace elem
