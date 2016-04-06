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

// Level 1 BLAS
// ============
void EL_BLAS(saxpy)
( const BlasInt* n, const float* alpha,
  const float* x, const BlasInt* incx,
  float* y, const BlasInt* incy );
void EL_BLAS(daxpy)
( const BlasInt* n, const double* alpha,
  const double* x, const BlasInt* incx,
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
( const BlasInt* n,
  const float* x, const BlasInt* incx,
        float* y, const BlasInt* incy );
void EL_BLAS(dcopy)
( const BlasInt* n,
  const double* x, const BlasInt* incx,
        double* y, const BlasInt* incy );
void EL_BLAS(ccopy)
( const BlasInt* n,
  const scomplex* x, const BlasInt* incx,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zcopy)
( const BlasInt* n,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* y, const BlasInt* incy );

float EL_BLAS(sdot)
( const BlasInt* n,
  const float* x, const BlasInt* incx,
  const float* y, const BlasInt* incy );
double EL_BLAS(ddot)
( const BlasInt* n,
  const double* x, const BlasInt* incx,
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
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* x, const BlasInt* incx,
  const float* beta,        float* y, const BlasInt* incy );
void EL_BLAS(dgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* x, const BlasInt* incx,
  const double* beta,        double* y, const BlasInt* incy );
void EL_BLAS(cgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(sger)
( const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* x, const BlasInt* incx,
                      const float* y, const BlasInt* incy,
                            float* A, const BlasInt* ALDim  );
void EL_BLAS(dger)
( const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* x, const BlasInt* incx,
                       const double* y, const BlasInt* incy,
                             double* A, const BlasInt* ALDim  );
void EL_BLAS(cgerc)
( const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* ALDim  );
void EL_BLAS(zgerc)
( const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* ALDim  );

void EL_BLAS(cgeru)
( const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* ALDim );
void EL_BLAS(zgeru)
( const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* ALDim );

void EL_BLAS(chemv)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zhemv)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(cher)
( const char* uplo, const BlasInt* m,
  const float* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* A, const BlasInt* ALDim );
void EL_BLAS(zher)
( const char* uplo, const BlasInt* m,
  const double* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* A, const BlasInt* ALDim );

void EL_BLAS(cher2)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* ALDim );
void EL_BLAS(zher2)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* ALDim );

void EL_BLAS(ssymv)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* x, const BlasInt* incx,
  const float* beta,        float* y, const BlasInt* incy );
void EL_BLAS(dsymv)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* x, const BlasInt* incx,
  const double* beta,        double* y, const BlasInt* incy );
// 'csymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csymv)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
// 'zsymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsymv)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(ssyr)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* x, const BlasInt* incx,
                            float* A, const BlasInt* ALDim  );
void EL_BLAS(dsyr)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* x, const BlasInt* incx,
                             double* A, const BlasInt* ALDim  );
// 'csyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csyr)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* A, const BlasInt* ALDim  );
// 'zsyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsyr)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* A, const BlasInt* ALDim  );

void EL_BLAS(ssyr2)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* x, const BlasInt* incx,
                      const float* y, const BlasInt* incy,
                            float* A, const BlasInt* ALDim  );
void EL_BLAS(dsyr2)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* x, const BlasInt* incx,
                       const double* y, const BlasInt* incy,
                             double* A, const BlasInt* ALDim  );

void EL_BLAS(strmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const float* A, const BlasInt* ALDim, float* x, const BlasInt* incx );
void EL_BLAS(dtrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const double* A, const BlasInt* ALDim, double* x, const BlasInt* incx );
void EL_BLAS(ctrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const scomplex* A, const BlasInt* ALDim,
        scomplex* x, const BlasInt* incx );
void EL_BLAS(ztrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const dcomplex* A, const BlasInt* ALDim,
        dcomplex* x, const BlasInt* incx );

void EL_BLAS(strsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const float* A, const BlasInt* ALDim, float* x, const BlasInt* incx );
void EL_BLAS(dtrsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const double* A, const BlasInt* ALDim, double* x, const BlasInt* incx );
void EL_BLAS(ctrsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const scomplex* A, const BlasInt* ALDim,
        scomplex* x, const BlasInt* incx );
void EL_BLAS(ztrsv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const dcomplex* A, const BlasInt* ALDim,
        dcomplex* x, const BlasInt* incx );

// Level 3 BLAS 
// ============
void EL_BLAS(sgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* B, const BlasInt* BLDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* B, const BlasInt* BLDim,
  const double* beta,        double* C, const BlasInt* CLDim );
void EL_BLAS(cgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(chemm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zhemm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(cher2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const float* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zher2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const double* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(cherk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const float* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zherk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const double* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(ssymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* B, const BlasInt* BLDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dsymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* B, const BlasInt* BLDim,
  const double* beta,        double* C, const BlasInt* CLDim );
void EL_BLAS(csymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zsymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(ssyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* B, const BlasInt* BLDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dsyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* B, const BlasInt* BLDim,
  const double* beta,        double* C, const BlasInt* CLDim );
void EL_BLAS(csyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zsyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(ssyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* ALDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dsyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* ALDim,
  const double* beta,        double* C, const BlasInt* CLDim );
void EL_BLAS(csyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zsyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(strmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* ALDim,
                            float* B, const BlasInt* BLDim );
void EL_BLAS(dtrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                             double* B, const BlasInt* BLDim );
void EL_BLAS(ctrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
        scomplex* B, const BlasInt* BLDim );
void EL_BLAS(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
        dcomplex* B, const BlasInt* BLDim );

void EL_BLAS(strsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* ALDim,
                            float* B, const BlasInt* BLDim );
void EL_BLAS(dtrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                             double* B, const BlasInt* BLDim );
void EL_BLAS(ctrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
        scomplex* B, const BlasInt* BLDim );
void EL_BLAS(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
        dcomplex* B, const BlasInt* BLDim );

} // extern "C"

namespace El {
namespace blas {

// Level 1 BLAS
// ============
template<typename T>
void Axpy
( BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma;
    for( BlasInt i=0; i<n; ++i )
    {
        gamma = alpha*x[i*incx];
        y[i*incy] += gamma;
    }
}
template void Axpy
( BlasInt n,
  const Int& alpha,
  const Int* x, BlasInt incx,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Axpy
( BlasInt n,
  const DoubleDouble& alpha, 
  const DoubleDouble* x, BlasInt incx,
        DoubleDouble* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const QuadDouble& alpha, 
  const QuadDouble* x, BlasInt incx,
        QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Axpy
( BlasInt n,
  const Quad& alpha, 
  const Quad* x, BlasInt incx,
        Quad* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const Complex<Quad>& alpha, 
  const Complex<Quad>* x, BlasInt incx,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Axpy
( BlasInt n,
  const BigInt& alpha, 
  const BigInt* x, BlasInt incx,
        BigInt* y, BlasInt incy );
template void Axpy
( BlasInt n,
  const BigFloat& alpha, 
  const BigFloat* x, BlasInt incx,
        BigFloat* y, BlasInt incy );
#endif

void Axpy
( BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
        float* y, BlasInt incy )
{ EL_BLAS(saxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
        double* y, BlasInt incy )
{ EL_BLAS(daxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
        scomplex* y, BlasInt incy )
{ EL_BLAS(caxpy)( &n, &alpha, x, &incx, y, &incy ); }
void Axpy
( BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }

template<typename T>
void Copy
( BlasInt n,
  const T* x, BlasInt incx,
        T* y, BlasInt incy )
{
    for( BlasInt i=0; i<n; ++i )
        y[i*incy] = x[i*incx];
}
template void Copy
( BlasInt n,
  const Int* x, BlasInt incx,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Copy
( BlasInt n,
  const DoubleDouble* x, BlasInt incx, 
        DoubleDouble* y, BlasInt incy );
template void Copy
( BlasInt n,
  const QuadDouble* x, BlasInt incx, 
        QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Copy
( BlasInt n,
  const Quad* x, BlasInt incx, 
        Quad* y, BlasInt incy );
template void Copy
( BlasInt n,
  const Complex<Quad>* x, BlasInt incx, 
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Copy
( BlasInt n,
  const BigInt* x, BlasInt incx, 
        BigInt* y, BlasInt incy );
template void Copy
( BlasInt n,
  const BigFloat* x, BlasInt incx, 
        BigFloat* y, BlasInt incy );
#endif

void Copy
( BlasInt n,
  const float* x, BlasInt incx,
        float* y, BlasInt incy )
{ EL_BLAS(scopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n,
  const double* x, BlasInt incx,
        double* y, BlasInt incy )
{ EL_BLAS(dcopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n,
  const scomplex* x, BlasInt incx,
        scomplex* y, BlasInt incy )
{ EL_BLAS(ccopy)( &n, x, &incx, y, &incy ); }
void Copy
( BlasInt n,
  const dcomplex* x, BlasInt incx,
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zcopy)( &n, x, &incx, y, &incy ); }

template<typename T>
T Dot
( BlasInt n,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma;
    T alpha = 0;
    for( BlasInt i=0; i<n; ++i )
    {
        Conj( x[i*incx], gamma );
        gamma *= y[i*incy];
        alpha += gamma;
    }
    return alpha;
}
template Int Dot
( BlasInt n,
  const Int* x, BlasInt incx,
  const Int* y, BlasInt incy );
template float Dot
( BlasInt n,
  const float* x, BlasInt incx,
  const float* y, BlasInt incy );
template scomplex Dot
( BlasInt n,
  const scomplex* x, BlasInt incx,
  const scomplex* y, BlasInt incy );
template dcomplex Dot
( BlasInt n,
  const dcomplex* x, BlasInt incx,
  const dcomplex* y, BlasInt incy );
#ifdef EL_HAVE_QD
template DoubleDouble Dot
( BlasInt n,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy );
template QuadDouble Dot
( BlasInt n,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template Quad Dot
( BlasInt n,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy );
template Complex<Quad> Dot
( BlasInt n,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template BigInt Dot
( BlasInt n,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy );
template BigFloat Dot
( BlasInt n,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy );
#endif

// NOTE: I am under the impression that it is generally unsafe to return 
//       anything except a double-precision float to C from Fortran
double Dot
( BlasInt n,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy )
{ return EL_BLAS(ddot)( &n, x, &incx, y, &incy ); }

template<typename T>
T Dotu
( BlasInt n,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma;
    T alpha = 0;
    for( BlasInt i=0; i<n; ++i )
    {
        gamma = x[i*incx];
        gamma *= y[i*incy];
        alpha += gamma;
    }
    return alpha;
}
template Int Dotu
( BlasInt n,
  const Int* x, BlasInt incx,
  const Int* y, BlasInt incy );
template float Dotu
( BlasInt n,
  const float* x, BlasInt incx,
  const float* y, BlasInt incy );
template scomplex Dotu
( BlasInt n,
  const scomplex* x, BlasInt incx,
  const scomplex* y, BlasInt incy );
template dcomplex Dotu
( BlasInt n,
  const dcomplex* x, BlasInt incx,
  const dcomplex* y, BlasInt incy );
#ifdef EL_HAVE_QD
template DoubleDouble Dotu
( BlasInt n,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy );
template QuadDouble Dotu
( BlasInt n,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template Quad Dotu
( BlasInt n,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy );
template Complex<Quad> Dotu
( BlasInt n,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template BigInt Dotu
( BlasInt n,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy );
template BigFloat Dotu
( BlasInt n,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy );
#endif

double Dotu
( BlasInt n,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy )
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
#ifdef EL_HAVE_QD
template DoubleDouble Nrm2( BlasInt n, const DoubleDouble* x, BlasInt incx );
template QuadDouble Nrm2( BlasInt n, const QuadDouble* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template Quad Nrm2( BlasInt n, const Quad* x, BlasInt incx );
template Quad Nrm2( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BigFloat Nrm2( BlasInt n, const BigFloat* x, BlasInt incx );
#endif

double Nrm2( BlasInt n, const double* x, BlasInt incx )
{ return EL_BLAS(dnrm2)( &n, x, &incx ); }
double Nrm2( BlasInt n, const dcomplex* x, BlasInt incx )
{ return EL_BLAS(dznrm2)( &n, x, &incx ); }

template<typename F>
BlasInt MaxInd( BlasInt n, const F* x, BlasInt incx )
{
    typedef Base<F> Real;
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    //       (A copy elision is assumed in the return value of Abs)
    Real absVal;
    Real maxAbsVal = -1;
    BlasInt maxAbsInd = -1;
    for( BlasInt i=0; i<n; ++i ) 
    {
        absVal = Abs(x[i*incx]);
        if( absVal > maxAbsVal )
        {
            maxAbsVal = absVal;
            maxAbsInd = i;
        }
    } 
    return maxAbsInd;
}
template BlasInt MaxInd( BlasInt n, const Int* x, BlasInt incx );
#ifdef EL_HAVE_QD
template BlasInt MaxInd( BlasInt n, const DoubleDouble* x, BlasInt incx );
template BlasInt MaxInd( BlasInt n, const QuadDouble* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template BlasInt MaxInd( BlasInt n, const Quad* x, BlasInt incx );
template BlasInt MaxInd( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BlasInt MaxInd( BlasInt n, const BigInt* x, BlasInt incx );
template BlasInt MaxInd( BlasInt n, const BigFloat* x, BlasInt incx );
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

template<typename Real>
Real Givens
( const Real& phi,
  const Real& gamma,
  Real* c,
  Real* s )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    Real zero(0), one(1);

    Real phiAbs = Abs(phi);
    if( phiAbs == zero )
    {
        *c = zero;
        *s = one;
        return gamma;
    }
    else
    {
        Real gammaAbs = Abs(gamma);
        Real scale = phiAbs + gammaAbs;
        Real norm = scale*lapack::SafeNorm(phi/scale,gamma/scale);
        *c = phiAbs/norm;
        *s = (phi/phiAbs)*gamma/norm;
        return *c*phi + *s*gamma;
    }
}
template<typename Real>
Complex<Real> Givens
( const Complex<Real>& phi,
  const Complex<Real>& gamma,
  Real* c,
  Complex<Real>* s )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    Real zero(0), one(1);

    Real phiAbs = Abs(phi);
    if( phiAbs == zero )
    {
        *c = zero;
        *s = one;
        return gamma;
    }
    else
    {
        Real gammaAbs = Abs(gamma);
        Real scale = phiAbs + gammaAbs;
        Real norm = scale*lapack::SafeNorm(Abs(phi/scale),Abs(gamma/scale));
        *c = phiAbs/norm;
        *s = (phi/phiAbs)*Conj(gamma)/norm;
        return *c*phi + *s*gamma;
    }
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
( const Quad& phi,
  const Quad& gamma,
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

float Givens
( const float& phi,
  const float& gamma,
  float* c,
  float* s )
{
    float phiCopy=phi, gammaCopy=gamma;
    EL_BLAS(srotg)( &phiCopy, &gammaCopy, c, s );
    return phi;
}
double Givens
( const double& phi,
  const double& gamma,
  double* c,
  double* s )
{
    double phiCopy=phi, gammaCopy=gamma;
    EL_BLAS(drotg)( &phiCopy, &gammaCopy, c, s );
    return phi;
}
scomplex Givens
( const scomplex& phi,
  const scomplex& gamma,
  float* c,
  scomplex* s )
{
    scomplex phiCopy=phi, gammaCopy=gamma;
    EL_BLAS(crotg)( &phiCopy, &gammaCopy, c, s );
    return phi;
}
dcomplex Givens
( const dcomplex& phi,
  const dcomplex& gamma,
  double* c,
  dcomplex* s )
{
    dcomplex phiCopy=phi, gammaCopy=gamma;
    EL_BLAS(zrotg)( &phiCopy, &gammaCopy, c, s );
    return phi;
}

template<typename F>
void Rot
( BlasInt n,
  F* x, BlasInt incx,
  F* y, BlasInt incy,
  const Base<F>& c,
  const F& s )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    F gamma, delta;
    for( BlasInt i=0; i<n; ++i )    
    {
        //gamma = c*x[i*incx] + s*y[i*incy];
        gamma = c;
        gamma *= x[i*incx];
        delta = s;
        delta *= y[i*incy];
        gamma += delta;

        //y[i*incy] = -Conj(s)*x[i*incx] + c*y[i*incy];
        y[i*incy] *= c;
        Conj( s, delta );
        delta *= x[i*incx];
        y[i*incy] -= delta;

        x[i*incx] = gamma;
    }
}
#ifdef EL_HAVE_QD
template void Rot
( BlasInt n,
  DoubleDouble* x, BlasInt incx,
  DoubleDouble* y, BlasInt incy,
  const DoubleDouble& c,
  const DoubleDouble& s );
template void Rot
( BlasInt n,
  QuadDouble* x, BlasInt incx,
  QuadDouble* y, BlasInt incy,
  const QuadDouble& c,
  const QuadDouble& s );
#endif
#ifdef EL_HAVE_QUAD
template void Rot
( BlasInt n,
  Quad* x, BlasInt incx,
  Quad* y, BlasInt incy,
  const Quad& c,
  const Quad& s );
template void Rot
( BlasInt n,
  Complex<Quad>* x, BlasInt incx,
  Complex<Quad>* y, BlasInt incy,
  const Quad& c,
  const Complex<Quad>& s );
#endif
#ifdef EL_HAVE_MPC
template void Rot
( BlasInt n,
  BigFloat* x, BlasInt incx,
  BigFloat* y, BlasInt incy,
  const BigFloat& c,
  const BigFloat& s );
#endif

void Rot
( BlasInt n,
  float* x, BlasInt incx, 
  float* y, BlasInt incy,
  const float& c,
  const float& s )
{ EL_BLAS(srot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n,
  double* x, BlasInt incx, 
  double* y, BlasInt incy,
  const double& c,
  const double& s )
{ EL_BLAS(drot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n,
  scomplex* x, BlasInt incx, 
  scomplex* y, BlasInt incy,
  const float& c,
  const scomplex& s )
{ EL_BLAS(crot)( &n, x, &incx, y, &incy, &c, &s ); }
void Rot
( BlasInt n,
  dcomplex* x, BlasInt incx, 
  dcomplex* y, BlasInt incy,
  const double& c,
  const dcomplex& s )
{ EL_BLAS(zrot)( &n, x, &incx, y, &incy, &c, &s ); }

template<typename T>
void Scal( BlasInt n, const T& alpha, T* x, BlasInt incx )
{
    for( BlasInt j=0; j<n; ++j )
        x[j*incx] *= alpha;
}
template void Scal
( BlasInt n,
  const Int& alpha,
  Int* x, BlasInt incx );
#ifdef EL_HAVE_QD
template void Scal
( BlasInt n,
  const DoubleDouble& alpha,
  DoubleDouble* x, BlasInt incx );
template void Scal
( BlasInt n,
  const QuadDouble& alpha,
  QuadDouble* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template void Scal
( BlasInt n,
  const Quad& alpha,
  Quad* x, BlasInt incx );
template void Scal
( BlasInt n,
  const Complex<Quad>& alpha,
  Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template void Scal
( BlasInt n,
  const BigInt& alpha,
  BigInt* x, BlasInt incx );
template void Scal
( BlasInt n,
  const BigFloat& alpha,
  BigFloat* x, BlasInt incx );
#endif

template<typename T>
void Scal( BlasInt n, const T& alpha, Complex<T>* x, BlasInt incx )
{
    for( BlasInt j=0; j<n; ++j )
        x[j*incx] *= alpha;
}
template void Scal
( BlasInt n,
  const float& alpha,
  Complex<float>* x, BlasInt incx );
template void Scal
( BlasInt n,
  const double& alpha,
  Complex<double>* x, BlasInt incx );
#ifdef EL_HAVE_QUAD
template void Scal
( BlasInt n,
  const Quad& alpha,
  Complex<Quad>* x, BlasInt incx );
#endif

void Scal( BlasInt n, const float& alpha, float* x, BlasInt incx )
{ EL_BLAS(sscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, const double& alpha, double* x, BlasInt incx )
{ EL_BLAS(dscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, const scomplex& alpha, scomplex* x, BlasInt incx )
{ EL_BLAS(cscal)( &n, &alpha, x, &incx ); }
void Scal( BlasInt n, const dcomplex& alpha, dcomplex* x, BlasInt incx )
{ EL_BLAS(zscal)( &n, &alpha, x, &incx ); }

// NOTE: 'nrm1' is not the official name but is consistent with 'nrm2'
template<typename F>
Base<F> Nrm1( BlasInt n, const F* x, BlasInt incx )
{
    // TODO: Avoid temporaries since constructing BigInt/BigFloat involves
    //       a memory allocation
    Base<F> sum=0;
    for( BlasInt i=0; i<n; ++i )
        sum += Abs(x[i*incx]);
    return sum;
}
template float Nrm1( BlasInt n, const float* x, BlasInt incx );
template float Nrm1( BlasInt n, const scomplex* x, BlasInt incx );
#ifdef EL_HAVE_QD
template DoubleDouble Nrm1( BlasInt n, const DoubleDouble* x, BlasInt incx );
template QuadDouble Nrm1( BlasInt n, const QuadDouble* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template Quad Nrm1( BlasInt n, const Quad* x, BlasInt incx );
template Quad Nrm1( BlasInt n, const Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BigInt Nrm1( BlasInt n, const BigInt* x, BlasInt incx );
template BigFloat Nrm1( BlasInt n, const BigFloat* x, BlasInt incx );
#endif

double Nrm1( BlasInt n, const double* x, BlasInt incx )
{ return EL_BLAS(dasum)( &n, x, &incx ); }
double Nrm1( BlasInt n, const dcomplex* x, BlasInt incx )
{ return EL_LAPACK(dzsum1)( &n, x, &incx ); }

template<typename T>
void Swap( BlasInt n, T* x, BlasInt incx, T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T temp;
    for( BlasInt i=0; i<n; ++i )
    {
        temp = x[i*incx];
        x[i*incx] = y[i*incy];
        y[i*incy] = temp;
    }
}
template void Swap( BlasInt n, Int* x, BlasInt incx, Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Swap
( BlasInt n, DoubleDouble* x, BlasInt incx, DoubleDouble* y, BlasInt incy );
template void Swap
( BlasInt n, QuadDouble* x, BlasInt incx, QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Swap( BlasInt n, Quad* x, BlasInt incx, Quad* y, BlasInt incy );
template void Swap
( BlasInt n, Complex<Quad>* x, BlasInt incx, Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Swap
( BlasInt n, BigInt* x, BlasInt incx, BigInt* y, BlasInt incy );
template void Swap
( BlasInt n, BigFloat* x, BlasInt incx, BigFloat* y, BlasInt incy );
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
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    // TODO: Special-case alpha=0, alpha=1, and alpha=-1?
    T gamma, delta;
    if( std::toupper(trans) == 'N' )
    {
        if( m > 0 && n == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<m; ++i )
                y[i*incy] = 0;
            return;
        }

        Scal( m, beta, y, incy );
        for( BlasInt j=0; j<n; ++j )
        {
            gamma = x[j*incx];
            gamma *= alpha;
            for( BlasInt i=0; i<m; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                delta = A[i+j*ALDim];
                delta *= gamma;
                y[i*incy] += delta;
            }
        }
    }
    else if( std::toupper(trans) == 'T' )
    {
        if( n > 0 && m == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<n; ++i )
                y[i*incy] = 0;
            return;
        }
        Scal( n, beta, y, incy );

        // Prescale x to avoid a factor of two more work than necessary
        vector<T> xAlpha(m);
        for( Int j=0; j<m; ++j )
        {
            xAlpha[j] = x[j*incx];
            xAlpha[j] *= alpha;
        }

        for( BlasInt i=0; i<n; ++i )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                // y[i*incy] += alpha*A[j+i*ALDim]*x[j*incx];
                gamma = A[j+i*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }
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

        // Prescale x to avoid a factor of two more work than necessary
        vector<T> xAlpha(m);
        for( Int j=0; j<m; ++j )
        {
            xAlpha[j] = x[j*incx];
            xAlpha[j] *= alpha;
        }

        for( BlasInt i=0; i<n; ++i )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                // y[i*incy] += alpha*Conj(A[j+i*ALDim])*x[j*incx];
                Conj( A[j+i*ALDim], gamma );
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }
    }
}
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Int& alpha,
  const Int* A, BlasInt ALDim,
  const Int* x, BlasInt incx, 
  const Int& beta,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble& beta,
        DoubleDouble* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble& beta,
        QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
  const Quad* x, BlasInt incx, 
  const Quad& beta,
        Quad* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>& beta,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
  const BigInt* x, BlasInt incx, 
  const BigInt& beta,
        BigInt* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
  const BigFloat* x, BlasInt incx, 
  const BigFloat& beta,
        BigFloat* y, BlasInt incy );
#endif

void Gemv
( char trans, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy )
{ EL_BLAS(cgemv)
  ( &trans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Gemv
( char trans, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zgemv)
  ( &trans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

template<typename T>
void Ger
( BlasInt m, BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    // TODO: Special-case alpha=0?
    T gamma, delta;
    for( BlasInt j=0; j<n; ++j )
    {
        Conj( y[j*incy], gamma );
        gamma *= alpha;
        for( BlasInt i=0; i<m; ++i )
        {
            // A[i+j*ALDim] += alpha*x[i*incx]*Conj(y[j*incy]);
            delta = x[i*incx];
            delta *= gamma;
            A[i+j*ALDim] += delta;
        }
    }
}
template void Ger
( BlasInt m, BlasInt n, 
  const Int& alpha,
  const Int* x, BlasInt incx, 
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Ger
( BlasInt m, BlasInt n, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Ger
( BlasInt m, BlasInt n, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Ger
( BlasInt m, BlasInt n, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Ger
( BlasInt m, BlasInt n, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
#endif

void Ger
( BlasInt m, BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Ger
( BlasInt m, BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim  )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Ger
( BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Ger
( BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

template<typename T>
void Geru
( BlasInt m, BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    // TODO: Special-case alpha=0?
    T gamma, delta;
    for( BlasInt j=0; j<n; ++j )
    {
        gamma = y[j*incy];
        gamma *= alpha;
        for( BlasInt i=0; i<m; ++i )
        {
            // A[i+j*ALDim] += alpha*x[i*incx]*y[j*incy];
            delta = x[i*incx];
            delta *= gamma;
            A[i+j*ALDim] += delta;
        }
    }
}
template void Geru
( BlasInt m, BlasInt n, 
  const Int& alpha,
  const Int* x, BlasInt incx, 
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Geru
( BlasInt m, BlasInt n, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Geru
( BlasInt m, BlasInt n, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Geru
( BlasInt m, BlasInt n, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Geru
( BlasInt m, BlasInt n, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
#endif

void Geru
( BlasInt m, BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Geru
( BlasInt m, BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim )
{ EL_BLAS(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Geru
( BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Geru
( BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &ALDim ); }

// TODO: Introduce some sort of blocking
template<typename T>
void Hemv
( char uplo, BlasInt m,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

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

    // Pre-scale x to avoid redundant computation
    vector<T> xAlpha(m);
    for( Int i=0; i<m; ++i )
    {
        xAlpha[i] = x[i*incx];
        xAlpha[i] *= alpha;
    }

    T gamma;
    if( uplo == LOWER )
    {
        // Multiply with the lower triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=j; i<m; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the adjoint of the strictly lower triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=j+1; i<m; ++i )
            {
                // y[j*incy] += alpha*Conj(A[i+j*ALDim])*x[i*incx];
                Conj( A[i+j*ALDim], gamma );
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
    else
    {
        // Multiply with the upper triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=0; i<=j; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the adjoint of the strictly upper triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=0; i<j; ++i )
            {
                // y[j*incy] += alpha*Conj(A[i+j*ALDim])*x[i*incx];
                Conj( A[i+j*ALDim], gamma );
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
}

template void Hemv
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* x, BlasInt incx, 
  const Int& beta,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Hemv
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble& beta,
        DoubleDouble* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble& beta,
        QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Hemv
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* x, BlasInt incx, 
  const Quad& beta,
        Quad* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>& beta,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Hemv
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* x, BlasInt incx, 
  const BigInt& beta,
        BigInt* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* x, BlasInt incx, 
  const BigFloat& beta,
        BigFloat* y, BlasInt incy );
#endif

void Hemv
( char uplo, BlasInt m,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy )
{ EL_BLAS(chemv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zhemv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

template<typename T>
void Her
( char uplo, BlasInt m,
  const Base<T>& alpha,
  const T* x, BlasInt incx,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta;
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( x[j*incx], gamma );
            gamma *= alpha;
            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*Conj(x[j*incx]);
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( x[j*incx], gamma );
            gamma *= alpha;
            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*Conj(x[j*incx]);
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
}
template void Her
( char uplo, BlasInt m,
  const Int& alpha,
  const Int* x, BlasInt incx, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Her
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx,
        DoubleDouble* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx,
        QuadDouble* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Her
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx,
        Quad* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Complex<Quad>* x, BlasInt incx, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Her
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx,
        BigInt* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx,
        BigFloat* A, BlasInt ALDim );
#endif

void Her
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
        float* A, BlasInt ALDim )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Her
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Her
( char uplo, BlasInt m,
  const float& alpha,
  const scomplex* x, BlasInt incx,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cher)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Her
( char uplo, BlasInt m,
  const double& alpha,
  const dcomplex* x, BlasInt incx,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zher)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

template<typename T>
void Her2
( char uplo, BlasInt m, 
  const T& alpha,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta, phi, psi, alphaConj;
    Conj( alpha, alphaConj );
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( y[j*incy], gamma );
            gamma *= alpha;
            Conj( x[j*incx], delta );
            delta *= alphaConj;
            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] += alpha*      x[i*incx]*Conj(y[j*incy]) + 
                //                 Conj(alpha)*y[i*incy]*Conj(x[j*incx]);
                phi = x[i*incx]; 
                phi *= gamma;
                psi = y[i*incy];
                psi *= delta;
                phi += psi;
                A[i+j*ALDim] += phi;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( y[j*incy], gamma );
            gamma *= alpha;
            Conj( x[j*incx], delta );
            delta *= alphaConj;
            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*      x[i*incx]*Conj(y[j*incy]) + 
                //                 Conj(alpha)*y[i*incy]*Conj(x[j*incx]);
                phi = x[i*incx]; 
                phi *= gamma;
                psi = y[i*incy];
                psi *= delta;
                phi += psi;
                A[i+j*ALDim] += phi;
            }
        }
    }
}
template void Her2
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* x, BlasInt incx, 
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Her2
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Her2
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Her2
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
#endif

void Her2
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Her2
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Her2
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Her2
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

// TODO: Introduce some sort of blocking
template<typename T>
void Symv
( char uplo, BlasInt m,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

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

    // Pre-scale x to avoid redundant computation
    vector<T> xAlpha(m);
    for( Int i=0; i<m; ++i )
    {
        xAlpha[i] = x[i*incx];
        xAlpha[i] *= alpha;
    }

    T gamma;
    if( uplo == LOWER )
    {
        // Multiply with the lower triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=j; i<m; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the transpose of the strictly lower triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=j+1; i<m; ++i )
            {
                // y[j*incy] += alpha*A[i+j*ALDim]*x[i*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
    else
    {
        // Multiply with the upper triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=0; i<=j; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the transpose of the strictly upper triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=0; i<j; ++i )
            {
                // y[j*incy] += alpha*A[i+j*ALDim]*x[i*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
}

template void Symv
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* x, BlasInt incx, 
  const Int& beta,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Symv
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble& beta,
        DoubleDouble* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble& beta,
        QuadDouble* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Symv
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* x, BlasInt incx, 
  const Quad& beta,
        Quad* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>& beta,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Symv
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* x, BlasInt incx, 
  const BigInt& beta,
        BigInt* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* x, BlasInt incx, 
  const BigFloat& beta,
        BigFloat* y, BlasInt incy );
#endif

void Symv
( char uplo, BlasInt m,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, BlasInt m,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    EL_LAPACK(csymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

void Symv
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    EL_LAPACK(zsymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

template<typename T>
void Syr
( char uplo, BlasInt m,
  const T& alpha,
  const T* x, BlasInt incx,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta;
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = x[j*incx];
            gamma *= alpha;
            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*x[j*incx];
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = x[j*incx];
            gamma *= alpha;
            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*x[j*incx];
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
}
template void Syr
( char uplo, BlasInt m,
  const Int& alpha,
  const Int* x, BlasInt incx, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Syr
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx,
        DoubleDouble* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx,
        QuadDouble* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syr
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx,
        Quad* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Syr
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx,
        BigInt* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx,
        BigFloat* A, BlasInt ALDim );
#endif

void Syr
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
        float* A, BlasInt ALDim  )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Syr
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Syr
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx,
        scomplex* A, BlasInt ALDim )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    EL_LAPACK(csyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); 
}

void Syr
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx,
        dcomplex* A, BlasInt ALDim )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    EL_LAPACK(zsyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); 
}

template<typename T>
void Syr2
( char uplo, BlasInt m, 
  const T& alpha,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta, phi, psi;
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = y[j*incy];
            gamma *= alpha;
            delta = x[j*incx];
            delta *= alpha;
            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] +=
                //   alpha*(x[i*incx]*y[j*incy]+y[i*incy]*x[j*incx])
                phi = x[i*incx]; 
                phi *= gamma;
                psi = y[i*incy];
                psi *= delta;
                phi += psi;
                A[i+j*ALDim] += phi;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = y[j*incy];
            gamma *= alpha;
            delta = x[j*incx]; 
            delta *= alpha;
            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*y[j*incy] + 
                //                 alpha*y[i*incy]*x[j*incx];
                phi = x[i*incx]; 
                phi *= gamma;
                psi = y[i*incy];
                psi *= delta;
                phi += psi;
                A[i+j*ALDim] += phi;
            }
        }
    }
}
template void Syr2
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* x, BlasInt incx,
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Syr2
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syr2
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Syr2
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
#endif

void Syr2
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Syr2
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Syr2
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{
    // csyr2 doesn't exist, so we route through csyr2k. However, csyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const BlasInt k = 1;
    const scomplex beta = 1.;
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &ALDim );
}

void Syr2
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{
    // zsyr2 doesn't exist, so we route through zsyr2k. However, zsyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const BlasInt k = 1;
    const dcomplex beta = 1.;
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &ALDim );
}

// TODO: Add correctness tests
template<typename T>
void Trmv
( char uplo, char trans, char diag,
  BlasInt m,
  const T* A, BlasInt ALDim,
        T* x, BlasInt incx )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    const bool lower = ( std::toupper(uplo) == 'L' );
    const bool conj = ( std::toupper(trans) == 'C' );
    const bool unitDiag = ( std::toupper(diag) == 'U' );
    const T zero(0);

    T gamma, delta;
    if( lower )
    {
        if( std::toupper(trans) == 'N' )
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                gamma = x[j*incx]; 
                if( gamma != zero )
                {
                    for( BlasInt i=m-1; i>j; --i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= gamma;
                        x[i*incx] += delta;
                    }
                    if( !unitDiag )
                        x[j*incx] *= A[j+j*ALDim];
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<m; ++j )
            {
                gamma = x[j*incx]; 
                if( conj )
                {
                    if( !unitDiag )
                    {
                        Conj( A[j+j*ALDim], delta );
                        gamma *= delta;
                    }
                    for( BlasInt i=j+1; i<m; ++i )
                    {
                        Conj( A[i+j*ALDim], delta );
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                else
                {
                    if( !unitDiag )
                        gamma *= A[j+j*ALDim];
                    for( BlasInt i=j+1; i<m; ++i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                x[j*incx] = gamma;
            }
        }
    }
    else
    {
        if( std::toupper(trans) == 'N' )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                gamma = x[j*incx]; 
                if( gamma != zero )
                {
                    for( BlasInt i=0; i<j; ++i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= gamma;
                        x[i*incx] += delta;
                    }
                    if( !unitDiag )
                        x[j*incx] *= A[j+j*ALDim];
                }
            }
        }
        else
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                gamma = x[j*incx]; 
                if( conj )
                {
                    if( !unitDiag )
                    {
                        Conj( A[j+j*ALDim], delta );
                        gamma *= delta;
                    }
                    for( BlasInt i=j-1; i>=0; --i )
                    {
                        Conj( A[i+j*ALDim], delta );
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                else
                {
                    if( !unitDiag )
                        gamma *= A[j+j*ALDim];
                    for( BlasInt i=j-1; i>=0; --i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                x[j*incx] = gamma;
            }
        }
    }
}
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Int* A, BlasInt ALDim,
        Int* x, BlasInt incx );
#ifdef EL_HAVE_QD
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const DoubleDouble* A, BlasInt ALDim,
        DoubleDouble* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const QuadDouble* A, BlasInt ALDim,
        QuadDouble* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Quad* A, BlasInt ALDim,
        Quad* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<Quad>* A, BlasInt ALDim,
        Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const BigInt* A, BlasInt ALDim,
        BigInt* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const BigFloat* A, BlasInt ALDim,
        BigFloat* x, BlasInt incx );
#endif

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt ALDim,
        float* x, BlasInt incx )
{ EL_BLAS(strmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt ALDim,
        double* x, BlasInt incx )
{ EL_BLAS(dtrmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt ALDim,
        scomplex* x, BlasInt incx )
{ EL_BLAS(ctrmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* x, BlasInt incx )
{ EL_BLAS(ztrmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

template<typename F>
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const F* A, BlasInt ALDim,
        F* x, BlasInt incx )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    const bool conj = ( std::toupper(trans) == 'C' );
    const bool lower = ( std::toupper(uplo) == 'L' );
    const bool unitDiag = ( std::toupper(diag) == 'U' );

    F gamma, delta;
    if( lower )
    {
        if( std::toupper(trans) == 'N' )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                if( !unitDiag ) 
                    x[j*incx] /= A[j+j*ALDim];
                gamma = x[j*incx];
                for( BlasInt i=j+1; i<m; ++i )
                {
                    delta = A[i+j*ALDim];
                    delta *= gamma;
                    x[i*incx] -= delta;
                }
            }
        }
        else
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                if( conj )
                {
                    if( !unitDiag )
                    {
                        Conj( A[j+j*ALDim], delta );
                        x[j*incx] /= delta;
                    }
                    gamma = x[j*incx];
                    for( BlasInt i=j-1; i>=0; --i )
                    {
                        Conj( A[j+i*ALDim], delta );
                        delta *= gamma;
                        x[i*incx] -= delta;
                    }
                }
                else
                {
                    if( !unitDiag )
                        x[j*incx] /= A[j+j*ALDim];
                    gamma = x[j*incx];
                    for( BlasInt i=j-1; i>=0; --i )
                    {
                        delta = A[j+i*ALDim];
                        delta *= gamma;
                        x[i*incx] -= delta;
                    }
                }
            }
        }
    }
    else
    {
        if( std::toupper(trans) == 'N' )
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                if( !unitDiag )
                    x[j*incx] /= A[j+j*ALDim];
                gamma = x[j*incx];
                for( BlasInt i=j-1; i>=0; --i )
                {
                    delta = A[i+j*ALDim];
                    delta *= gamma;
                    x[i*incx] -= delta;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<m; ++j )
            {
                if( conj )
                {
                    if( !unitDiag ) 
                        x[j*incx] /= A[j+j*ALDim];
                    gamma = x[j*incx];
                    for( BlasInt i=j+1; i<m; ++i )
                    {
                        delta = A[j+i*ALDim];
                        delta *= gamma;
                        x[i*incx] -= delta;
                    }
                }
                else
                {
                    if( !unitDiag ) 
                    {
                        Conj( A[j+j*ALDim], delta );
                        x[j*incx] /= delta;
                    }
                    gamma = x[j*incx];
                    for( BlasInt i=j+1; i<m; ++i )
                    {
                        Conj( A[j+i*ALDim], delta );
                        delta *= gamma;
                        x[i*incx] -= delta;
                    }
                }
            }
        }
    }
}
#ifdef EL_HAVE_QD
template void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const DoubleDouble* A, BlasInt ALDim,
        DoubleDouble* x, BlasInt incx );
template void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const QuadDouble* A, BlasInt ALDim,
        QuadDouble* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const Quad* A, BlasInt ALDim,
        Quad* x, BlasInt incx );
template void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<Quad>* A, BlasInt ALDim,
        Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const BigFloat* A, BlasInt ALDim,
        BigFloat* x, BlasInt incx );
#endif

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt ALDim,
        float* x, BlasInt incx )
{ EL_BLAS(strsv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt ALDim,
        double* x, BlasInt incx )
{ EL_BLAS(dtrsv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt ALDim,
        scomplex* x, BlasInt incx )
{ EL_BLAS(ctrsv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* x, BlasInt incx )
{ EL_BLAS(ztrsv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

// Level 3 BLAS
// ============
template<typename T>
void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k,
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( m > 0 && n > 0 && k == 0 && beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] = 0;
        return;
    }

    // Scale C
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] *= beta;
    }

    // Naive implementation
    T gamma, delta;
    if( std::toupper(transA) == 'N' && std::toupper(transB) == 'N' )
    {
        // C := alpha A B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt l=0; l<k; ++l )
            {
                gamma = alpha*B[l+j*BLDim];
                for( BlasInt i=0; i<m; ++i )
                {
                    delta = A[i+l*ALDim];
                    delta *= gamma;
                    C[i+j*CLDim] += delta;
                }
            }
        }
    }
    else if( std::toupper(transA) == 'N' )
    {
        if( std::toupper(transB) == 'T' )
        {
            // C := alpha A B^T + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt l=0; l<k; ++l )
                {
                    gamma = alpha*B[j+l*BLDim];
                    for( BlasInt i=0; i<m; ++i )
                    {
                        delta = A[i+l*ALDim];
                        delta *= gamma;
                        C[i+j*CLDim] += delta;
                    }
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
                    gamma = alpha*Conj(B[j+l*BLDim]);
                    for( BlasInt i=0; i<m; ++i )
                    {
                        delta = A[i+l*ALDim];
                        delta *= gamma;
                        C[i+j*CLDim] += delta;
                    }
                }
            }
        }
    }
    else if( std::toupper(transB) == 'N' )
    {
        if( std::toupper(transA) == 'T' )
        {
            // C := alpha A^T B + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[l+i*ALDim];
                        delta *= B[l+j*BLDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
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
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );
                        delta *= B[l+j*BLDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        if( std::toupper(transA) == 'T' && std::toupper(transB) == 'T' )
        {
            // C := alpha A^T B^T + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[l+i*ALDim];
                        delta *= B[j+l*BLDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else if( std::toupper(transA) == 'T' )
        {
            // C := alpha A^T B^H + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( B[j+l*BLDim], delta );
                        delta *= A[l+i*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else if( std::toupper(transB) == 'T' )
        {
            // C := alpha A^H B^T + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );
                        delta *= B[j+l*BLDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            // C := alpha A^H B^H + C
            T phi;
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );       
                        Conj( B[j+l*BLDim], phi );
                        delta *= phi;
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim,
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
#endif

void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");
      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    const char fixedTransA = ( std::toupper(transA) == 'C' ? 'T' : transA );
    const char fixedTransB = ( std::toupper(transB) == 'C' ? 'T' : transB );
    EL_BLAS(sgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");

      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }      

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    const char fixedTransA = ( std::toupper(transA) == 'C' ? 'T' : transA );
    const char fixedTransB = ( std::toupper(transB) == 'C' ? 'T' : transB );
    EL_BLAS(dgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");

      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }      

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    EL_BLAS(cgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    DEBUG_ONLY(
      CSE cse("blas::Gemm");

      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }      

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    EL_BLAS(zgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

template<typename T>
void Hemm
( char side, char uplo,
  BlasInt m, BlasInt n, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

    // Scale C
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] = 0;

    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] *= beta;
    }

    // Naive implementation
    T gamma, delta;
    if( std::toupper(side) == 'L' && std::toupper(uplo) == 'L' )
    {
        // C := alpha tril(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=i; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha tril(A,-1)^H B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i+1; l<m; ++l )
                {
                    Conj( A[l+i*ALDim], delta );
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'L' && std::toupper(uplo) == 'U' )
    {
        // C := alpha triu(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i; l<m; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha triu(A,-)^H B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<i; ++l )
                {
                    Conj( A[l+i*ALDim], delta );
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'L' )
    {
        // C := alpha B tril(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j; l<n; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B tril(A,-1)^H + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<j; ++l )
                {
                    Conj( A[j+l*ALDim], delta );
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'U' )
    {
        // C := alpha B triu(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=j; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B triu(A,1)^H + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j+1; l<n; ++l )
                {
                    Conj( A[j+l*ALDim], delta );
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    DEBUG_ONLY(
      else
          LogicError("Unsuported Hemm option");
    )
}
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
#endif

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(chemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zhemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

template<typename T>
void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const Base<T>& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( beta == Base<T>(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] = 0;
    }
    else if( beta != Base<T>(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] *= beta;
    }

    const T alphaConj = Conj(alpha);
    const bool normal = ( std::toupper(trans) == 'N' );
    const bool lower = ( std::toupper(uplo) == 'L' );

    T gamma, delta, phi;
    if( normal )
    {
        // C := alpha A B^H + Conj(alpha) B A^H + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( B[j+l*BLDim], phi );
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        Conj( A[j+l*ALDim], phi );
                        phi *= B[i+l*BLDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( B[j+l*BLDim], phi );
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        Conj( A[j+l*ALDim], phi );
                        phi *= B[i+l*BLDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^H B + Conj(alpha) B^H A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], phi );
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        Conj( B[l+i*BLDim], phi );
                        phi *= A[l+j*ALDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
                    C[i+j*CLDim] += gamma; 
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], phi );
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        Conj( B[l+i*BLDim], phi );
                        phi *= A[l+j*ALDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Quad& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
#endif

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(ssyr2k)
    ( &uplo, &transFixed, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dsyr2k)
    ( &uplo, &transFixed, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const float& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(cher2k)
    ( &uplo, &trans, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const double& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zher2k)
    ( &uplo, &trans, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

template<typename T>
void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Base<T>& alpha,
  const T* A, BlasInt ALDim, 
  const Base<T>& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( beta == Base<T>(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] = 0;
    }
    else if( beta != Base<T>(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] *= beta;
    }

    const bool normal = ( std::toupper(trans) == 'N' );
    const bool lower = ( std::toupper(uplo) == 'L' );

    T gamma, delta;
    if( normal )
    {
        // C := alpha A A^H + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[j+l*ALDim], delta );
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[j+l*ALDim], delta );
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^H A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Quad& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Quad& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
#endif

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    const char transFixed = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(ssyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    const char transFixed = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(dsyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const scomplex* A, BlasInt ALDim,
  const float& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(cherk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const dcomplex* A, BlasInt ALDim,
  const double& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zherk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

template<typename T>
void Symm
( char side, char uplo,
  BlasInt m, BlasInt n, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

    // Scale C
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] = 0;

    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] *= beta;
    }

    // Naive implementation
    T gamma, delta;
    if( std::toupper(side) == 'L' && std::toupper(uplo) == 'L' )
    {
        // C := alpha tril(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=i; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha tril(A,-1)^T B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i+1; l<m; ++l )
                {
                    delta = A[l+i*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'L' && std::toupper(uplo) == 'U' )
    {
        // C := alpha triu(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i; l<m; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha triu(A,-)^T B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<i; ++l )
                {
                    delta = A[l+i*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'L' )
    {
        // C := alpha B tril(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j; l<n; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B tril(A,-1)^T + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<j; ++l )
                {
                    delta = A[j+l*ALDim];
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'U' )
    {
        // C := alpha B triu(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=j; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B triu(A,1)^T + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j+1; l<n; ++l )
                {
                    delta = A[j+l*ALDim];
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    DEBUG_ONLY(
      else
          LogicError("Unsuported Symm option");
    )
}
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
#endif

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(csymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

template<typename T>
void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] *= beta;
    }

    const bool normal = ( std::toupper(trans) == 'N' );
    const bool lower = ( std::toupper(uplo) == 'L' );

    T gamma, phi;
    if( normal )
    {
        // C := alpha A B^T + alpha B A^T + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        phi = B[j+l*BLDim];
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        phi = A[j+l*ALDim];
                        phi *= B[i+l*BLDim];
                        gamma += phi;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        phi = B[j+l*BLDim];
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        phi = A[j+l*ALDim];
                        phi *= B[i+l*BLDim];
                        gamma += phi;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^T B + alpha B^T A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        phi = A[l+i*ALDim];
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        phi = B[l+i*BLDim];
                        phi *= A[l+j*ALDim];
                        gamma += phi;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma; 
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        phi = A[l+i*ALDim];
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        phi = B[l+i*BLDim];
                        phi *= A[l+j*ALDim];
                        gamma += phi;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
#endif

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

template<typename T>
void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] *= beta;
    }

    const bool normal = ( std::toupper(trans) == 'N' );
    const bool lower = ( std::toupper(uplo) == 'L' );

    T gamma, delta;
    if( normal )
    {
        // C := alpha A A^T + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[j+l*ALDim];
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[j+l*ALDim];
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^T A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[l+i*ALDim];
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[l+i*ALDim];
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
#endif

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(csyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zsyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

template<typename T>
void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim,
        T* B, BlasInt BLDim )
{
    const bool onLeft = ( std::toupper(side) == 'L' );
    const bool conjugate = ( std::toupper(trans) == 'C' );

    // Scale B
    for( BlasInt j=0; j<n; ++j )
        for( BlasInt i=0; i<m; ++i )
            B[i+j*BLDim] *= alpha;

    // TODO: Legitimate blocked implementations...it seems offensive to
    //       repeatedly stream all of the triangular matrix through memory
    //       for each row/column of B

    if( onLeft )
    {
        for( BlasInt j=0; j<n; ++j )
            Trmv( uplo, trans, unit, m, A, ALDim, &B[j*BLDim], 1 );
    }
    else
    {
        char newTrans = ( std::toupper(trans) == 'N' ? 'T' : 'N' );
        for( BlasInt i=0; i<m; ++i )
        {
            if( conjugate )
                for( Int j=0; j<n; ++j )
                    B[i+j*BLDim] = Conj(B[i+j*BLDim]);
            Trmv( uplo, newTrans, unit, n, A, ALDim, &B[i], BLDim );
            // TODO: Avoid temporaries here
            if( conjugate )
                for( Int j=0; j<n; ++j )
                    B[i+j*BLDim] = Conj(B[i+j*BLDim]);
        }
    }
}
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Int& alpha,
  const Int* A, BlasInt ALDim,
        Int* B, BlasInt BLDim );
#ifdef EL_HAVE_QD
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
        DoubleDouble* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
        QuadDouble* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
        Quad* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim,
        Complex<Quad>* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_MPC
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
        BigInt* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
        BigFloat* B, BlasInt BLDim );
#endif

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
        float* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );    
    EL_BLAS(strmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
        double* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );    
    EL_BLAS(dtrmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
        scomplex* B, BlasInt BLDim )
{
    EL_BLAS(ctrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* B, BlasInt BLDim )
{
    EL_BLAS(ztrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

template<typename F>
void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const F& alpha,
  const F* A, BlasInt ALDim,
        F* B, BlasInt BLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    const bool onLeft = ( std::toupper(side) == 'L' );
    const bool lower = ( std::toupper(uplo) == 'L' );
    const bool conjugate = ( std::toupper(trans) == 'C' );
    const bool unitDiag = ( std::toupper(unit) == 'U' );

    // Scale B
    for( BlasInt j=0; j<n; ++j )
        for( BlasInt i=0; i<m; ++i )
            B[i+j*BLDim] *= alpha;

    F alpha11, alpha11Conj;
    if( onLeft )
    {
        if( std::toupper(trans) == 'N' )
        {
            if( lower )
            {
                for( BlasInt k=0; k<m; ++k )     
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )      
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( m-(k+1), n,
                      F(-1), &A[(k+1)+k*ALDim], 1,
                             &B[k],             BLDim,
                             &B[k+1],           BLDim );
                }
            }
            else
            {
                for( BlasInt k=m-1; k>=0; --k )     
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )      
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( k, n,
                      F(-1), &A[k*ALDim], 1,
                             &B[k],       BLDim,
                              B,          BLDim );
                }
            }
        }
        else if( conjugate )
        {
            if( lower )
            {
               std::vector<F> aRow(m);
               for( BlasInt k=m-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11Conj;
                    }
                    for( Int s=0; s<k; ++s )
                        Conj( A[k+s*ALDim], aRow[s] );
                    Geru
                    ( k, n,
                      F(-1), aRow.data(), 1,
                             &B[k],       BLDim,
                              B,          BLDim );
                }
            }
            else
            {
                for( BlasInt k=0; k<m; ++k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11Conj;
                    }
                    std::vector<F> aRow(m);
                    for( Int s=0; s<m-(k+1); ++s )
                        Conj( A[k+(k+1+s)*ALDim], aRow[s] );
                    Geru
                    ( m-(k+1), n,
                      F(-1), aRow.data(), 1,
                             &B[k],       BLDim,
                             &B[k+1],     BLDim );
                }
            }
        }
        else
        {
            if( lower )
            {
                for( BlasInt k=m-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( k, n,
                      F(-1), &A[k], ALDim,
                             &B[k], BLDim,
                              B,    BLDim );
                }
            }
            else
            {
                for( BlasInt k=0; k<m; ++k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( m-(k+1), n,
                      F(-1), &A[k+(k+1)*ALDim], ALDim,
                             &B[k],             BLDim,
                             &B[k+1],           BLDim );
                }
            }
        }
    }
    else
    {
        if( std::toupper(trans) == 'N' )
        {
            if( lower )
            {
                for( BlasInt k=n-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, k,
                      F(-1), &B[k*BLDim], 1,
                             &A[k],       ALDim,
                              B,          BLDim );
                }
            }
            else
            {
                for( BlasInt k=0; k<n; ++k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, n-(k+1),
                      F(-1), &B[k*BLDim],       1,
                             &A[k+(k+1)*ALDim], ALDim,
                             &B[(k+1)*BLDim],   BLDim );
                }
            }
        }
        else if( conjugate )
        {
            if( lower )
            {
                for( BlasInt k=0; k<n; ++k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11Conj;
                    }
                    blas::Ger
                    ( m, n-(k+1),
                      F(-1), &B[k*BLDim],       1,
                             &A[(k+1)+k*ALDim], 1,
                             &B[(k+1)*BLDim],   BLDim );
                }
            }
            else
            {
                for( BlasInt k=n-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11Conj;
                    }
                    blas::Ger
                    ( m, k,
                      F(-1), &B[k*BLDim], 1,
                             &A[k*ALDim], 1,
                              B,          BLDim );
                }
            }
        }
        else
        {
            if( lower )
            {
                for( BlasInt k=0; k<n; ++k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, n-(k+1),
                      F(-1), &B[k*BLDim],       1,
                             &A[(k+1)+k*ALDim], 1,
                             &B[(k+1)*BLDim],   BLDim );
                }
            }
            else
            {
                for( BlasInt k=n-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, k,
                      F(-1), &B[k*BLDim], 1,
                             &A[k*ALDim], 1,
                              B,          BLDim );
                }
            }
        }
    }
}
#ifdef EL_HAVE_QD
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
        DoubleDouble* B, BlasInt BLDim );
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
        QuadDouble* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
        Quad* B, BlasInt BLDim );
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim,
        Complex<Quad>* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_MPC
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
        BigFloat* B, BlasInt BLDim );
#endif

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
        float* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(strsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
        double* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(dtrsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
        scomplex* B, BlasInt BLDim )
{
    EL_BLAS(ctrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* B, BlasInt BLDim )
{
    EL_BLAS(ztrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

} // namespace blas
} // namespace El
