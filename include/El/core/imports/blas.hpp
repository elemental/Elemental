/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_BLAS_HPP
#define EL_IMPORTS_BLAS_HPP

namespace El {

#ifdef EL_USE_64BIT_BLAS_INTS
typedef long long int BlasInt;
#else
typedef int BlasInt;
#endif

namespace blas {

// NOTE: templated routines are custom and not wrappers

// Level 1 BLAS 
// ============
template<typename T>
void Axpy( BlasInt n, T alpha, const T* x, BlasInt incx, T* y, BlasInt incy );

void Axpy
( BlasInt n, float    alpha, 
  const float   * x, BlasInt incx, float   * y, BlasInt incy );
void Axpy
( BlasInt n, double   alpha, 
  const double  * x, BlasInt incx, double  * y, BlasInt incy );
void Axpy
( BlasInt n, scomplex alpha, 
  const scomplex* x, BlasInt incx, scomplex* y, BlasInt incy );
void Axpy
( BlasInt n, dcomplex alpha, 
  const dcomplex* x, BlasInt incx, dcomplex* y, BlasInt incy );
template<typename T>
void Axpy( BlasInt n, T alpha, const T* x, BlasInt incx, T* y, BlasInt incy );

template<typename T>
void Copy( BlasInt n, const T* x, BlasInt incx, T* y, BlasInt incy );

void Copy
( BlasInt n, const float   * x, BlasInt incx, float   * y, BlasInt incy );
void Copy
( BlasInt n, const double  * x, BlasInt incx, double  * y, BlasInt incy );
void Copy
( BlasInt n, const scomplex* x, BlasInt incx, scomplex* y, BlasInt incy );
void Copy
( BlasInt n, const dcomplex* x, BlasInt incx, dcomplex* y, BlasInt incy );
template<typename T>
void Copy( BlasInt n, const T* x, BlasInt incx, T* y, BlasInt incy );

float    Dot
( BlasInt n, const float   * x, BlasInt incx, const float   * y, BlasInt incy );
double   Dot
( BlasInt n, const double  * x, BlasInt incx, const double  * y, BlasInt incy );
scomplex Dot
( BlasInt n, const scomplex* x, BlasInt incx, const scomplex* y, BlasInt incy );
dcomplex Dot
( BlasInt n, const dcomplex* x, BlasInt incx, const dcomplex* y, BlasInt incy );
template<typename T>
T Dot( BlasInt n, const T* x, BlasInt incx, const T* y, BlasInt incy );

float    Dotc
( BlasInt n, const float   * x, BlasInt incx, const float   * y, BlasInt incy );
double   Dotc
( BlasInt n, const double  * x, BlasInt incx, const double  * y, BlasInt incy );
scomplex Dotc
( BlasInt n, const scomplex* x, BlasInt incx, const scomplex* y, BlasInt incy );
dcomplex Dotc
( BlasInt n, const dcomplex* x, BlasInt incx, const dcomplex* y, BlasInt incy );
template<typename T>
T Dotc( BlasInt n, const T* x, BlasInt incx, const T* y, BlasInt incy );

float    Dotu
( BlasInt n, const float   * x, BlasInt incx, const float   * y, BlasInt incy );
double   Dotu
( BlasInt n, const double  * x, BlasInt incx, const double  * y, BlasInt incy );
scomplex Dotu
( BlasInt n, const scomplex* x, BlasInt incx, const scomplex* y, BlasInt incy );
dcomplex Dotu
( BlasInt n, const dcomplex* x, BlasInt incx, const dcomplex* y, BlasInt incy );
template<typename T>
T Dotu( BlasInt n, const T* x, BlasInt incx, const T* y, BlasInt incy );

template<typename F>
Base<F> Nrm2( BlasInt n, const F* x, BlasInt incx );

float  Nrm2( BlasInt n, const float   * x, BlasInt incx );
double Nrm2( BlasInt n, const double  * x, BlasInt incx );
float  Nrm2( BlasInt n, const scomplex* x, BlasInt incx );
double Nrm2( BlasInt n, const dcomplex* x, BlasInt incx );

float Givens
( float alpha, float beta, float* c, float* s );
double Givens
( double alpha, double beta, double* c, double* s );
scomplex Givens
( scomplex alpha, scomplex beta, float* c, scomplex* s );
dcomplex Givens
( dcomplex alpha, dcomplex beta, double* c, dcomplex* s );

void Rot
( BlasInt n, float   * x, BlasInt incx, 
             float   * y, BlasInt incy, float  c, float    s );
void Rot
( BlasInt n, double  * x, BlasInt incx, 
             double  * y, BlasInt incy, double c, double   s );
void Rot
( BlasInt n, scomplex* x, BlasInt incx, 
             scomplex* y, BlasInt incy, float  c, scomplex s );
void Rot
( BlasInt n, dcomplex* x, BlasInt incx, 
             dcomplex* y, BlasInt incy, double c, dcomplex s );

template<typename T> 
void Scal( BlasInt n, T alpha,         T*  x, BlasInt incx );
template<typename T> 
void Scal( BlasInt n, T alpha, Complex<T>* x, BlasInt incx );
void Scal( BlasInt n, float    alpha, float   * x, BlasInt incx );
void Scal( BlasInt n, double   alpha, double  * x, BlasInt incx );
void Scal( BlasInt n, scomplex alpha, scomplex* x, BlasInt incx );
void Scal( BlasInt n, dcomplex alpha, dcomplex* x, BlasInt incx );

// NOTE: Nrm1 is not the official name but is consistent with Nrm2
float  Nrm1( BlasInt n, const float   * x, BlasInt incx );
double Nrm1( BlasInt n, const double  * x, BlasInt incx );
float  Nrm1( BlasInt n, const scomplex* x, BlasInt incx );
double Nrm1( BlasInt n, const dcomplex* x, BlasInt incx );

void Swap( BlasInt n, float   * x, BlasInt incx, float   * y, BlasInt incy );
void Swap( BlasInt n, double  * x, BlasInt incx, double  * y, BlasInt incy );
void Swap( BlasInt n, scomplex* x, BlasInt incx, scomplex* y, BlasInt incy );
void Swap( BlasInt n, dcomplex* x, BlasInt incx, dcomplex* y, BlasInt incy );
template<typename T> 
void Swap( BlasInt n, T* x, BlasInt incx, T* y, BlasInt incy );
            
// Level 2 BLAS
// ============
template<typename T>
void Gemv
( char trans, BlasInt m, BlasInt n,
  T alpha, const T* A, BlasInt lda, 
           const T* x, BlasInt incx,
  T beta,        T* y, BlasInt incy );

void Gemv
( char trans, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, 
               const float* x, BlasInt incx,
  float beta,        float* y, BlasInt incy );
void Gemv
( char trans, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, 
                const double* x, BlasInt incx,
  double beta,        double* y, BlasInt incy );
void Gemv
( char trans, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* x, BlasInt incx,
  scomplex beta,        scomplex* y, BlasInt incy );
void Gemv
( char trans, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* x, BlasInt incx,
  dcomplex beta,        dcomplex* y, BlasInt incy );

template<typename T>
void Ger
( BlasInt m, BlasInt n,
  T alpha, const T* x, BlasInt incx, 
           const T* y, BlasInt incy,
                 T* A, BlasInt lda );

void Ger
( BlasInt m, BlasInt n,
  float alpha, const float* x, BlasInt incx, 
               const float* y, BlasInt incy,
                     float* A, BlasInt lda );
void Ger
( BlasInt m, BlasInt n,
  double alpha, const double* x, BlasInt incx, 
                const double* y, BlasInt incy,
                      double* A, BlasInt lda );
void Ger
( BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* x, BlasInt incx, 
                  const scomplex* y, BlasInt incy,
                        scomplex* A, BlasInt lda );
void Ger
( BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* x, BlasInt incx, 
                  const dcomplex* y, BlasInt incy,
                        dcomplex* A, BlasInt lda );

template<typename T>
void Geru
( BlasInt m, BlasInt n,
  T alpha, const T* x, BlasInt incx, 
           const T* y, BlasInt incy,
                 T* A, BlasInt lda );

void Geru
( BlasInt m, BlasInt n,
  float alpha, const float* x, BlasInt incx, 
               const float* y, BlasInt incy,
                     float* A, BlasInt lda );
void Geru
( BlasInt m, BlasInt n,
  double alpha, const double* x, BlasInt incx, 
                const double* y, BlasInt incy,
                      double* A, BlasInt lda );
void Geru
( BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* x, BlasInt incx, 
                  const scomplex* y, BlasInt incy,
                        scomplex* A, BlasInt lda );
void Geru
( BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* x, BlasInt incx, 
                  const dcomplex* y, BlasInt incy,
                        dcomplex* A, BlasInt lda );

template<typename T>
void Hemv
( char uplo, BlasInt m,
  T alpha, const T* A, BlasInt lda, 
           const T* x, BlasInt incx,
  T beta,        T* y, BlasInt incy );

void Hemv
( char uplo, BlasInt m,
  float alpha, const float* A, BlasInt lda, 
               const float* x, BlasInt incx,
  float beta,        float* y, BlasInt incy );
void Hemv
( char uplo, BlasInt m,
  double alpha, const double* A, BlasInt lda, 
                const double* x, BlasInt incx,
  double beta,        double* y, BlasInt incy );
void Hemv
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* x, BlasInt incx,
  scomplex beta,        scomplex* y, BlasInt incy );
void Hemv
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* x, BlasInt incx,
  dcomplex beta,        dcomplex* y, BlasInt incy );

template<typename T>
void Her
( char uplo, BlasInt m,
  Base<T> alpha, const T* x, BlasInt incx, T* A, BlasInt lda );

void Her
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, float* A, BlasInt lda );
void Her
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, double* A, BlasInt lda );
void Her
( char uplo, BlasInt m,
  float alpha, const scomplex* x, BlasInt incx, scomplex* A, BlasInt lda );
void Her
( char uplo, BlasInt m,
  double alpha, const dcomplex* x, BlasInt incx, dcomplex* A, BlasInt lda );

template<typename T>
void Her2
( char uplo, BlasInt m,
  T alpha, const T* x, BlasInt incx, 
           const T* y, BlasInt incy,
                 T* A, BlasInt lda );

void Her2
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, 
               const float* y, BlasInt incy,
                     float* A, BlasInt lda );
void Her2
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, 
                const double* y, BlasInt incy,
                      double* A, BlasInt lda );
void Her2
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* x, BlasInt incx, 
                  const scomplex* y, BlasInt incy,
                        scomplex* A, BlasInt lda );
void Her2
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* x, BlasInt incx, 
                  const dcomplex* y, BlasInt incy,
                        dcomplex* A, BlasInt lda );

template<typename T>
void Symv
( char uplo, BlasInt m,
  T alpha, const T* A, BlasInt lda, 
           const T* x, BlasInt incx,
  T beta,        T* y, BlasInt incy );

void Symv
( char uplo, BlasInt m,
  float alpha, const float* A, BlasInt lda, 
               const float* x, BlasInt incx,
  float beta,        float* y, BlasInt incy );
void Symv
( char uplo, BlasInt m, 
  double alpha, const double* A, BlasInt lda, 
                const double* x, BlasInt incx,
  double beta,        double* y, BlasInt incy );
void Symv
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* x, BlasInt incx,
  scomplex beta,        scomplex* y, BlasInt incy );
void Symv
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* x, BlasInt incx,
  dcomplex beta,        dcomplex* y, BlasInt incy );

template<typename T>
void Syr
( char uplo, BlasInt m, T alpha, const T* x, BlasInt incx, T* A, BlasInt lda );

void Syr
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, float* A, BlasInt lda );
void Syr
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, double* A, BlasInt lda );
void Syr
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* x, BlasInt incx, scomplex* A, BlasInt lda ); 
void Syr
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* x, BlasInt incx, dcomplex* A, BlasInt lda );

template<typename T>
void Syr2
( char uplo, BlasInt m,
  T alpha, const T* x, BlasInt incx, 
           const T* y, BlasInt incy,
                 T* A, BlasInt lda );

void Syr2
( char uplo, BlasInt m,
  float alpha, const float* x, BlasInt incx, 
               const float* y, BlasInt incy,
                     float* A, BlasInt lda );
void Syr2
( char uplo, BlasInt m,
  double alpha, const double* x, BlasInt incx, 
                const double* y, BlasInt incy,
                      double* A, BlasInt lda );
void Syr2
( char uplo, BlasInt m,
  scomplex alpha, const scomplex* x, BlasInt incx, 
                  const scomplex* y, BlasInt incy,
                        scomplex* A, BlasInt lda );
void Syr2
( char uplo, BlasInt m,
  dcomplex alpha, const dcomplex* x, BlasInt incx, 
                  const dcomplex* y, BlasInt incy,
                        dcomplex* A, BlasInt lda );

template<typename T>
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const T* A, BlasInt lda, T* x, BlasInt incx );

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt lda, float* x, BlasInt incx );
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt lda, double* x, BlasInt incx );
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt lda, scomplex* x, BlasInt incx );
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt lda, dcomplex* x, BlasInt incx );

template<typename F>
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const F* A, BlasInt lda, F* x, BlasInt incx );

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt lda, float* x, BlasInt incx );
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt lda, double* x, BlasInt incx );
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt lda, scomplex* x, BlasInt incx );
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt lda, dcomplex* x, BlasInt incx );

// Level 3 BLAS
// ============
template<typename T>
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  T alpha, const T* A, BlasInt lda, 
           const T* B, BlasInt ldb,
  T beta,        T* C, BlasInt ldc );

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc );
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc );
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc );
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc );

template<typename T>
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  T alpha, const T* A, BlasInt lda, 
           const T* B, BlasInt ldb,
  T beta,        T* C, BlasInt ldc );

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc );
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc );
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc );
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc );

// TODO: Templated Her2k

void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc );
void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc );
void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  float beta,           scomplex* C, BlasInt ldc );
void Her2k
( char uplo, char trans, BlasInt n, BlasInt k,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  double beta,          dcomplex* C, BlasInt ldc );

template<typename T>
void Herk
( char uplo, BlasInt n, BlasInt k, 
  Base<T> alpha, const T* A, BlasInt lda, 
  Base<T> beta,        T* C, BlasInt ldc );

void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda, 
  float beta,        float* C, BlasInt ldc );
void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda, 
  double beta,        double* C, BlasInt ldc );
void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const scomplex* A, BlasInt lda,
  float beta,        scomplex* C, BlasInt ldc );
void Herk
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const dcomplex* A, BlasInt lda,
  double beta,        dcomplex* C, BlasInt ldc );

template<typename T>
void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  T alpha, const T* A, BlasInt lda, 
           const T* B, BlasInt ldb,
  T beta,        T* C, BlasInt ldc );

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc );
void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc );
void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc );
void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc );

// TODO: Templated Syr2k

void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda, 
               const float* B, BlasInt ldb,
  float beta,        float* C, BlasInt ldc );
void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda, 
                const double* B, BlasInt ldb,
  double beta,        double* C, BlasInt ldc );
void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  scomplex alpha, const scomplex* A, BlasInt lda, 
                  const scomplex* B, BlasInt ldb,
  scomplex beta,        scomplex* C, BlasInt ldc );
void Syr2k
( char uplo, char trans, BlasInt n, BlasInt k,
  dcomplex alpha, const dcomplex* A, BlasInt lda, 
                  const dcomplex* B, BlasInt ldb,
  dcomplex beta,        dcomplex* C, BlasInt ldc );

template<typename T>
void Syrk
( char uplo, BlasInt n, BlasInt k, 
  T alpha, const T* A, BlasInt lda, 
  T beta,        T* C, BlasInt ldc );

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  float alpha, const float* A, BlasInt lda,
  float beta,        float* C, BlasInt ldc );
void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  double alpha, const double* A, BlasInt lda,
  double beta,        double* C, BlasInt ldc );
void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  scomplex alpha, const scomplex* A, BlasInt lda,
  scomplex beta,        scomplex* C, BlasInt ldc );
void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  dcomplex alpha, const dcomplex* A, BlasInt lda,
  dcomplex beta,        dcomplex* C, BlasInt ldc );

// TODO: Templated Trmm

void Trmm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, float* B, BlasInt ldb );
void Trmm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, double* B, BlasInt ldb );
void Trmm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb );
void Trmm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb );

// TODO: Templated Trsm

void Trsm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, float* B, BlasInt ldb );
void Trsm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, double* B, BlasInt ldb );
void Trsm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb );
void Trsm
( char side,  char uplo, char trans, char unit, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb );

} // namespace blas
} // namespace El

#endif // ifndef EL_IMPORTS_BLAS_DECL_HPP
