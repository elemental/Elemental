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
namespace blas {

// NOTE: templated routines are custom and not wrappers

// Level 1 BLAS 
// ============
template<typename T>
void Axpy( int n, T alpha, const T* x, int incx, T* y, int incy );

void Axpy
( int n, float    alpha, const float   * x, int incx, float   * y, int incy );
void Axpy
( int n, double   alpha, const double  * x, int incx, double  * y, int incy );
void Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy );
void Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy );
template<typename T>
void Axpy( int n, T alpha, const T* x, int incx, T* y, int incy );

template<typename T>
void Copy( int n, const T* x, int incx, T* y, int incy );

void Copy( int n, const float   * x, int incx, float   * y, int incy );
void Copy( int n, const double  * x, int incx, double  * y, int incy );
void Copy( int n, const scomplex* x, int incx, scomplex* y, int incy );
void Copy( int n, const dcomplex* x, int incx, dcomplex* y, int incy );
template<typename T>
void Copy( int n, const T* x, int incx, T* y, int incy );

float    Dot( int n, const float   * x, int incx, const float   * y, int incy );
double   Dot( int n, const double  * x, int incx, const double  * y, int incy );
scomplex Dot( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dot( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dot( int n, const T* x, int incx, const T* y, int incy );

float    Dotc
( int n, const float   * x, int incx, const float   * y, int incy );
double   Dotc
( int n, const double  * x, int incx, const double  * y, int incy );
scomplex Dotc
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotc
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dotc( int n, const T* x, int incx, const T* y, int incy );

float    Dotu
( int n, const float   * x, int incx, const float   * y, int incy );
double   Dotu
( int n, const double  * x, int incx, const double  * y, int incy );
scomplex Dotu
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotu
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dotu( int n, const T* x, int incx, const T* y, int incy );

template<typename F>
Base<F> Nrm2( int n, const F* x, int incx );

float  Nrm2( int n, const float   * x, int incx );
double Nrm2( int n, const double  * x, int incx );
float  Nrm2( int n, const scomplex* x, int incx );
double Nrm2( int n, const dcomplex* x, int incx );

float Givens
( float alpha, float beta, float* c, float* s );
double Givens
( double alpha, double beta, double* c, double* s );
scomplex Givens
( scomplex alpha, scomplex beta, float* c, scomplex* s );
dcomplex Givens
( dcomplex alpha, dcomplex beta, double* c, dcomplex* s );

void Rot
( int n, float   * x, int incx, float   * y, int incy, float  c, float    s );
void Rot
( int n, double  * x, int incx, double  * y, int incy, double c, double   s );
void Rot
( int n, scomplex* x, int incx, scomplex* y, int incy, float  c, scomplex s );
void Rot
( int n, dcomplex* x, int incx, dcomplex* y, int incy, double c, dcomplex s );

template<typename T> void Scal( int n, T alpha,         T*  x, int incx );
template<typename T> void Scal( int n, T alpha, Complex<T>* x, int incx );
void Scal( int n, float    alpha, float   * x, int incx );
void Scal( int n, double   alpha, double  * x, int incx );
void Scal( int n, scomplex alpha, scomplex* x, int incx );
void Scal( int n, dcomplex alpha, dcomplex* x, int incx );

// NOTE: Nrm1 is not the official name but is consistent with Nrm2
float  Nrm1( int n, const float   * x, int incx );
double Nrm1( int n, const double  * x, int incx );
float  Nrm1( int n, const scomplex* x, int incx );
double Nrm1( int n, const dcomplex* x, int incx );

void Swap( int n, float   * x, int incx, float   * y, int incy );
void Swap( int n, double  * x, int incx, double  * y, int incy );
void Swap( int n, scomplex* x, int incx, scomplex* y, int incy );
void Swap( int n, dcomplex* x, int incx, dcomplex* y, int incy );
template<typename T> void Swap( int n, T* x, int incx, T* y, int incy );
            
// Level 2 BLAS
// ============
template<typename T>
void Gemv
( char trans, int m, int n,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy );

void Gemv
( char trans, int m, int n,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );
void Gemv
( char trans, int m, int n,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );
void Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );
void Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );

template<typename T>
void Ger
( int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
                                       T* A, int lda );

void Ger
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Ger
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Ger
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Ger
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );

template<typename T>
void Geru
( int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
                                       T* A, int lda );

void Geru
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Geru
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Geru
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Geru
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );

template<typename T>
void Hemv
( char uplo, int m,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,                             T* y, int incy );

void Hemv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );
void Hemv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );
void Hemv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );
void Hemv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );

template<typename T>
void Her
( char uplo, int m,
  Base<T> alpha, const T* x, int incx, T* A, int lda );

void Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda );
void Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda );
void Her
( char uplo, int m,
  float alpha, const scomplex* x, int incx, scomplex* A, int lda );
void Her
( char uplo, int m,
  double alpha, const dcomplex* x, int incx, dcomplex* A, int lda );

template<typename T>
void Her2
( char uplo, int m,
  T alpha, const T* x, int incx, 
           const T* y, int incy,
                 T* A, int lda );

void Her2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Her2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Her2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Her2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );

template<typename T>
void Symv
( char uplo, int m,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,                             T* y, int incy );

void Symv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );
void Symv
( char uplo, int m, 
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );
void Symv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );
void Symv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );

template<typename T>
void Syr
( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda );

void Syr
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda );
void Syr
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda );
void Syr
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda ); 
void Syr
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda );

template<typename T>
void Syr2
( char uplo, int m,
  T alpha, const T* x, int incx, 
           const T* y, int incy,
                 T* A, int lda );

void Syr2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Syr2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Syr2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Syr2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );

template<typename T>
void Trmv
( char uplo, char trans, char diag, int m,
  const T* A, int lda, T* x, int incx );

void Trmv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx );
void Trmv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx );
void Trmv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx );
void Trmv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx );

template<typename F>
void Trsv
( char uplo, char trans, char diag, int m,
  const F* A, int lda, F* x, int incx );

void Trsv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx );
void Trsv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx );
void Trsv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx );
void Trsv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx );

// Level 3 BLAS
// ============
template<typename T>
void Gemm
( char transA, char transB, int m, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Gemm
( char transA, char transB, int m, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Gemm
( char transA, char transB, int m, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Gemm
( char transA, char transB, int m, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Gemm
( char transA, char transB, int m, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );

template<typename T>
void Hemm
( char side, char uplo, int m, int n,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Hemm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Hemm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Hemm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Hemm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );

// TODO: Templated Her2k

void Her2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Her2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Her2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  float beta,           scomplex* C, int ldc );
void Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  double beta,          dcomplex* C, int ldc );

template<typename T>
void Herk
( char uplo, int n, int k, 
  Base<T> alpha, const T* A, int lda, 
  Base<T> beta,        T* C, int ldc );

void Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, float beta, float* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, double beta, double* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  float alpha, const scomplex* A, int lda,
  float beta,        scomplex* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  double alpha, const dcomplex* A, int lda,
  double beta,        dcomplex* C, int ldc );

template<typename T>
void Symm
( char side, char uplo, int m, int n,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Symm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Symm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Symm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Symm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );

// TODO: Templated Syr2k

void Syr2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Syr2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Syr2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Syr2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );

template<typename T>
void Syrk
( char uplo, int n, int k, 
  T alpha, const T* A, int lda, 
  T beta,        T* C, int ldc );

void Syrk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc );
void Syrk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc );
void Syrk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc );
void Syrk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc );

// TODO: Templated Trmm

void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb );
void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb );
void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb );
void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb );

// TODO: Templated Trsm

void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb );
void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb );
void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb );
void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb );

} // namespace blas
} // namespace El

#endif // ifndef EL_IMPORTS_BLAS_DECL_HPP
