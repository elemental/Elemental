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
#ifndef ELEMENTAL_IMPORTS_BLAS_HPP
#define ELEMENTAL_IMPORTS_BLAS_HPP 1

namespace elem {
namespace blas {

//
// NOTE: templated routines are custom and not wrappers
//

//----------------------------------------------------------------//
// Level 1 BLAS                                                   //
//----------------------------------------------------------------//
void Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy );
void Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy );
void Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy );
void Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy );
template<typename T>
void Axpy( int n, T alpha, const T* x, int incx, T* y, int incy );

float Dot( int n, const float* x, int incx, const float* y, int incy );
double Dot( int n, const double* x, int incx, const double* y, int incy );
scomplex Dot( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dot( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dot( int n, const T* x, int incx, const T* y, int incy );

float Dotc
( int n, const float* x, int incx, const float* y, int incy );
double Dotc
( int n, const double* x, int incx, const double* y, int incy );
scomplex Dotc
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotc
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dotc( int n, const T* x, int incx, const T* y, int incy );

float Dotu
( int n, const float* x, int incx, const float* y, int incy );
double Dotu
( int n, const double* x, int incx, const double* y, int incy );
scomplex Dotu
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotu
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dotu( int n, const T* x, int incx, const T* y, int incy );

float Nrm2( int n, const float* x, int incx );
double Nrm2( int n, const double* x, int incx );
float Nrm2( int n, const scomplex* x, int incx );
double Nrm2( int n, const dcomplex* x, int incx );
template<typename F> F Nrm2( int n, const F* x, int incx );

void Scal( int n, float alpha, float* x, int incx );
void Scal( int n, double alpha, double* x, int incx );
void Scal( int n, scomplex alpha, scomplex* x, int incx );
void Scal( int n, dcomplex alpha, dcomplex* x, int incx );
template<typename F> void Scal( int n, F alpha, F* x, int incx );
            
//----------------------------------------------------------------//
// Level 2 BLAS                                                   //
//----------------------------------------------------------------//
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
void Gemv
( char trans, int m, int n,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy );

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
void Ger
( char trans, int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
  T beta,        T* A, int lda );

void Gerc
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Gerc
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Gerc
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Gerc
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
template<typename T>
void Gerc
( char trans, int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
  T beta,        T* A, int lda );

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
void Geru
( char trans, int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
  T beta,        T* A, int lda );

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
void Hemv
( char uplo, int m,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy );

void Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda );
void Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda );
void Her
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda );
void Her
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda );
template<typename T>
void Hemv( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda );

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
void Her2
( char uplo, int m,
  T alpha, const T* x, int incx, const T* y, int incy, 
                 T* A, int lda );

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
void Symv
( char uplo, int m,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy );

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
void Syr( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda );

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
void Syr2
( char uplo, int m,
  T alpha, const T* x, int incx, const T* y, int incy,
                 T* A, int lda );

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
template<typename T>
void Trmv
( char uplo, char trans, char diag, int m,
  const T* A, int lda, T* x, int incx );

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
template<typename T>
void Trsv
( char uplo, char trans, char diag, int m,
  const T* A, int lda, T* x, int incx );

//----------------------------------------------------------------//
// Level 3 BLAS                                                   //
//----------------------------------------------------------------//
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
void Gemm
( char transA, char transB, int m, int n, int k,
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
template<typename T>
void Hemm
( char side, char uplo, int m, int n,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

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
  scomplex beta,        scomplex* C, int ldc );
void Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Her2k
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, float beta, float* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, double beta, double* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Herk
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda,
  T beta,        T* C, int ldc );

// NOTE: This is the only non-standard naming convention of an existing BLAS
//       routine. The routines for forming U := U' U and L := L L' are in
//       LAPACK and are called ?lauum. I am instead labeling it Hetrmm to 
//       match the BLAS naming conventions, and it stands for 
//       'HErmitian TRiangular Matrix-Matrix multiplication'
void Hetrmm( char uplo, int n, float* A, int lda );
void Hetrmm( char uplo, int n, double* A, int lda );
void Hetrmm( char uplo, int n, scomplex* A, int lda );
void Hetrmm( char uplo, int n, dcomplex* A, int lda );
template<typename T> void Hetrmm( char uplo, int n, T* A, int lda );

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
template<typename T>
void Symm
( char side, char uplo, int m, int n,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

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
void Syr2k
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
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
template<typename T>
void Syrk
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda,
  T beta,        T* C, int ldc );

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
template<typename T>
void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  T alpha, const T* A, int lda, T* B, int ldb );

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
template<typename T>
void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  T alpha, const T* A, int lda, T* B, int ldb );

} // namespace blas
} // namespace elem

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

void LAPACK(slauum)( char* uplo, int* n, float* A, int* lda, int* info );
void LAPACK(dlauum)( char* uplo, int* n, double* A, int* lda, int* info );
void LAPACK(clauum)
( char* uplo, int* n, elem::scomplex* A, int* lda, int* info );
void LAPACK(zlauum)
( char* uplo, int* n, elem::dcomplex* A, int* lda, int* info );

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

//
// Templated wrappers
//

namespace elem {
namespace blas {

//
// Level 1 BLAS
//

template<typename T>
inline void Axpy
( int n, T alpha, const T* x, int incx, T* y, int incy )
{
    for( int i=0; i<n; ++i )
        y[i*incy] += alpha*x[i*incx];
}

template<typename T>
inline T Dot( int n, const T* x, int incx, const T* y, int incy )
{ Dotc( n, x, incx, y, incy ); }

template<typename T>
inline T Dotc( int n, const T* x, int incx, const T* y, int incy )
{
    T alpha = 0;
    for( int i=0; i<n; ++i )
        alpha += Conj(x[i*incx])*y[i*incy];
    return alpha;
}

template<typename T>
inline T Dotu( int n, const T* x, int incx, const T* y, int incy )
{
    T alpha = 0;
    for( int i=0; i<n; ++i )
        alpha += x[i*incx]*y[i*incy];
    return alpha;
}

// TODO: templated Nrm2

template<typename T>
void Scal( int n, T alpha, T* x, int incx )
{
    for( int i=0; i<n; ++i )
        x[i*incx] *= alpha;
}

// 
// Level 2 BLAS
//

template<typename T>
void Gemv
( char trans, int m, int n,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy )
{
    if( trans == 'N' )
    {
        if( m > 0 && n == 0 && beta == 0 )
        {
            for( int i=0; i<m; ++i )
                y[i*incy] = 0;   
            return;
        }
        Scal( m, beta, y, incy );
        for( int i=0; i<m; ++i ) 
            for( int j=0; j<n; ++j )
                y[i*incy] += alpha*A[i+j*lda]*x[j*incx];
    }
    else if( trans == 'T' ) 
    {
        if( n > 0 && m == 0 && beta == 0 )
        {
            for( int i=0; i<n; ++i )
                y[i*incy] = 0;   
            return;
        }
        Scal( n, beta, y, incy );
        for( int i=0; i<n; ++i ) 
            for( int j=0; j<m; ++j )
                y[i*incy] += alpha*A[j+i*lda]*x[j*incx];
    }
    else
    {
        if( n > 0 && m == 0 && beta == 0 )
        {
            for( int i=0; i<n; ++i )
                y[i*incy] = 0;   
            return;
        }
        Scal( n, beta, y, incy );
        for( int i=0; i<n; ++i ) 
            for( int j=0; j<m; ++j )
                y[i*incy] += alpha*Conj(A[j+i*lda])*x[j*incx];
    }
}

// TODO: templated Ger
// TODO: templated Gerc
// TODO: templated Geru
// TODO: templated Hemv
// TODO: templated Her
// TODO: templated Her2
// TODO: templated Symv
// TODO: templated Syr
// TODO: templated Syr2
// TODO: templated Trmv
// TODO: templated Trsv

//
// Level 3 BLAS
//

template<typename T>
void Gemm
( char transA, char transB, int m, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc )
{
    if( m > 0 && n > 0 && k == 0 && beta == 0 )
    {
        for( int j=0; j<n; ++j )
            for( int i=0; i<m; ++i )
                C[i+j*ldc] = 0;
        return;
    }

    // Scale C
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            C[i+j*ldc] *= beta;

    // Naive implementation
    if( transA == 'N' && transB == 'N' )
    {
        // C := alpha A B + C
        for( int j=0; j<n; ++j )
            for( int i=0; i<m; ++i )
                for( int l=0; l<k; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*B[l+j*ldb];
    }
    else if( transA == 'N' )
    {
        if( transB == 'T' )
        {
            // C := alpha A B^T + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[i+l*lda]*B[j+l*ldb];
        }
        else
        {
            // C := alpha A B^H + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[i+l*lda]*Conj(B[j+l*ldb]);
        }
    }
    else if( transB == 'N' )
    {
        if( transA == 'T' )
        {
            // C := alpha A^T B + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*B[l+j*ldb];
        }
        else
        {
            // C := alpha A^H B + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*B[l+j*ldb];
        }
    }
    else
    {
        if( transA == 'T' && transB == 'T' )
        {
            // C := alpha A^T B^T + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*B[j+l*ldb];
        }
        else if( transA == 'T' )
        {
            // C := alpha A^T B^H + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*Conj(B[j+l*ldb]);
        }
        else if( transB == 'T' )
        {
            // C := alpha A^H B^T + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*B[j+l*ldb];
        }
        else
        {
            // C := alpha A^H B^H + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*Conj(B[j+l*ldb]);
        }
    }
}

// TODO: templated Hemm
// TODO: templated Her2k
// TODO: templated Herk
// TODO: templated Hetrmm
// TODO: templated Symm
// TODO: templated Syr2k
// TODO: templated Syrk
// TODO: templated Trmm
// TODO: templated Trsm

} // namespace blas
} // namespace elem

#endif /* ELEMENTAL_IMPORTS_BLAS_HPP */

