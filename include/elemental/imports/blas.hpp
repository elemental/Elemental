/*
   Copyright (c) 2009-2011, Jack Poulson
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

namespace elemental {
namespace blas {
//----------------------------------------------------------------//
// Level 1 BLAS                                                   //
//----------------------------------------------------------------//
void Axpy
( int n, int alpha, const int* x, int incx, int* y, int incy );
void Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy );
void Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy );
void Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy );
void Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy );

float Dot( int n, const float* x, int incx, const float* y, int incy );
double Dot( int n, const double* x, int incx, const double* y, int incy );
scomplex Dot( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dot( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );

float Dotc
( int n, const float* x, int incx, const float* y, int incy );
double Dotc
( int n, const double* x, int incx, const double* y, int incy );
scomplex Dotc
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotc
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );

float Dotu
( int n, const float* x, int incx, const float* y, int incy );
double Dotu
( int n, const double* x, int incx, const double* y, int incy );
scomplex Dotu
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotu
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );

float Nrm2( int n, const float* x, int incx );
double Nrm2( int n, const double* x, int incx );
float Nrm2( int n, const scomplex* x, int incx );
double Nrm2( int n, const dcomplex* x, int incx );

void Scal( int n, float alpha, float* x, int incx );
void Scal( int n, double alpha, double* x, int incx );
void Scal( int n, scomplex alpha, scomplex* x, int incx );
void Scal( int n, dcomplex alpha, dcomplex* x, int incx );
            
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

// NOTE: This is the only non-standard naming convention of an existing BLAS
//       routine. The routines for forming U := U' U and L := L L' are in
//       LAPACK and are called ?lauum. I am instead labeling it Hetrmm to 
//       match the BLAS naming conventions, and it stands for 
//       'HErmitian TRiangular Matrix-Matrix multiplication'
void Hetrmm( char uplo, int n, float* A, int lda );
void Hetrmm( char uplo, int n, double* A, int lda );
void Hetrmm( char uplo, int n, scomplex* A, int lda );
void Hetrmm( char uplo, int n, dcomplex* A, int lda );

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

} // blas
} // elemental

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
        elemental::scomplex* y, const int* incy );
void BLAS(zaxpy)
( const int* n, 
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
        elemental::dcomplex* y, const int* incy );

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
( const int* n, const elemental::scomplex* x, const int* incx );
double BLAS(dznrm2)
( const int* n, const elemental::dcomplex* x, const int* incx );

void BLAS(sscal)
( const int* n, const float* alpha, float* x, const int* incx );
void BLAS(dscal)
( const int* n, const double* alpha, double* x, const int* incx );
void BLAS(cscal)
( const int* n, const elemental::scomplex* alpha, elemental::scomplex* x, 
  const int* incx );
void BLAS(zscal)
( const int* n, const elemental::dcomplex* alpha, elemental::dcomplex* x, 
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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* beta,        
        elemental::scomplex* y, const int* incy );
void BLAS(zgemv)
( const char* trans, const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* y, const int* incy );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* y, const int* incy,
        elemental::scomplex* A, const int* lda  );
void BLAS(zgerc)
( const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* y, const int* incy,
        elemental::dcomplex* A, const int* lda  );

void BLAS(cgeru)
( const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* y, const int* incy,
        elemental::scomplex* A, const int* lda  );
void BLAS(zgeru)
( const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* y, const int* incy,
        elemental::dcomplex* A, const int* lda  );

void BLAS(chemv)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* beta,        
        elemental::scomplex* y, const int* incy );
void BLAS(zhemv)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* y, const int* incy );

void BLAS(cher)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
        elemental::scomplex* A, const int* lda  );
void BLAS(zher)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
        elemental::dcomplex* A, const int* lda  );

void BLAS(cher2)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* y, const int* incy,
        elemental::scomplex* A, const int* lda  );
void BLAS(zher2)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* y, const int* incy,
        elemental::dcomplex* A, const int* lda  );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* beta,        
        elemental::scomplex* y, const int* incy );
// 'zsymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void LAPACK(zsymv)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* y, const int* incy );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
        elemental::scomplex* A, const int* lda  );
// 'zsyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void LAPACK(zsyr)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
        elemental::dcomplex* A, const int* lda  );

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
  const elemental::scomplex* A, const int* lda, 
        elemental::scomplex* x, const int* incx );
void BLAS(ztrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elemental::dcomplex* A, const int* lda, 
        elemental::dcomplex* x, const int* incx );

void BLAS(strsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const float* A, const int* lda, float* x, const int* incx );
void BLAS(dtrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const double* A, const int* lda, double* x, const int* incx );
void BLAS(ctrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elemental::scomplex* A, const int* lda, 
        elemental::scomplex* x, const int* incx );
void BLAS(ztrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elemental::dcomplex* A, const int* lda, 
        elemental::dcomplex* x, const int* incx );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );
void BLAS(zgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

void BLAS(chemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );
void BLAS(zhemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

void BLAS(cher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );
void BLAS(zher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

void BLAS(cherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );
void BLAS(zherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

void LAPACK(slauum)( char* uplo, int* n, float* A, int* lda, int* info );
void LAPACK(dlauum)( char* uplo, int* n, double* A, int* lda, int* info );
void LAPACK(clauum)
( char* uplo, int* n, elemental::scomplex* A, int* lda, int* info );
void LAPACK(zlauum)
( char* uplo, int* n, elemental::dcomplex* A, int* lda, int* info );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );
void BLAS(zsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );
void BLAS(zsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );
void BLAS(zsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
        elemental::scomplex* B, const int* ldb );
void BLAS(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
        elemental::dcomplex* B, const int* ldb );

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
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
        elemental::scomplex* B, const int* ldb );
void BLAS(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
        elemental::dcomplex* B, const int* ldb );
} // extern "C"

#endif /* ELEMENTAL_IMPORTS_BLAS_HPP */

