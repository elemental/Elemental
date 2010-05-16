/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_WRAPPERS_BLAS_HPP
#define ELEMENTAL_WRAPPERS_BLAS_HPP 1

#include "elemental/environment.hpp"

#ifdef FUNDERSCORE
#define C2F(name) name ## _
#else
#define C2F(name) name
#endif

namespace elemental {
namespace wrappers {
namespace blas {
//----------------------------------------------------------------//
// Level 1 BLAS                                                   //
//----------------------------------------------------------------//
void
Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy );

void
Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy );

#ifndef WITHOUT_COMPLEX
void
Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy );

void
Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy );
#endif

float
Dot
( int n, const float* x, int incx, const float* y, int incy );

double
Dot
( int n, const double* x, int incx, const double* y, int incy );

#ifndef WITHOUT_COMPLEX
scomplex
Dot
( int n, const scomplex* x, int incx, const scomplex* y, int incy );

dcomplex
Dot
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
#endif

float
Dotc
( int n, const float* x, int incx, const float* y, int incy );

double
Dotc
( int n, const double* x, int incx, const double* y, int incy );

#ifndef WITHOUT_COMPLEX
scomplex
Dotc
( int n, const scomplex* x, int incx, const scomplex* y, int incy );

dcomplex
Dotc
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
#endif

float
Dotu
( int n, const float* x, int incx, const float* y, int incy );

double
Dotu
( int n, const double* x, int incx, const double* y, int incy );

#ifndef WITHOUT_COMPLEX
scomplex
Dotu
( int n, const scomplex* x, int incx, const scomplex* y, int incy );

dcomplex
Dotu
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
#endif

float
Nrm2
( int n, const float* x, int incx );

double
Nrm2
( int n, const double* x, int incx );

#ifndef WITHOUT_COMPLEX
float
Nrm2
( int n, const scomplex* x, int incx );

double
Nrm2
( int n, const dcomplex* x, int incx );
#endif

void
Scal
( int n, float alpha, float* x, int incx );

void
Scal
( int n, double alpha, double* x, int incx );

#ifndef WITHOUT_COMPLEX
void
Scal
( int n, scomplex alpha, scomplex* x, int incx );

void
Scal
( int n, dcomplex alpha, dcomplex* x, int incx );
#endif
            
//----------------------------------------------------------------//
// Level 2 BLAS                                                   //
//----------------------------------------------------------------//
void
Gemv
( char trans, int m, int n,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );

void
Gemv
( char trans, int m, int n,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );

#ifndef WITHOUT_COMPLEX
void
Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );

void
Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );
#endif

void
Ger
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );

void
Ger
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Ger
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );

void
Ger
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
#endif

void
Gerc
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );

void
Gerc
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Gerc
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );

void
Gerc
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
#endif

void
Geru
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );

void
Geru
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Geru
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );

void
Geru
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
#endif

void
Hemv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );

void
Hemv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );

#ifndef WITHOUT_COMPLEX
void
Hemv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );

void
Hemv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );
#endif

void
Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda );

void
Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Her
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda );

void
Her
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda );
#endif

void
Her2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );

void
Her2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Her2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );

void
Her2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
#endif

void
Symv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );

void
Symv
( char uplo, int m, 
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );

#ifndef WITHOUT_COMPLEX
void
Symv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );

void
Symv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );
#endif

void
Syr
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda );

void
Syr
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Syr
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda ); 

void
Syr
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda );
#endif

void
Syr2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );

void
Syr2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Syr2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );

void
Syr2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
#endif

void
Trmv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx );

void
Trmv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx );

#ifndef WITHOUT_COMPLEX
void
Trmv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx );

void
Trmv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx );
#endif

void
Trsv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx );

void
Trsv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx );

#ifndef WITHOUT_COMPLEX
void
Trsv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx );

void
Trsv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx );
#endif

//----------------------------------------------------------------//
// Level 3 BLAS                                                   //
//----------------------------------------------------------------//
void
Gemm
( char transA, char transB, int m, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );

void
Gemm
( char transA, char transB, int m, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );

#ifndef WITHOUT_COMPLEX
void
Gemm
( char transA, char transB, int m, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );

void
Gemm
( char transA, char transB, int m, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
#endif

void
Hemm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );

void
Hemm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );

#ifndef WITHOUT_COMPLEX
void
Hemm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );

void
Hemm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
#endif

void
Her2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );

void
Her2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );

#ifndef WITHOUT_COMPLEX
void
Her2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );

void
Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
#endif

void
Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, float beta, float* C, int ldc );

void
Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, double beta, double* C, int ldc );

#ifndef WITHOUT_COMPLEX
void
Herk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc );

void
Herk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc );
#endif

void
Symm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );

void
Symm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );

#ifndef WITHOUT_COMPLEX
void
Symm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );

void
Symm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
#endif

void
Syr2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );

void
Syr2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );

#ifndef WITHOUT_COMPLEX
void
Syr2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );

void
Syr2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
#endif

void
Syrk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc );

void
Syrk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc );

#ifndef WITHOUT_COMPLEX
void
Syrk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc );

void
Syrk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc );
#endif

void
Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* X, int ldb );

void
Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* X, int ldb );

#ifndef WITHOUT_COMPLEX
void
Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* X, int ldb );

void
Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* X, int ldb );
#endif

void
Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* X, int ldb );

void
Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* X, int ldb );

#ifndef WITHOUT_COMPLEX
void
Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* X, int ldb );

void
Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* X, int ldb );
#endif
} // blas
} // wrappers
} // elemental

extern "C" {
//------------------------------------------------------------------------//
// Level 1 BLAS                                                           //
//------------------------------------------------------------------------//
void C2F(saxpy)
( const int* n, const float* alpha, const float* x, const int* incx,
                                          float* y, const int* incy );

void C2F(daxpy)
( const int* n, const double* alpha, const double* x, const int* incx,
                                           double* y, const int* incy );

#ifndef WITHOUT_COMPLEX
void C2F(caxpy)
( const int* n, 
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
        elemental::scomplex* y, const int* incy );
    
void C2F(zaxpy)
( const int* n, 
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
        elemental::dcomplex* y, const int* incy );
#endif

float C2F(sdot)
( const int* n, const float* x, const int* incx,
                const float* y, const int* incy );

double C2F(ddot)
( const int* n, const double* x, const int* incx,
                const double* y, const int* incy );

#ifndef WITHOUT_COMPLEX
elemental::scomplex C2F(cdotu)
( const int* n, const elemental::scomplex* x, const int* incx,
                const elemental::scomplex* y, const int* incy );

elemental::dcomplex C2F(zdotu)
( const int* n, const elemental::dcomplex* x, const int* incx,
                const elemental::dcomplex* y, const int* incy );

elemental::scomplex C2F(cdotc)
( const int* n, const elemental::scomplex* x, const int* incx,
                const elemental::scomplex* y, const int* incy );

elemental::dcomplex C2F(zdotc)
( const int* n, const elemental::dcomplex* x, const int* incx,
                const elemental::dcomplex* y, const int* incy );
#endif

float C2F(snrm2)
( const int* n, const float* x, const int* incx );

double C2F(dnrm2)
( const int* n, const double* x, const int* incx );

#ifndef WITHOUT_COMPLEX
float C2F(scnrm2)
( const int* n, const elemental::scomplex* x, const int* incx );

double C2F(dznrm2)
( const int* n, const elemental::dcomplex* x, const int* incx );
#endif

void C2F(sscal)
( const int* n, const float* alpha, float* x, const int* incx );

void C2F(dscal)
( const int* n, const double* alpha, double* x, const int* incx );
    
#ifndef WITHOUT_COMPLEX
void C2F(cscal)
( const int* n, const elemental::scomplex* alpha, elemental::scomplex* x, 
  const int* incx );
    
void C2F(zscal)
( const int* n, const elemental::dcomplex* alpha, elemental::dcomplex* x, 
  const int* incx );
#endif

//------------------------------------------------------------------------//
// Level 2 BLAS                                                           //
//------------------------------------------------------------------------//
void C2F(sgemv)
( const char* trans, const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                      const float* x, const int* incx,
  const float* beta,        float* y, const int* incy );

void C2F(dgemv)
( const char* trans, const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                       const double* x, const int* incx,
  const double* beta,        double* y, const int* incy );

#ifndef WITHOUT_COMPLEX
void C2F(cgemv)
( const char* trans, const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* beta,        
        elemental::scomplex* y, const int* incy );

void C2F(zgemv)
( const char* trans, const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* y, const int* incy );
#endif

void C2F(sger)
( const int* m, const int* n,
  const float* alpha, const float* x, const int* incx,
                      const float* y, const int* incy,
                            float* A, const int* lda  );

void C2F(dger)
( const int* m, const int* n,
  const double* alpha, const double* x, const int* incx,
                       const double* y, const int* incy,
                             double* A, const int* lda  );

#ifndef WITHOUT_COMPLEX
void C2F(cgerc)
( const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* y, const int* incy,
        elemental::scomplex* A, const int* lda  );

void C2F(zgerc)
( const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* y, const int* incy,
        elemental::dcomplex* A, const int* lda  );

void C2F(cgeru)
( const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* y, const int* incy,
        elemental::scomplex* A, const int* lda  );

void C2F(zgeru)
( const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* y, const int* incy,
        elemental::dcomplex* A, const int* lda  );

void C2F(chemv)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* beta,        
        elemental::scomplex* y, const int* incy );

void C2F(zhemv)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* y, const int* incy );

void C2F(cher)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
        elemental::scomplex* A, const int* lda  );

void C2F(zher)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
        elemental::dcomplex* A, const int* lda  );

void C2F(cher2)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* y, const int* incy,
        elemental::scomplex* A, const int* lda  );

void C2F(zher2)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* y, const int* incy,
        elemental::dcomplex* A, const int* lda  );
#endif

void C2F(ssymv)
( const char* uplo, const int* m,
  const float* alpha, const float* A, const int* lda,
                      const float* x, const int* incx,
  const float* beta,        float* y, const int* incy );

void C2F(dsymv)
( const char* uplo, const int* m,
  const double* alpha, const double* A, const int* lda,
                       const double* x, const int* incx,
  const double* beta,        double* y, const int* incy );

#ifndef WITHOUT_COMPLEX
// Recall that 'csymv' is an auxiliary LAPACK routine
void C2F(csymv)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* x, const int* incx,
  const elemental::scomplex* beta,        
        elemental::scomplex* y, const int* incy );

// Recall that 'zsymv' is an auxiliary LAPACK routine
void C2F(zsymv)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* x, const int* incx,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* y, const int* incy );
#endif

void C2F(ssyr)
( const char* uplo, const int* m,
  const float* alpha, const float* x, const int* incx,
                            float* A, const int* lda  );

void C2F(dsyr)
( const char* uplo, const int* m,
  const double* alpha, const double* x, const int* incx,
                             double* A, const int* lda  );

#ifndef WITHOUT_COMPLEX
// Recall that 'csyr' is an auxiliary LAPACK routine
void C2F(csyr)
( const char* uplo, const int* m,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* x, const int* incx,
        elemental::scomplex* A, const int* lda  );

// Recall that 'zsyr' is an auxiliary LAPACK routine
void C2F(zsyr)
( const char* uplo, const int* m,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* x, const int* incx,
        elemental::dcomplex* A, const int* lda  );
#endif

void C2F(ssyr2)
( const char* uplo, const int* m,
  const float* alpha, const float* x, const int* incx,
                      const float* y, const int* incy,
                            float* A, const int* lda  );

void C2F(dsyr2)
( const char* uplo, const int* m,
  const double* alpha, const double* x, const int* incx,
                       const double* y, const int* incy,
                             double* A, const int* lda  );

void C2F(strmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const float* A, const int* lda, float* x, const int* incx );

void C2F(dtrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const double* A, const int* lda, double* x, const int* incx );

#ifndef WITHOUT_COMPLEX
void C2F(ctrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elemental::scomplex* A, const int* lda, 
        elemental::scomplex* x, const int* incx );

void C2F(ztrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elemental::dcomplex* A, const int* lda, 
        elemental::dcomplex* x, const int* incx );
#endif

void C2F(strsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const float* A, const int* lda, float* x, const int* incx );

void C2F(dtrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const double* A, const int* lda, double* x, const int* incx );

#ifndef WITHOUT_COMPLEX
void C2F(ctrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elemental::scomplex* A, const int* lda, 
        elemental::scomplex* x, const int* incx );

void C2F(ztrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const elemental::dcomplex* A, const int* lda, 
        elemental::dcomplex* x, const int* incx );
#endif

//------------------------------------------------------------------------//
// Level 3 BLAS                                                           //
//------------------------------------------------------------------------//
void C2F(sgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );
    
void C2F(dgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );
    
#ifndef WITHOUT_COMPLEX
void C2F(cgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );

void C2F(zgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

void C2F(chemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );

void C2F(zhemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

void C2F(cher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );

void C2F(zher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );

void C2F(cherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );

void C2F(zherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );
#endif

void C2F(ssymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );

void C2F(dsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );

#ifndef WITHOUT_COMPLEX
void C2F(csymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );

void C2F(zsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );
#endif

void C2F(ssyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
                      const float* B, const int* ldb,
  const float* beta,        float* C, const int* ldc );

void C2F(dsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
                       const double* B, const int* ldb,
  const double* beta,        double* C, const int* ldc );

#ifndef WITHOUT_COMPLEX
void C2F(csyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* B, const int* ldb,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );

void C2F(zsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );
#endif

void C2F(ssyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const float* alpha, const float* A, const int* lda,
  const float* beta,        float* C, const int* ldc );

void C2F(dsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const double* alpha, const double* A, const int* lda,
  const double* beta,        double* C, const int* ldc );

#ifndef WITHOUT_COMPLEX
void C2F(csyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
  const elemental::scomplex* beta,        
        elemental::scomplex* C, const int* ldc );

void C2F(zsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* beta,        
        elemental::dcomplex* C, const int* ldc );
#endif

void C2F(strmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, 
  const float* alpha, const float* A, const int* lda,
                            float* B, const int* ldb );

void C2F(dtrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                             double* B, const int* ldb );

#ifndef WITHOUT_COMPLEX
void C2F(ctrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
        elemental::scomplex* B, const int* ldb );

void C2F(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
        elemental::dcomplex* B, const int* ldb );
#endif

void C2F(strsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n, 
  const float* alpha, const float* A, const int* lda,
                            float* B, const int* ldb );

void C2F(dtrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const double* alpha, const double* A, const int* lda,
                             double* B, const int* ldb );

#ifndef WITHOUT_COMPLEX
void C2F(ctrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const elemental::scomplex* alpha, 
  const elemental::scomplex* A, const int* lda,
        elemental::scomplex* B, const int* ldb );

void C2F(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const elemental::dcomplex* alpha, 
  const elemental::dcomplex* A, const int* lda,
        elemental::dcomplex* B, const int* ldb );
#endif
} // extern "C"

//----------------------------------------------------------------------------//
// Implementations begin here                                                 //
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Level 1 BLAS                                                               //
//----------------------------------------------------------------------------//
inline void
elemental::wrappers::blas::Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy )
{ C2F(saxpy)( &n, &alpha, x, &incx, y, &incy ); }

inline void
elemental::wrappers::blas::Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy )
{ C2F(daxpy)( &n, &alpha, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy )
{ C2F(caxpy)( &n, &alpha, x, &incx, y, &incy ); }

inline void
elemental::wrappers::blas::Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy )
{ C2F(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }
#endif

inline float
elemental::wrappers::blas::Dot
( int n, const float* x, int incx, const float* y, int incy )
{ return C2F(sdot)( &n, x, &incx, y, &incy ); }

inline double
elemental::wrappers::blas::Dot
( int n, const double* x, int incx, const double* y, int incy )
{ return C2F(ddot)( &n, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline elemental::scomplex
elemental::wrappers::blas::Dot
( int n, const elemental::scomplex* x, int incx,
         const elemental::scomplex* y, int incy )
{ return C2F(cdotc)( &n, x, &incx, y, &incy ); }

inline elemental::dcomplex
elemental::wrappers::blas::Dot
( int n, const elemental::dcomplex* x, int incx,
         const elemental::dcomplex* y, int incy )
{ return C2F(zdotc)( &n, x, &incx, y, &incy ); }
#endif

inline float
elemental::wrappers::blas::Dotc
( int n, const float* x, int incx, const float* y, int incy )
{ return C2F(sdot)( &n, x, &incx, y, &incy ); }

inline double
elemental::wrappers::blas::Dotc
( int n, const double* x, int incx, const double* y, int incy )
{ return C2F(ddot)( &n, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline elemental::scomplex
elemental::wrappers::blas::Dotc
( int n, const elemental::scomplex* x, int incx,
         const elemental::scomplex* y, int incy )
{ return C2F(cdotc)( &n, x, &incx, y, &incy ); }

inline elemental::dcomplex
elemental::wrappers::blas::Dotc
( int n, const elemental::dcomplex* x, int incx,
         const elemental::dcomplex* y, int incy )
{ return C2F(zdotc)( &n, x, &incx, y, &incy ); }
#endif

inline float
elemental::wrappers::blas::Dotu
( int n, const float* x, int incx, const float* y, int incy )
{ return C2F(sdot)( &n, x, &incx, y, &incy ); }

inline double
elemental::wrappers::blas::Dotu
( int n, const double* x, int incx, const double* y, int incy )
{ return C2F(ddot)( &n, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline elemental::scomplex
elemental::wrappers::blas::Dotu
( int n, const elemental::scomplex* x, int incx,
         const elemental::scomplex* y, int incy )
{ return C2F(cdotu)( &n, x, &incx, y, &incy ); }

inline elemental::dcomplex
elemental::wrappers::blas::Dotu
( int n, const elemental::dcomplex* x, int incx,
         const elemental::dcomplex* y, int incy )
{ return C2F(zdotu)( &n, x, &incx, y, &incy ); }
#endif

inline float
elemental::wrappers::blas::Nrm2
( int n, const float* x, int incx )
{ return C2F(snrm2)( &n, x, &incx ); }

inline double
elemental::wrappers::blas::Nrm2
( int n, const double* x, int incx )
{ return C2F(dnrm2)( &n, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline float
elemental::wrappers::blas::Nrm2
( int n, const scomplex* x, int incx )
{ return C2F(scnrm2)( &n, x, &incx ); }

inline double
elemental::wrappers::blas::Nrm2
( int n, const dcomplex* x, int incx )
{ return C2F(dznrm2)( &n, x, &incx ); }
#endif

inline void
elemental::wrappers::blas::Scal
( int n, float alpha, float* x, int incx )
{ C2F(sscal)( &n, &alpha, x, &incx ); }

inline void
elemental::wrappers::blas::Scal
( int n, double alpha, double* x, int incx )
{ C2F(dscal)( &n, &alpha, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Scal
( int n, scomplex alpha, scomplex* x, int incx )
{ C2F(cscal)( &n, &alpha, x, &incx ); }

inline void
elemental::wrappers::blas::Scal
( int n, dcomplex alpha, dcomplex* x, int incx )
{ C2F(zscal)( &n, &alpha, x, &incx ); }
#endif

//----------------------------------------------------------------------------//
// Level 2 BLAS                                                               //
//----------------------------------------------------------------------------//
inline void
elemental::wrappers::blas::Gemv
( char trans, int m, int n,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

inline void
elemental::wrappers::blas::Gemv
( char trans, int m, int n,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ C2F(cgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
elemental::wrappers::blas::Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ C2F(zgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }
#endif

inline void
elemental::wrappers::blas::Ger
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ C2F(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Ger
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda  )
{ C2F(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Ger
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ C2F(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Ger
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ C2F(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
elemental::wrappers::blas::Gerc
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ C2F(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Gerc
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ C2F(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Gerc
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ C2F(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Gerc
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ C2F(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
elemental::wrappers::blas::Geru
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ C2F(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Geru
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ C2F(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Geru
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ C2F(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Geru
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ C2F(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
elemental::wrappers::blas::Hemv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ C2F(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
elemental::wrappers::blas::Hemv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ C2F(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Hemv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{ C2F(chemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
elemental::wrappers::blas::Hemv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{ C2F(zhemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }
#endif

inline void
elemental::wrappers::blas::Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda )
{ C2F(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

inline void
elemental::wrappers::blas::Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ C2F(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Her
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{ C2F(cher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

inline void
elemental::wrappers::blas::Her
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{ C2F(zher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }
#endif

inline void
elemental::wrappers::blas::Her2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ C2F(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Her2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ C2F(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Her2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda )
{ C2F(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Her2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda )
{ C2F(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
elemental::wrappers::blas::Symv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy )
{ C2F(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
elemental::wrappers::blas::Symv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy )
{ C2F(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Symv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    C2F(csymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

inline void
elemental::wrappers::blas::Symv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    C2F(zsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}
#endif

inline void
elemental::wrappers::blas::Syr
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda  )
{ C2F(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

inline void
elemental::wrappers::blas::Syr
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda )
{ C2F(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Syr
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    C2F(csyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

inline void
elemental::wrappers::blas::Syr
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    C2F(zsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}
#endif

inline void
elemental::wrappers::blas::Syr2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda )
{ C2F(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
elemental::wrappers::blas::Syr2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda )
{ C2F(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Syr2
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
    C2F(csyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}

inline void
elemental::wrappers::blas::Syr2
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
    C2F(zsyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &lda );
}
#endif

inline void
elemental::wrappers::blas::Trmv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ C2F(strmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
elemental::wrappers::blas::Trmv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ C2F(dtrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Trmv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ C2F(ctrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
elemental::wrappers::blas::Trmv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ C2F(ztrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }
#endif

inline void
elemental::wrappers::blas::Trsv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx )
{ C2F(strsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
elemental::wrappers::blas::Trsv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx )
{ C2F(dtrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Trsv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx )
{ C2F(ctrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
elemental::wrappers::blas::Trsv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx )
{ C2F(ztrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }
#endif

//----------------------------------------------------------------------------//
// Level 3 BLAS                                                               //
//----------------------------------------------------------------------------//
inline void
elemental::wrappers::blas::Gemm
( char transA, char transB, int m, int n, int k, 
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    C2F(sgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Gemm
( char transA, char transB,
  int m, int n, int k, 
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    C2F(dgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Gemm
( char transA, char transB, int m, int n, int k, 
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    C2F(cgemm)( &transA, &transB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Gemm
( char transA, char transB, int m, int n, int k, 
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    C2F(zgemm)( &transA, &transB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
elemental::wrappers::blas::Hemm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    C2F(ssymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Hemm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    C2F(dsymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Hemm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    C2F(chemm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Hemm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    C2F(zhemm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
elemental::wrappers::blas::Her2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(ssyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Her2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(dsyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Her2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    C2F(cher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    C2F(zher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
elemental::wrappers::blas::Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(ssyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(dsyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Herk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ C2F(cherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

inline void
elemental::wrappers::blas::Herk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ C2F(zherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }
#endif

inline void
elemental::wrappers::blas::Symm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    C2F(ssymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Symm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    C2F(dsymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Symm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    C2F(csymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Symm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    C2F(zsymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
elemental::wrappers::blas::Syr2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc )
{
    C2F(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Syr2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc )
{
    C2F(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Syr2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc )
{
    C2F(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
elemental::wrappers::blas::Syr2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc )
{
    C2F(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
elemental::wrappers::blas::Syrk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc )
{ C2F(ssyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

inline void
elemental::wrappers::blas::Syrk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc )
{ C2F(dsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Syrk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc )
{ C2F(csyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

inline void
elemental::wrappers::blas::Syrk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc )
{ C2F(zsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }
#endif

inline void
elemental::wrappers::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    C2F(strmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}

inline void
elemental::wrappers::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    C2F(dtrmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    C2F(ctrmm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}

inline void
elemental::wrappers::blas::Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    C2F(ztrmm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}
#endif

inline void
elemental::wrappers::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(strsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 

inline void
elemental::wrappers::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(dtrsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb )
{
    C2F(ctrsm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 

inline void
elemental::wrappers::blas::Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb )
{
    C2F(ztrsm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 
#endif

#endif /* ELEMENTAL_WRAPPERS_BLAS_HPP */

