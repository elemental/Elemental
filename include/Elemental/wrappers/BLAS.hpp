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

#include "Elemental/Environment.hpp"

#ifdef FUNDERSCORE
#define C2F(name) name ## _
#else
#define C2F(name) name
#endif

namespace Elemental {
namespace wrappers {
namespace BLAS {
//----------------------------------------------------------------//
// Level 1 BLAS                                                   //
//----------------------------------------------------------------//
void
Axpy
( const int n, const float alpha, 
  const float* x, const int incx,
        float* y, const int incy );

void
Axpy
( const int n, const double alpha, 
  const double* x, const int incx,
        double* y, const int incy );

#ifndef WITHOUT_COMPLEX
void
Axpy
( const int n, const scomplex alpha,
  const scomplex* x, const int incx,
        scomplex* y, const int incy );

void
Axpy
( const int n, const dcomplex alpha,
  const dcomplex* x, const int incx,
        dcomplex* y, const int incy );
#endif

float
Dot
( const int n, const float* x, const int incx,
               const float* y, const int incy );

double
Dot
( const int n, const double* x, const int incx,
               const double* y, const int incy );

#ifndef WITHOUT_COMPLEX
scomplex
Dot
( const int n, const scomplex* x, const int incx,
               const scomplex* y, const int incy );

dcomplex
Dot
( const int n, const dcomplex* x, const int incx,
               const dcomplex* y, const int incy );
#endif

float
Dotc
( const int n, const float* x, const int incx,
               const float* y, const int incy );

double
Dotc
( const int n, const double* x, const int incx,
               const double* y, const int incy );

#ifndef WITHOUT_COMPLEX
scomplex
Dotc
( const int n, const scomplex* x, const int incx,
               const scomplex* y, const int incy );

dcomplex
Dotc
( const int n, const dcomplex* x, const int incx,
               const dcomplex* y, const int incy );
#endif

float
Dotu
( const int n, const float* x, const int incx,
               const float* y, const int incy );

double
Dotu
( const int n, const double* x, const int incx,
               const double* y, const int incy );

#ifndef WITHOUT_COMPLEX
scomplex
Dotu
( const int n, const scomplex* x, const int incx,
               const scomplex* y, const int incy );

dcomplex
Dotu
( const int n, const dcomplex* x, const int incx,
               const dcomplex* y, const int incy );
#endif

float
Nrm2
( const int n, const float* x, const int incx );

double
Nrm2
( const int n, const double* x, const int incx );

#ifndef WITHOUT_COMPLEX
float
Nrm2
( const int n, const scomplex* x, const int incx );

double
Nrm2
( const int n, const dcomplex* x, const int incx );
#endif

void
Scal
( const int n, const float alpha, float* x, const int incx );

void
Scal
( const int n, const double alpha, double* x, const int incx );

#ifndef WITHOUT_COMPLEX
void
Scal
( const int n, const scomplex alpha, scomplex* x, const int incx );

void
Scal
( const int n, const dcomplex alpha, dcomplex* x, const int incx );
#endif
            
//----------------------------------------------------------------//
// Level 2 BLAS                                                   //
//----------------------------------------------------------------//
void
Gemv
( const char trans, const int m, const int n,
  const float alpha, const float* A, const int lda,
                     const float* x, const int incx,
  const float beta,        float* y, const int incy );

void
Gemv
( const char trans, const int m, const int n,
  const double alpha, const double* A, const int lda,
                      const double* x, const int incx,
  const double beta,        double* y, const int incy );

#ifndef WITHOUT_COMPLEX
void
Gemv
( const char trans, const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* x, const int incx,
  const scomplex beta,        scomplex* y, const int incy );

void
Gemv
( const char trans, const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* x, const int incx,
  const dcomplex beta,        dcomplex* y, const int incy );
#endif

void
Ger
( const int m, const int n,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  );

void
Ger
( const int m, const int n,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  );

#ifndef WITHOUT_COMPLEX
void
Ger
( const int m, const int n,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  );

void
Ger
( const int m, const int n,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  );
#endif

void
Gerc
( const int m, const int n,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  );

void
Gerc
( const int m, const int n,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  );

#ifndef WITHOUT_COMPLEX
void
Gerc
( const int m, const int n,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  );

void
Gerc
( const int m, const int n,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  );
#endif

void
Geru
( const int m, const int n,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  );

void
Geru
( const int m, const int n,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  );

#ifndef WITHOUT_COMPLEX
void
Geru
( const int m, const int n,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  );

void
Geru
( const int m, const int n,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  );
#endif

void
Hemv
( const char uplo, const int m,
  const float alpha, const float* A, const int lda,
                     const float* x, const int incx,
  const float beta,        float* y, const int incy );

void
Hemv
( const char uplo, const int m,
  const double alpha, const double* A, const int lda,
                      const double* x, const int incx,
  const double beta,        double* y, const int incy );

#ifndef WITHOUT_COMPLEX
void
Hemv
( const char uplo, const int m,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* x, const int incx,
  const scomplex beta,        scomplex* y, const int incy );

void
Hemv
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* x, const int incx,
  const dcomplex beta,        dcomplex* y, const int incy );
#endif

void
Her
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                           float* A, const int lda  );

void
Her
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                            double* A, const int lda  );

#ifndef WITHOUT_COMPLEX
void
Her
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                              scomplex* A, const int lda  );

void
Her
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                              dcomplex* A, const int lda  );
#endif

void
Her2
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  );

void
Her2
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  );

#ifndef WITHOUT_COMPLEX
void
Her2
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  );

void
Her2
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  );
#endif

void
Symv
( const char uplo, const int m,
  const float alpha, const float* A, const int lda,
                     const float* x, const int incx,
  const float beta,        float* y, const int incy );

void
Symv
( const char uplo, const int m, 
  const double alpha, const double* A, const int lda,
                      const double* x, const int incx,
  const double beta,        double* y, const int incy );

#ifndef WITHOUT_COMPLEX
void
Symv
( const char uplo, const int m,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* x, const int incx,
  const scomplex beta,        scomplex* y, const int incy );

void
Symv
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* x, const int incx,
  const dcomplex beta,        dcomplex* y, const int incy );
#endif

void
Syr
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                           float* A, const int lda  );

void
Syr
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                            double* A, const int lda  );

#ifndef WITHOUT_COMPLEX
void
Syr
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                              scomplex* A, const int lda  ); 

void
Syr
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                              dcomplex* A, const int lda  );
#endif

void
Syr2
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  );

void
Syr2
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  );

#ifndef WITHOUT_COMPLEX
void
Syr2
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  );

void
Syr2
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  );
#endif

void
Trmv
( const char uplo, const char trans, const char diag, const int m,
  const float* A, const int lda, float* x, const int incx         );

void
Trmv
( const char uplo, const char trans, const char diag, const int m,
  const double* A, const int lda, double* x, const int incx       );

#ifndef WITHOUT_COMPLEX
void
Trmv
( const char uplo, const char trans, const char diag, const int m,
  const scomplex* A, const int lda, scomplex* x, const int incx   );

void
Trmv
( const char uplo, const char trans, const char diag, const int m,
  const dcomplex* A, const int lda, dcomplex* x, const int incx   );
#endif

void
Trsv
( const char uplo, const char trans, const char diag, const int m,
  const float* A, const int lda, float* x, const int incx         );

void
Trsv
( const char uplo, const char trans, const char diag, const int m,
  const double* A, const int lda, double* x, const int incx       );

#ifndef WITHOUT_COMPLEX
void
Trsv
( const char uplo, const char trans, const char diag, const int m,
  const scomplex* A, const int lda, scomplex* x, const int incx   );

void
Trsv
( const char uplo, const char trans, const char diag, const int m,
  const dcomplex* A, const int lda, dcomplex* x, const int incx   );
#endif

//----------------------------------------------------------------//
// Level 3 BLAS                                                   //
//----------------------------------------------------------------//
void
Gemm
( const char transA, const char transB,
  const int m, const int n, const int k,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc );

void
Gemm
( const char transA, const char transB,
  const int m, const int n, const int k,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc );

#ifndef WITHOUT_COMPLEX
void
Gemm
( const char transA, const char transB,
  const int m, const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc );

void
Gemm
( const char transA, const char transB,
  const int m, const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc );
#endif

void
Hemm
( const char side, const char uplo,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc );

void
Hemm
( const char side, const char uplo,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc );

#ifndef WITHOUT_COMPLEX
void
Hemm
( const char side, const char uplo,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc );

void
Hemm
( const char side, const char uplo,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc );
#endif

void
Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc );

void
Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc );

#ifndef WITHOUT_COMPLEX
void
Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc );

void
Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc );
#endif

void
Herk
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
  const float beta,        float* C, const int ldc );

void
Herk
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
  const double beta,        double* C, const int ldc );

#ifndef WITHOUT_COMPLEX
void
Herk
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
  const scomplex beta,        scomplex* C, const int ldc );

void
Herk
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
  const dcomplex beta,        dcomplex* C, const int ldc );
#endif

void
Symm
( const char side, const char uplo,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc );

void
Symm
( const char side, const char uplo,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc );

#ifndef WITHOUT_COMPLEX
void
Symm
( const char side, const char uplo,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc );

void
Symm
( const char side, const char uplo,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc );
#endif

void
Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc );

void
Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc );

#ifndef WITHOUT_COMPLEX
void
Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc );

void
Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc );
#endif

void
Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
  const float beta,        float* C, const int ldc );

void
Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
  const double beta,        double* C, const int ldc );

#ifndef WITHOUT_COMPLEX
void
Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
  const scomplex beta,        scomplex* C, const int ldc );

void
Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
  const dcomplex beta,        dcomplex* C, const int ldc );
#endif

void
Trmm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                           float* X, const int ldb );

void
Trmm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                            double* X, const int ldb );

#ifndef WITHOUT_COMPLEX
void
Trmm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                              scomplex* X, const int ldb );

void
Trmm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                              dcomplex* X, const int ldb );
#endif

void
Trsm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                           float* X, const int ldb );

void
Trsm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                            double* X, const int ldb );

#ifndef WITHOUT_COMPLEX
void
Trsm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                              scomplex* X, const int ldb );

void
Trsm
( const char side,  const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                              dcomplex* X, const int ldb );
#endif
} // BLAS
} // wrappers
} // Elemental

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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* x, const int* incx,
        Elemental::scomplex* y, const int* incy );
    
void C2F(zaxpy)
( const int* n, 
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* x, const int* incx,
        Elemental::dcomplex* y, const int* incy );
#endif

float C2F(sdot)
( const int* n, const float* x, const int* incx,
                const float* y, const int* incy );

double C2F(ddot)
( const int* n, const double* x, const int* incx,
                const double* y, const int* incy );

#ifndef WITHOUT_COMPLEX
Elemental::scomplex C2F(cdotu)
( const int* n, const Elemental::scomplex* x, const int* incx,
                const Elemental::scomplex* y, const int* incy );

Elemental::dcomplex C2F(zdotu)
( const int* n, const Elemental::dcomplex* x, const int* incx,
                const Elemental::dcomplex* y, const int* incy );

Elemental::scomplex C2F(cdotc)
( const int* n, const Elemental::scomplex* x, const int* incx,
                const Elemental::scomplex* y, const int* incy );

Elemental::dcomplex C2F(zdotc)
( const int* n, const Elemental::dcomplex* x, const int* incx,
                const Elemental::dcomplex* y, const int* incy );
#endif

float C2F(snrm2)
( const int* n, const float* x, const int* incx );

double C2F(dnrm2)
( const int* n, const double* x, const int* incx );

#ifndef WITHOUT_COMPLEX
float C2F(scnrm2)
( const int* n, const Elemental::scomplex* x, const int* incx );

double C2F(dznrm2)
( const int* n, const Elemental::dcomplex* x, const int* incx );
#endif

void C2F(sscal)
( const int* n, const float* alpha, float* x, const int* incx );

void C2F(dscal)
( const int* n, const double* alpha, double* x, const int* incx );
    
#ifndef WITHOUT_COMPLEX
void C2F(cscal)
( const int* n, const Elemental::scomplex* alpha, Elemental::scomplex* x, 
  const int* incx );
    
void C2F(zscal)
( const int* n, const Elemental::dcomplex* alpha, Elemental::dcomplex* x, 
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* x, const int* incx,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* y, const int* incy );

void C2F(zgemv)
( const char* trans, const int* m, const int* n,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* x, const int* incx,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* y, const int* incy );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* x, const int* incx,
  const Elemental::scomplex* y, const int* incy,
        Elemental::scomplex* A, const int* lda  );

void C2F(zgerc)
( const int* m, const int* n,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* x, const int* incx,
  const Elemental::dcomplex* y, const int* incy,
        Elemental::dcomplex* A, const int* lda  );

void C2F(cgeru)
( const int* m, const int* n,
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* x, const int* incx,
  const Elemental::scomplex* y, const int* incy,
        Elemental::scomplex* A, const int* lda  );

void C2F(zgeru)
( const int* m, const int* n,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* x, const int* incx,
  const Elemental::dcomplex* y, const int* incy,
        Elemental::dcomplex* A, const int* lda  );

void C2F(chemv)
( const char* uplo, const int* m,
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* x, const int* incx,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* y, const int* incy );

void C2F(zhemv)
( const char* uplo, const int* m,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* x, const int* incx,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* y, const int* incy );

void C2F(cher)
( const char* uplo, const int* m,
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* x, const int* incx,
        Elemental::scomplex* A, const int* lda  );

void C2F(zher)
( const char* uplo, const int* m,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* x, const int* incx,
        Elemental::dcomplex* A, const int* lda  );

void C2F(cher2)
( const char* uplo, const int* m,
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* x, const int* incx,
  const Elemental::scomplex* y, const int* incy,
        Elemental::scomplex* A, const int* lda  );

void C2F(zher2)
( const char* uplo, const int* m,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* x, const int* incx,
  const Elemental::dcomplex* y, const int* incy,
        Elemental::dcomplex* A, const int* lda  );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* x, const int* incx,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* y, const int* incy );

// Recall that 'zsymv' is an auxiliary LAPACK routine
void C2F(zsymv)
( const char* uplo, const int* m,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* x, const int* incx,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* y, const int* incy );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* x, const int* incx,
        Elemental::scomplex* A, const int* lda  );

// Recall that 'zsyr' is an auxiliary LAPACK routine
void C2F(zsyr)
( const char* uplo, const int* m,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* x, const int* incx,
        Elemental::dcomplex* A, const int* lda  );
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
  const Elemental::scomplex* A, const int* lda, 
        Elemental::scomplex* x, const int* incx );

void C2F(ztrmv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const Elemental::dcomplex* A, const int* lda, 
        Elemental::dcomplex* x, const int* incx );
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
  const Elemental::scomplex* A, const int* lda, 
        Elemental::scomplex* x, const int* incx );

void C2F(ztrsv)
( const char* uplo, const char* trans, const char* diag, const int* m,
  const Elemental::dcomplex* A, const int* lda, 
        Elemental::dcomplex* x, const int* incx );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* B, const int* ldb,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* C, const int* ldc );

void C2F(zgemm)
( const char* transA, const char* transB,
  const int* m, const int* n, const int* k,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* B, const int* ldb,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* C, const int* ldc );

void C2F(chemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* B, const int* ldb,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* C, const int* ldc );

void C2F(zhemm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* B, const int* ldb,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* C, const int* ldc );

void C2F(cher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* B, const int* ldb,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* C, const int* ldc );

void C2F(zher2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* B, const int* ldb,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* C, const int* ldc );

void C2F(cherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* C, const int* ldc );

void C2F(zherk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* C, const int* ldc );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* B, const int* ldb,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* C, const int* ldc );

void C2F(zsymm)
( const char* side, const char* uplo,
  const int* m, const int* n,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* B, const int* ldb,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* C, const int* ldc );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* B, const int* ldb,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* C, const int* ldc );

void C2F(zsyr2k)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* B, const int* ldb,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* C, const int* ldc );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
  const Elemental::scomplex* beta,        
        Elemental::scomplex* C, const int* ldc );

void C2F(zsyrk)
( const char* uplo, const char* trans,
  const int* n, const int* k,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
  const Elemental::dcomplex* beta,        
        Elemental::dcomplex* C, const int* ldc );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
        Elemental::scomplex* B, const int* ldb );

void C2F(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
        Elemental::dcomplex* B, const int* ldb );
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
  const Elemental::scomplex* alpha, 
  const Elemental::scomplex* A, const int* lda,
        Elemental::scomplex* B, const int* ldb );

void C2F(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const int* m, const int* n,
  const Elemental::dcomplex* alpha, 
  const Elemental::dcomplex* A, const int* lda,
        Elemental::dcomplex* B, const int* ldb );
#endif
} // extern "C"

//----------------------------------------------------------------------------//
// Implementations begin here                                                 //
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
// Level 1 BLAS                                                               //
//----------------------------------------------------------------------------//
inline void
Elemental::wrappers::BLAS::Axpy
( const int n, const float alpha, const float* x, const int incx,
                                        float* y, const int incy )
{ C2F(saxpy)( &n, &alpha, x, &incx, y, &incy ); }

inline void
Elemental::wrappers::BLAS::Axpy
( const int n, const double alpha, const double* x, const int incx,
                                         double* y, const int incy )
{ C2F(daxpy)( &n, &alpha, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Axpy
( const int n, const scomplex alpha, const scomplex* x, const int incx,
                                           scomplex* y, const int incy )
{ C2F(caxpy)( &n, &alpha, x, &incx, y, &incy ); }

inline void
Elemental::wrappers::BLAS::Axpy
( const int n, const dcomplex alpha, const dcomplex* x, const int incx,
                                           dcomplex* y, const int incy )
{ C2F(zaxpy)( &n, &alpha, x, &incx, y, &incy ); }
#endif

inline float
Elemental::wrappers::BLAS::Dot
( const int n, const float* x, const int incx,
               const float* y, const int incy )
{ return C2F(sdot)( &n, x, &incx, y, &incy ); }

inline double
Elemental::wrappers::BLAS::Dot
( const int n, const double* x, const int incx,
               const double* y, const int incy )
{ return C2F(ddot)( &n, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline Elemental::scomplex
Elemental::wrappers::BLAS::Dot
( const int n, const Elemental::scomplex* x, const int incx,
               const Elemental::scomplex* y, const int incy )
{ return C2F(cdotc)( &n, x, &incx, y, &incy ); }

inline Elemental::dcomplex
Elemental::wrappers::BLAS::Dot
( const int n, const Elemental::dcomplex* x, const int incx,
               const Elemental::dcomplex* y, const int incy )
{ return C2F(zdotc)( &n, x, &incx, y, &incy ); }
#endif

inline float
Elemental::wrappers::BLAS::Dotc
( const int n, const float* x, const int incx,
               const float* y, const int incy )
{ return C2F(sdot)( &n, x, &incx, y, &incy ); }

inline double
Elemental::wrappers::BLAS::Dotc
( const int n, const double* x, const int incx,
               const double* y, const int incy )
{ return C2F(ddot)( &n, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline Elemental::scomplex
Elemental::wrappers::BLAS::Dotc
( const int n, const Elemental::scomplex* x, const int incx,
               const Elemental::scomplex* y, const int incy )
{ return C2F(cdotc)( &n, x, &incx, y, &incy ); }

inline Elemental::dcomplex
Elemental::wrappers::BLAS::Dotc
( const int n, const Elemental::dcomplex* x, const int incx,
               const Elemental::dcomplex* y, const int incy )
{ return C2F(zdotc)( &n, x, &incx, y, &incy ); }
#endif

inline float
Elemental::wrappers::BLAS::Dotu
( const int n, const float* x, const int incx,
               const float* y, const int incy )
{ return C2F(sdot)( &n, x, &incx, y, &incy ); }

inline double
Elemental::wrappers::BLAS::Dotu
( const int n, const double* x, const int incx,
               const double* y, const int incy )
{ return C2F(ddot)( &n, x, &incx, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline Elemental::scomplex
Elemental::wrappers::BLAS::Dotu
( const int n, const Elemental::scomplex* x, const int incx,
               const Elemental::scomplex* y, const int incy )
{ return C2F(cdotu)( &n, x, &incx, y, &incy ); }

inline Elemental::dcomplex
Elemental::wrappers::BLAS::Dotu
( const int n, const Elemental::dcomplex* x, const int incx,
               const Elemental::dcomplex* y, const int incy )
{ return C2F(zdotu)( &n, x, &incx, y, &incy ); }
#endif

inline float
Elemental::wrappers::BLAS::Nrm2
( const int n, const float* x, const int incx )
{ return C2F(snrm2)( &n, x, &incx ); }

inline double
Elemental::wrappers::BLAS::Nrm2
( const int n, const double* x, const int incx )
{ return C2F(dnrm2)( &n, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline float
Elemental::wrappers::BLAS::Nrm2
( const int n, const scomplex* x, const int incx )
{ return C2F(scnrm2)( &n, x, &incx ); }

inline double
Elemental::wrappers::BLAS::Nrm2
( const int n, const dcomplex* x, const int incx )
{ return C2F(dznrm2)( &n, x, &incx ); }
#endif

inline void
Elemental::wrappers::BLAS::Scal
( const int n, const float alpha, float* x, const int incx )
{ C2F(sscal)( &n, &alpha, x, &incx ); }

inline void
Elemental::wrappers::BLAS::Scal
( const int n, const double alpha, double* x, const int incx )
{ C2F(dscal)( &n, &alpha, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Scal
( const int n, const scomplex alpha, scomplex* x, const int incx )
{ C2F(cscal)( &n, &alpha, x, &incx ); }

inline void
Elemental::wrappers::BLAS::Scal
( const int n, const dcomplex alpha, dcomplex* x, const int incx )
{ C2F(zscal)( &n, &alpha, x, &incx ); }
#endif

//----------------------------------------------------------------------------//
// Level 2 BLAS                                                               //
//----------------------------------------------------------------------------//
inline void
Elemental::wrappers::BLAS::Gemv
( const char trans, const int m, const int n,
  const float alpha, const float* A, const int lda,
                     const float* x, const int incx,
  const float beta,        float* y, const int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

inline void
Elemental::wrappers::BLAS::Gemv
( const char trans, const int m, const int n,
  const double alpha, const double* A, const int lda,
                      const double* x, const int incx,
  const double beta,        double* y, const int incy )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Gemv
( const char trans, const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* x, const int incx,
  const scomplex beta,        scomplex* y, const int incy )
{ C2F(cgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
Elemental::wrappers::BLAS::Gemv
( const char trans, const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* x, const int incx,
  const dcomplex beta,        dcomplex* y, const int incy )
{ C2F(zgemv)( &trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }
#endif

inline void
Elemental::wrappers::BLAS::Ger
( const int m, const int n,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  )
{ C2F(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Ger
( const int m, const int n,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  )
{ C2F(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Ger
( const int m, const int n,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  )
{ C2F(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Ger
( const int m, const int n,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  )
{ C2F(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
Elemental::wrappers::BLAS::Gerc
( const int m, const int n,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  )
{ C2F(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Gerc
( const int m, const int n,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  )
{ C2F(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Gerc
( const int m, const int n,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  )
{ C2F(cgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Gerc
( const int m, const int n,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  )
{ C2F(zgerc)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
Elemental::wrappers::BLAS::Geru
( const int m, const int n,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  )
{ C2F(sger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Geru
( const int m, const int n,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  )
{ C2F(dger)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Geru
( const int m, const int n,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  )
{ C2F(cgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Geru
( const int m, const int n,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  )
{ C2F(zgeru)( &m, &n, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
Elemental::wrappers::BLAS::Hemv
( const char uplo, const int m,
  const float alpha, const float* A, const int lda,
                     const float* x, const int incx,
  const float beta,        float* y, const int incy )
{ C2F(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
Elemental::wrappers::BLAS::Hemv
( const char uplo, const int m,
  const double alpha, const double* A, const int lda,
                      const double* x, const int incx,
  const double beta,        double* y, const int incy )
{ C2F(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Hemv
( const char uplo, const int m,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* x, const int incx,
  const scomplex beta,        scomplex* y, const int incy )
{ C2F(chemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
Elemental::wrappers::BLAS::Hemv
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* x, const int incx,
  const dcomplex beta,        dcomplex* y, const int incy )
{ C2F(zhemv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }
#endif

inline void
Elemental::wrappers::BLAS::Her
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                           float* A, const int lda  )
{ C2F(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Her
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                            double* A, const int lda  )
{ C2F(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Her
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                              scomplex* A, const int lda  )
{ C2F(cher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Her
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                              dcomplex* A, const int lda  )
{ C2F(zher)( &uplo, &m, &alpha, x, &incx, A, &lda ); }
#endif

inline void
Elemental::wrappers::BLAS::Her2
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  )
{ C2F(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Her2
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  )
{ C2F(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Her2
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  )
{ C2F(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Her2
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  )
{ C2F(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }
#endif

inline void
Elemental::wrappers::BLAS::Symv
( const char uplo, const int m,
  const float alpha, const float* A, const int lda,
                     const float* x, const int incx,
  const float beta,        float* y, const int incy )
{ C2F(ssymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

inline void
Elemental::wrappers::BLAS::Symv
( const char uplo, const int m,
  const double alpha, const double* A, const int lda,
                      const double* x, const int incx,
  const double beta,        double* y, const int incy )
{ C2F(dsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Symv
( const char uplo, const int m,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* x, const int incx,
  const scomplex beta,        scomplex* y, const int incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    C2F(csymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}

inline void
Elemental::wrappers::BLAS::Symv
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* x, const int incx,
  const dcomplex beta,        dcomplex* y, const int incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    C2F(zsymv)( &uplo, &m, &alpha, A, &lda, x, &incx, &beta, y, &incy );
}
#endif

inline void
Elemental::wrappers::BLAS::Syr
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                           float* A, const int lda  )
{ C2F(ssyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Syr
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                            double* A, const int lda  )
{ C2F(dsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Syr
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                              scomplex* A, const int lda  )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    C2F(csyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}

inline void
Elemental::wrappers::BLAS::Syr
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                              dcomplex* A, const int lda  )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    C2F(zsyr)( &uplo, &m, &alpha, x, &incx, A, &lda ); 
}
#endif

inline void
Elemental::wrappers::BLAS::Syr2
( const char uplo, const int m,
  const float alpha, const float* x, const int incx,
                     const float* y, const int incy,
                           float* A, const int lda  )
{ C2F(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

inline void
Elemental::wrappers::BLAS::Syr2
( const char uplo, const int m,
  const double alpha, const double* x, const int incx,
                      const double* y, const int incy,
                            double* A, const int lda  )
{ C2F(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &lda ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Syr2
( const char uplo, const int m,
  const scomplex alpha, const scomplex* x, const int incx,
                        const scomplex* y, const int incy,
                              scomplex* A, const int lda  )
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
Elemental::wrappers::BLAS::Syr2
( const char uplo, const int m,
  const dcomplex alpha, const dcomplex* x, const int incx,
                        const dcomplex* y, const int incy,
                              dcomplex* A, const int lda  )
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
Elemental::wrappers::BLAS::Trmv
( const char uplo, const char trans, const char diag, const int m,
  const float* A, const int lda, float* x, const int incx         )
{ C2F(strmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
Elemental::wrappers::BLAS::Trmv
( const char uplo, const char trans, const char diag, const int m,
  const double* A, const int lda, double* x, const int incx       )
{ C2F(dtrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Trmv
( const char uplo, const char trans, const char diag, const int m,
  const scomplex* A, const int lda, scomplex* x, const int incx   )
{ C2F(ctrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
Elemental::wrappers::BLAS::Trmv
( const char uplo, const char trans, const char diag, const int m,
  const dcomplex* A, const int lda, dcomplex* x, const int incx   )
{ C2F(ztrmv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }
#endif

inline void
Elemental::wrappers::BLAS::Trsv
( const char uplo, const char trans, const char diag, const int m,
  const float* A, const int lda, float* x, const int incx         )
{ C2F(strsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
Elemental::wrappers::BLAS::Trsv
( const char uplo, const char trans, const char diag, const int m,
  const double* A, const int lda, double* x, const int incx       )
{ C2F(dtrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Trsv
( const char uplo, const char trans, const char diag, const int m,
  const scomplex* A, const int lda, scomplex* x, const int incx   )
{ C2F(ctrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }

inline void
Elemental::wrappers::BLAS::Trsv
( const char uplo, const char trans, const char diag, const int m,
  const dcomplex* A, const int lda, dcomplex* x, const int incx   )
{ C2F(ztrsv)( &uplo, &trans, &diag, &m, A, &lda, x, &incx ); }
#endif

//----------------------------------------------------------------------------//
// Level 3 BLAS                                                               //
//----------------------------------------------------------------------------//
inline void
Elemental::wrappers::BLAS::Gemm
( const char transA, const char transB,
  const int m, const int n, const int k, 
  const float alpha, const float* A, const int lda, 
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    C2F(sgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Gemm
( const char transA, const char transB,
  const int m, const int n, const int k, 
  const double alpha, const double* A, const int lda, 
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc )
{
    const char fixedTransA = ( transA == 'C' ? 'T' : transA );
    const char fixedTransB = ( transB == 'C' ? 'T' : transB );
    C2F(dgemm)( &fixedTransA, &fixedTransB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Gemm
( const char transA, const char transB,
  const int m, const int n, const int k, 
  const scomplex alpha, const scomplex* A, const int lda, 
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc )
{
    C2F(cgemm)( &transA, &transB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Gemm
( const char transA, const char transB,
  const int m, const int n, const int k, 
  const dcomplex alpha, const dcomplex* A, const int lda, 
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc )
{
    C2F(zgemm)( &transA, &transB, &m, &n, &k,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
Elemental::wrappers::BLAS::Hemm
( const char side, const char uplo,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc )
{
    C2F(ssymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Hemm
( const char side, const char uplo,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc )
{
    C2F(dsymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Hemm
( const char side, const char uplo,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc )
{
    C2F(chemm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Hemm
( const char side, const char uplo,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc )
{
    C2F(zhemm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
Elemental::wrappers::BLAS::Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(ssyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(dsyr2k)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc )
{
    C2F(cher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Her2k
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc )
{
    C2F(zher2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
Elemental::wrappers::BLAS::Herk
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
  const float beta,        float* C, const int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(ssyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Herk
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
  const double beta,        double* C, const int ldc )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    C2F(dsyrk)( &uplo, &transFixed, &n, &k, &alpha, A, &lda, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Herk
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
  const scomplex beta,        scomplex* C, const int ldc )
{ C2F(cherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

inline void
Elemental::wrappers::BLAS::Herk
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
  const dcomplex beta,        dcomplex* C, const int ldc )
{ C2F(zherk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }
#endif

inline void
Elemental::wrappers::BLAS::Symm
( const char side, const char uplo,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc )
{
    C2F(ssymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Symm
( const char side, const char uplo,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc )
{
    C2F(dsymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Symm
( const char side, const char uplo,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc )
{
    C2F(csymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Symm
( const char side, const char uplo,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc )
{
    C2F(zsymm)( &side, &uplo, &m, &n,
                &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
Elemental::wrappers::BLAS::Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
                     const float* B, const int ldb,
  const float beta,        float* C, const int ldc )
{
    C2F(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
                      const double* B, const int ldb,
  const double beta,        double* C, const int ldc )
{
    C2F(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
                        const scomplex* B, const int ldb,
  const scomplex beta,        scomplex* C, const int ldc )
{
    C2F(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}

inline void
Elemental::wrappers::BLAS::Syr2k
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
                        const dcomplex* B, const int ldb,
  const dcomplex beta,        dcomplex* C, const int ldc )
{
    C2F(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
}
#endif

inline void
Elemental::wrappers::BLAS::Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const float alpha, const float* A, const int lda,
  const float beta,        float* C, const int ldc )
{ C2F(ssyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

inline void
Elemental::wrappers::BLAS::Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const double alpha, const double* A, const int lda,
  const double beta,        double* C, const int ldc )
{ C2F(dsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const scomplex alpha, const scomplex* A, const int lda,
  const scomplex beta,        scomplex* C, const int ldc )
{ C2F(csyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }

inline void
Elemental::wrappers::BLAS::Syrk
( const char uplo, const char trans,
  const int n, const int k,
  const dcomplex alpha, const dcomplex* A, const int lda,
  const dcomplex beta,        dcomplex* C, const int ldc )
{ C2F(zsyrk)( &uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc ); }
#endif

inline void
Elemental::wrappers::BLAS::Trmm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                           float* B, const int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    C2F(strmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}

inline void
Elemental::wrappers::BLAS::Trmm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                            double* B, const int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );    
    C2F(dtrmm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Trmm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                              scomplex* B, const int ldb )
{
    C2F(ctrmm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}

inline void
Elemental::wrappers::BLAS::Trmm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                              dcomplex* B, const int ldb )
{
    C2F(ztrmm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
}
#endif

inline void
Elemental::wrappers::BLAS::Trsm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const float alpha, const float* A, const int lda,
                           float* B, const int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(strsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 

inline void
Elemental::wrappers::BLAS::Trsm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const double alpha, const double* A, const int lda,
                            double* B, const int ldb )
{
    const char fixedTrans = ( trans == 'C' ? 'T' : trans );
    C2F(dtrsm)( &side, &uplo, &fixedTrans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::BLAS::Trsm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const scomplex alpha, const scomplex* A, const int lda,
                              scomplex* B, const int ldb )
{
    C2F(ctrsm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 

inline void
Elemental::wrappers::BLAS::Trsm
( const char side, const char uplo, 
  const char trans, const char unit,
  const int m, const int n,
  const dcomplex alpha, const dcomplex* A, const int lda,
                              dcomplex* B, const int ldb )
{
    C2F(ztrsm)( &side, &uplo, &trans, &unit, &m, &n,
                &alpha, A, &lda, B, &ldb );
} 
#endif

#endif /* ELEMENTAL_WRAPPERS_BLAS_HPP */

