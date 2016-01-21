/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_SCALAPACK_PBLAS_HPP
#define EL_IMPORTS_SCALAPACK_PBLAS_HPP

#ifdef EL_HAVE_SCALAPACK

namespace El {
namespace pblas {

// Level 2
// =======

// Gemv
// ----
void Gemv
( char trans, int m, int n, 
  float alpha, const float* A, const int* descA, 
               const float* x, const int* descx, int incx, 
  float beta,        float* y, const int* descy, int incy );
void Gemv
( char trans, int m, int n, 
  double alpha, const double* A, const int* descA, 
                const double* x, const int* descx, int incx, 
  double beta,        double* y, const int* descy, int incy );
void Gemv
( char trans, int m, int n, 
  scomplex alpha, const scomplex* A, const int* descA, 
                  const scomplex* x, const int* descx, int incx, 
  scomplex beta,        scomplex* y, const int* descy, int incy );
void Gemv
( char trans, int m, int n, 
  dcomplex alpha, const dcomplex* A, const int* descA, 
                  const dcomplex* x, const int* descx, int incx, 
  dcomplex beta,        dcomplex* y, const int* descy, int incy );

// Hemv
// ----
void Hemv
( char uplo, int n,
  scomplex alpha, const scomplex* A, const int* descA, 
                  const scomplex* x, const int* descx, int incx,
  scomplex beta,        scomplex* y, const int* descy, int incy );
void Hemv
( char uplo, int n,
  dcomplex alpha, const dcomplex* A, const int* descA, 
                  const dcomplex* x, const int* descx, int incx,
  dcomplex beta,        dcomplex* y, const int* descy, int incy );

// Symv
// ----
void Symv
( char uplo, int n,
  float alpha, const float* A, const int* descA,
               const float* x, const int* descx, int incx,
  float beta,        float* y, const int* descy, int incy );
void Symv
( char uplo, int n,
  double alpha, const double* A, const int* descA,
                const double* x, const int* descx, int incx,
  double beta,        double* y, const int* descy, int incy );

// Trmv
// ----
void Trmv
( char uplo, char trans, char diag, int n,
  const float* A, const int* descA, 
        float* x, const int* descx, int incx );
void Trmv
( char uplo, char trans, char diag, int n,
  const double* A, const int* descA, 
        double* x, const int* descx, int incx );
void Trmv
( char uplo, char trans, char diag, int n,
  const scomplex* A, const int* descA, 
        scomplex* x, const int* descx, int incx );
void Trmv
( char uplo, char trans, char diag, int n,
  const dcomplex* A, const int* descA, 
        dcomplex* x, const int* descx, int incx );

// Trsv
// ----
void Trsv
( char uplo, char trans, char diag, int n,
  const float* A, const int* descA, 
        float* x, const int* descx, int incx );
void Trsv
( char uplo, char trans, char diag, int n,
  const double* A, const int* descA, 
        double* x, const int* descx, int incx );
void Trsv
( char uplo, char trans, char diag, int n,
  const scomplex* A, const int* descA, 
        scomplex* x, const int* descx, int incx );
void Trsv
( char uplo, char trans, char diag, int n,
  const dcomplex* A, const int* descA, 
        dcomplex* x, const int* descx, int incx );

// Level 3
// =======

// Gemm
// ----
void Gemm
( char transa, char transb, int m, int n, int k,
  float alpha, const float* A, const int* descA,
               const float* B, const int* descB,
  float beta,        float* C, const int* descC );
void Gemm
( char transa, char transb, int m, int n, int k,
  double alpha, const double* A, const int* descA,
                const double* B, const int* descB,
  double beta,        double* C, const int* descC );
void Gemm
( char transa, char transb, int m, int n, int k,
  scomplex alpha, const scomplex* A, const int* descA,
                  const scomplex* B, const int* descB,
  scomplex beta,        scomplex* C, const int* descC );
void Gemm
( char transa, char transb, int m, int n, int k,
  dcomplex alpha, const dcomplex* A, const int* descA,
                  const dcomplex* B, const int* descB,
  dcomplex beta,        dcomplex* C, const int* descC );

// Trmm
// ----
void Trmm
( char side, char uplo, char trans, char diag, int m, int n,
  float alpha, const float* A, const int* descA,
                     float* B, const int* descB );
void Trmm
( char side, char uplo, char trans, char diag, int m, int n,
  double alpha, const double* A, const int* descA,
                      double* B, const int* descB );
void Trmm
( char side, char uplo, char trans, char diag, int m, int n,
  scomplex alpha, const scomplex* A, const int* descA,
                        scomplex* B, const int* descB );
void Trmm
( char side, char uplo, char trans, char diag, int m, int n,
  dcomplex alpha, const dcomplex* A, const int* descA,
                        dcomplex* B, const int* descB );

// Trsm
// ----
void Trsm
( char side, char uplo, char trans, char diag, int m, int n,
  float alpha, const float* A, const int* descA,
                     float* B, const int* descB );
void Trsm
( char side, char uplo, char trans, char diag, int m, int n,
  double alpha, const double* A, const int* descA,
                      double* B, const int* descB );
void Trsm
( char side, char uplo, char trans, char diag, int m, int n,
  scomplex alpha, const scomplex* A, const int* descA,
                        scomplex* B, const int* descB );
void Trsm
( char side, char uplo, char trans, char diag, int m, int n,
  dcomplex alpha, const dcomplex* A, const int* descA,
                        dcomplex* B, const int* descB );

} // namespace pblas
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK

#endif // ifndef EL_IMPORTS_SCALAPACK_PBLAS_HPP
