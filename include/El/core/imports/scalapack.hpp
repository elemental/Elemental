/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_SCALAPACK_HPP
#define EL_IMPORTS_SCALAPACK_HPP

#ifdef EL_HAVE_SCALAPACK

namespace El {

namespace blacs {

int Handle( MPI_Comm comm );
int GridInit( int bhandle, bool colMajor, int gridHeight, int gridWidth );
int GridHeight( int context );
int GridWidth( int context );
int GridRow( int context );
int GridCol( int context );
void FreeHandle( int bhandle );
void FreeGrid( int context );
void Exit( bool finished=false );

typedef typename std::array<int,9> Desc;

} // namespace blacs

namespace scalapack {

// NOTE: The vast majority of these routines are for benchmarking purposes,
//       but the Hessenberg QR algorithm is actively used by Elemental's
//       Pseudospectrum routine.

// BLAS 2
// ======

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

// BLAS 3
// ======

// Gemm
// ----
void Gemm
( char trans, char transb, int m, int n, int k,
  float alpha, const float* A, const int* descA,
               const float* B, const int* descB,
  float beta,        float* C, const int* descC );
void Gemm
( char trans, char transb, int m, int n, int k,
  double alpha, const double* A, const int* descA,
                const double* B, const int* descB,
  double beta,        double* C, const int* descC );
void Gemm
( char trans, char transb, int m, int n, int k,
  scomplex alpha, const scomplex* A, const int* descA,
                  const scomplex* B, const int* descB,
  scomplex beta,        scomplex* C, const int* descC );
void Gemm
( char trans, char transb, int m, int n, int k,
  dcomplex alpha, const dcomplex* A, const int* descA,
                  const dcomplex* B, const int* descB,
  dcomplex beta,        dcomplex* C, const int* descC );

// Factorizations
// ==============

// Cholesky decomposition
// ----------------------
void Cholesky( char uplo, int n, float* A, const int* desca );
void Cholesky( char uplo, int n, double* A, const int* desca );
void Cholesky( char uplo, int n, scomplex* A, const int* desca );
void Cholesky( char uplo, int n, dcomplex* A, const int* desca );

// Spectral analysis
// =================

// Hermitian eigenvalue decomposition
// ----------------------------------

// Compute eigenvalues
// ^^^^^^^^^^^^^^^^^^^

// All eigenvalues
// """""""""""""""
void HermitianEig
( char uplo, int n, float* A, const int* desca, float* w, 
  float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, const int* desca, double* w, 
  double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, const int* desca, float* w, 
  float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, const int* desca, double* w, 
  double abstol=0 );

// Floating-point range
// """"""""""""""""""""
int HermitianEig
( char uplo, int n, float* A, const int* desca, float* w,
  float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, double* A, const int* desca, double* w,
  double vl, double vu, double abstol=0 );
int HermitianEig
( char uplo, int n, scomplex* A, const int* desca, float* w,
  float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, dcomplex* A, const int* desca, double* w,
  double vl, double vu, double abstol=0 );

// Index range
// """""""""""
void HermitianEig
( char uplo, int n, float* A, const int* desca, float* w,
  int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, const int* desca, double* w,
  int il, int iu, double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, const int* desca, float* w,
  int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, const int* desca, double* w,
  int il, int iu, double abstol=0 );

// Compute eigenpairs
// ^^^^^^^^^^^^^^^^^^

// All eigenpairs
// """"""""""""""
void HermitianEig
( char uplo, int n, float* A, const int* desca, float* w, 
  float* Z, const int* descz, float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, const int* desca, double* w, 
  double* Z, const int* descz, double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, const int* desca, float* w, 
  scomplex* Z, const int* descz, float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, const int* desca, double* w, 
  dcomplex* Z, const int* descz, double abstol=0 );

// Floating-point range
// """"""""""""""""""""
int HermitianEig
( char uplo, int n, float* A, const int* desca, float* w, 
  float* Z, const int* descz, float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, double* A, const int* desca, double* w, 
  double* Z, const int* descz, double vl, double vu, double abstol=0 );
int HermitianEig
( char uplo, int n, scomplex* A, const int* desca, float* w, 
  scomplex* Z, const int* descz, float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, dcomplex* A, const int* desca, double* w, 
  dcomplex* Z, const int* descz, double vl, double vu, double abstol=0 );

// Index range
// """""""""""
void HermitianEig
( char uplo, int n, float* A, const int* desca, float* w, 
  float* Z, const int* descz, int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, const int* desca, double* w, 
  double* Z, const int* descz, int il, int iu, double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, const int* desca, float* w, 
  scomplex* Z, const int* descz, int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, const int* desca, double* w, 
  dcomplex* Z, const int* descz, int il, int iu, double abstol=0 );

// Hessenberg Schur decomposition via the QR algorithm
// ---------------------------------------------------
// NOTE: In all of these routines, the matrix needs to be explicitly 
//       upper-Hessenberg before the call, otherwise behavior is unpredictable

void HessenbergSchur
( int n, float* H, const int* desch, scomplex* w, bool fullTriangle=false, 
  bool aed=false );
void HessenbergSchur
( int n, double* H, const int* desch, dcomplex* w, bool fullTriangle=false,
  bool aed=false );
void HessenbergSchur
( int n, scomplex* H, const int* desch, scomplex* w, bool fullTriangle=false,
  bool aed=false );
void HessenbergSchur
( int n, dcomplex* H, const int* desch, dcomplex* w, bool fullTriangle=false,
  bool aed=false );

void HessenbergSchur
( int n, float* H, const int* desch, scomplex* w, 
  float* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false,
  bool aed=false );
void HessenbergSchur
( int n, double* H, const int* desch, dcomplex* w, 
  double* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false,
  bool aed=false );
void HessenbergSchur
( int n, scomplex* H, const int* desch, scomplex* w, 
  scomplex* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false,
  bool aed=false );
void HessenbergSchur
( int n, dcomplex* H, const int* desch, dcomplex* w, 
  dcomplex* Q, const int* descq, bool fullTriangle=true, bool multiplyQ=false,
  bool aed=false );

// Hessenberg eigenvalues/pairs
// ----------------------------
void HessenbergEig( int n, float* H, const int* desch, scomplex* w );
void HessenbergEig( int n, double* H, const int* desch, dcomplex* w );
void HessenbergEig( int n, scomplex* H, const int* desch, scomplex* w );
void HessenbergEig( int n, dcomplex* H, const int* desch, dcomplex* w );

// TODO: Compute the eigenvectors

} // namespace scalapack
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK

#endif // ifndef EL_IMPORTS_SCALAPACK_HPP
