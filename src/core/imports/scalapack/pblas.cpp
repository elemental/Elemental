/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#ifdef EL_HAVE_SCALAPACK

using El::scomplex;
using El::dcomplex;

extern "C" {

// Level 1
// =======

// Level 2
// =======

// Gemv
// ----
void EL_SCALAPACK(psgemv)
( const char* trans, 
  const int* m, const int* n,
  const float* alpha, 
  const float* A, const int* iA, const int* jA, const int* descA,
  const float* x, const int* ix, const int* jx, const int* descx, 
  const int* incx, 
  const float* beta,
  const float* y, const int* iy, const int* jy, const int* descy,
  const int* incy );
void EL_SCALAPACK(pdgemv)
( const char* trans, 
  const int* m, const int* n,
  const double* alpha, 
  const double* A, const int* iA, const int* jA, const int* descA,
  const double* x, const int* ix, const int* jx, const int* descx, 
  const int* incx, 
  const double* beta,
  const double* y, const int* iy, const int* jy, const int* descy,
  const int* incy );
void EL_SCALAPACK(pcgemv)
( const char* trans, 
  const int* m, const int* n,
  const scomplex* alpha, 
  const scomplex* A, const int* iA, const int* jA, const int* descA,
  const scomplex* x, const int* ix, const int* jx, const int* descx, 
  const int* incx, 
  const scomplex* beta,
  const scomplex* y, const int* iy, const int* jy, const int* descy,
  const int* incy );
void EL_SCALAPACK(pzgemv)
( const char* trans, 
  const int* m, const int* n,
  const dcomplex* alpha, 
  const dcomplex* A, const int* iA, const int* jA, const int* descA,
  const dcomplex* x, const int* ix, const int* jx, const int* descx, 
  const int* incx, 
  const dcomplex* beta,
  const dcomplex* y, const int* iy, const int* jy, const int* descy,
  const int* incy );

// Hemv
// ----
void EL_SCALAPACK(pchemv)
( const char* uplo, const int* n,
  const scomplex* alpha,
  const scomplex* A, const int* iA, const int* jA, const int* descA,
  const scomplex* x, const int* ix, const int* jx, const int* descx, 
  const int* incx,
  const scomplex* beta,
        scomplex* y, const int* iy, const int* jy, const int* descy,
  const int* incy );
void EL_SCALAPACK(pzhemv)
( const char* uplo, const int* n,
  const dcomplex* alpha,
  const dcomplex* A, const int* iA, const int* jA, const int* descA,
  const dcomplex* x, const int* ix, const int* jx, const int* descx, 
  const int* incx,
  const dcomplex* beta,
        dcomplex* y, const int* iy, const int* jy, const int* descy,
  const int* incy );

// Symv
// ----
void EL_SCALAPACK(pssymv)
( const char* uplo, const int* n,
  const float* alpha,
  const float* A, const int* iA, const int* jA, const int* descA,
  const float* x, const int* ix, const int* jx, const int* descx, 
  const int* incx,
  const float* beta,
        float* y, const int* iy, const int* jy, const int* descy,
  const int* incy );
void EL_SCALAPACK(pdsymv)
( const char* uplo, const int* n,
  const double* alpha,
  const double* A, const int* iA, const int* jA, const int* descA,
  const double* x, const int* ix, const int* jx, const int* descx, 
  const int* incx,
  const double* beta,
        double* y, const int* iy, const int* jy, const int* descy,
  const int* incy );

// Trmv
// ----
void EL_SCALAPACK(pstrmv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const float* A, const int* iA, const int* jA, const int* descA,
        float* x, const int* ix, const int* jx, const int* descx,
  const int* incx );
void EL_SCALAPACK(pdtrmv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const double* A, const int* iA, const int* jA, const int* descA,
        double* x, const int* ix, const int* jx, const int* descx,
  const int* incx );
void EL_SCALAPACK(pctrmv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const scomplex* A, const int* iA, const int* jA, const int* descA,
        scomplex* x, const int* ix, const int* jx, const int* descx,
  const int* incx );
void EL_SCALAPACK(pztrmv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const dcomplex* A, const int* iA, const int* jA, const int* descA,
        dcomplex* x, const int* ix, const int* jx, const int* descx,
  const int* incx );

// Trsv
// ----
void EL_SCALAPACK(pstrsv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const float* A, const int* iA, const int* jA, const int* descA,
        float* x, const int* ix, const int* jx, const int* descx,
  const int* incx );
void EL_SCALAPACK(pdtrsv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const double* A, const int* iA, const int* jA, const int* descA,
        double* x, const int* ix, const int* jx, const int* descx,
  const int* incx );
void EL_SCALAPACK(pctrsv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const scomplex* A, const int* iA, const int* jA, const int* descA,
        scomplex* x, const int* ix, const int* jx, const int* descx,
  const int* incx );
void EL_SCALAPACK(pztrsv)
( const char* uplo, const char* trans, const char* diag, const int* n,
  const dcomplex* A, const int* iA, const int* jA, const int* descA,
        dcomplex* x, const int* ix, const int* jx, const int* descx,
  const int* incx );

// Level 3
// =======

// Gemm
// ----
void EL_SCALAPACK(psgemm)
( const char* transa, const char* transb,
  const int* m, const int* n, const int* k,
  const float* alpha, 
  const float* A, const int* iA, const int* jA, const int* descA,
  const float* B, const int* iB, const int* jB, const int* descB,
  const float* beta,
  const float* C, const int* iC, const int* jC, const int* descC );
void EL_SCALAPACK(pdgemm)
( const char* transa, const char* transb,
  const int* m, const int* n, const int* k,
  const double* alpha, 
  const double* A, const int* iA, const int* jA, const int* descA,
  const double* B, const int* iB, const int* jB, const int* descB,
  const double* beta,
  const double* C, const int* iC, const int* jC, const int* descC );
 void EL_SCALAPACK(pcgemm)
( const char* transa, const char* transb,
  const int* m, const int* n, const int* k,
  const scomplex* alpha, 
  const scomplex* A, const int* iA, const int* jA, const int* descA,
  const scomplex* B, const int* iB, const int* jB, const int* descB,
  const scomplex* beta,
  const scomplex* C, const int* iC, const int* jC, const int* descC );
void EL_SCALAPACK(pzgemm)
( const char* transa, const char* transb,
  const int* m, const int* n, const int* k,
  const dcomplex* alpha, 
  const dcomplex* A, const int* iA, const int* jA, const int* descA,
  const dcomplex* B, const int* iB, const int* jB, const int* descB,
  const dcomplex* beta,
  const dcomplex* C, const int* iC, const int* jC, const int* descC );

// Trmm
// ----
void EL_SCALAPACK(pstrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const float* alpha, 
  const float* A, const int* iA, const int* jA, const int* descA,
        float* B, const int* iB, const int* jB, const int* descB );
void EL_SCALAPACK(pdtrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const double* alpha, 
  const double* A, const int* iA, const int* jA, const int* descA,
        double* B, const int* iB, const int* jB, const int* descB );
void EL_SCALAPACK(pctrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const scomplex* alpha, 
  const scomplex* A, const int* iA, const int* jA, const int* descA,
        scomplex* B, const int* iB, const int* jB, const int* descB );
void EL_SCALAPACK(pztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const dcomplex* alpha, 
  const dcomplex* A, const int* iA, const int* jA, const int* descA,
        dcomplex* B, const int* iB, const int* jB, const int* descB );

// Trsm
// ----
void EL_SCALAPACK(pstrsm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const float* alpha, 
  const float* A, const int* iA, const int* jA, const int* descA,
        float* B, const int* iB, const int* jB, const int* descB );
void EL_SCALAPACK(pdtrsm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const double* alpha, 
  const double* A, const int* iA, const int* jA, const int* descA,
        double* B, const int* iB, const int* jB, const int* descB );
void EL_SCALAPACK(pctrsm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const scomplex* alpha, 
  const scomplex* A, const int* iA, const int* jA, const int* descA,
        scomplex* B, const int* iB, const int* jB, const int* descB );
void EL_SCALAPACK(pztrsm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const int* m, const int* n, const dcomplex* alpha, 
  const dcomplex* A, const int* iA, const int* jA, const int* descA,
        dcomplex* B, const int* iB, const int* jB, const int* descB );

} // extern "C"

namespace El {
namespace pblas {

// Level 1
// =======

// Level 2
// =======

// Gemv
// ----
void Gemv
( char trans, int m, int n,
  float alpha, const float* A, const int* descA,
               const float* x, const int* descx, const int incx,
  float beta,        float* y, const int* descy, const int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(psgemv)
    ( &trans, &m, &n,
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

void Gemv
( char trans, int m, int n,
  double alpha, const double* A, const int* descA,
                const double* x, const int* descx, const int incx,
  double beta,        double* y, const int* descy, const int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(pdgemv)
    ( &trans, &m, &n,
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

void Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, const int* descA,
                  const scomplex* x, const int* descx, const int incx,
  scomplex beta,        scomplex* y, const int* descy, const int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(pcgemv)
    ( &trans, &m, &n,
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

void Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, const int* descA,
                  const dcomplex* x, const int* descx, const int incx,
  dcomplex beta,        dcomplex* y, const int* descy, const int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(pzgemv)
    ( &trans, &m, &n,
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

// Hemv
// ----
void Hemv
( char uplo, int n,
  scomplex alpha, const scomplex* A, const int* descA,
                  const scomplex* x, const int* descx, int incx,
  scomplex beta,        scomplex* y, const int* descy, int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(pchemv)
    ( &uplo, &n, 
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

void Hemv
( char uplo, int n,
  dcomplex alpha, const dcomplex* A, const int* descA,
                  const dcomplex* x, const int* descx, int incx,
  dcomplex beta,        dcomplex* y, const int* descy, int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(pzhemv)
    ( &uplo, &n, 
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

// Symv
// ----
void Symv
( char uplo, int n,
  float alpha, const float* A, const int* descA,
               const float* x, const int* descx, int incx,
  float beta,        float* y, const int* descy, int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(pssymv)
    ( &uplo, &n, 
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

void Symv
( char uplo, int n,
  double alpha, const double* A, const int* descA,
                const double* x, const int* descx, int incx,
  double beta,        double* y, const int* descy, int incy )
{
    int iA=1, jA=1, ix=1, jx=1, iy=1, jy=1;
    EL_SCALAPACK(pdsymv)
    ( &uplo, &n, 
      &alpha, A, &iA, &jA, descA, x, &ix, &jx, descx, &incx,
      &beta,                      y, &iy, &jy, descy, &incy );
}

// Trmv
// ----
void Trmv
( char uplo, char trans, char diag, int n,
  const float* A, const int* descA,
        float* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pstrmv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

void Trmv
( char uplo, char trans, char diag, int n,
  const double* A, const int* descA,
        double* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pdtrmv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

void Trmv
( char uplo, char trans, char diag, int n,
  const scomplex* A, const int* descA,
        scomplex* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pctrmv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

void Trmv
( char uplo, char trans, char diag, int n,
  const dcomplex* A, const int* descA,
        dcomplex* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pztrmv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

// Trsv
// ----
void Trsv
( char uplo, char trans, char diag, int n,
  const float* A, const int* descA,
        float* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pstrsv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

void Trsv
( char uplo, char trans, char diag, int n,
  const double* A, const int* descA,
        double* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pdtrsv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

void Trsv
( char uplo, char trans, char diag, int n,
  const scomplex* A, const int* descA,
        scomplex* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pctrsv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

void Trsv
( char uplo, char trans, char diag, int n,
  const dcomplex* A, const int* descA,
        dcomplex* x, const int* descx, const int incx )
{
    int iA=1, jA=1, ix=1, jx=1;
    EL_SCALAPACK(pztrsv)
    ( &uplo, &trans, &diag, &n,
      A, &iA, &jA, descA, x, &ix, &jx, descx, &incx );
}

// Level 3
// =======

// Gemm
// ----
void Gemm
( char transa, char transb, int m, int n, int k,
  float alpha, const float* A, const int* descA,
               const float* B, const int* descB,
  float beta,        float* C, const int* descC )
{
    int iA=1, jA=1, iB=1, jB=1, iC=1, jC=1;
    EL_SCALAPACK(psgemm)
    ( &transa, &transb, &m, &n, &k,
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB,
      &beta,  C, &iC, &jC, descC );
}

void Gemm
( char transa, char transb, int m, int n, int k,
  double alpha, const double* A, const int* descA,
                const double* B, const int* descB,
  double beta,        double* C, const int* descC )
{
    int iA=1, jA=1, iB=1, jB=1, iC=1, jC=1;
    EL_SCALAPACK(pdgemm)
    ( &transa, &transb, &m, &n, &k,
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB,
      &beta,  C, &iC, &jC, descC );
}

void Gemm
( char transa, char transb, int m, int n, int k,
  scomplex alpha, const scomplex* A, const int* descA,
                  const scomplex* B, const int* descB,
  scomplex beta,        scomplex* C, const int* descC )
{
    int iA=1, jA=1, iB=1, jB=1, iC=1, jC=1;
    EL_SCALAPACK(pcgemm)
    ( &transa, &transb, &m, &n, &k,
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB,
      &beta,  C, &iC, &jC, descC );
}

void Gemm
( char transa, char transb, int m, int n, int k,
  dcomplex alpha, const dcomplex* A, const int* descA,
                  const dcomplex* B, const int* descB,
  dcomplex beta,        dcomplex* C, const int* descC )
{
    int iA=1, jA=1, iB=1, jB=1, iC=1, jC=1;
    EL_SCALAPACK(pzgemm)
    ( &transa, &transb, &m, &n, &k,
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB,
      &beta,  C, &iC, &jC, descC );
}

// Trmm
// ----
void Trmm
( char side, char uplo, char trans, char diag, int m, int n, 
  float alpha, const float* A, const int* descA, 
                     float* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pstrmm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

void Trmm
( char side, char uplo, char trans, char diag, int m, int n, 
  double alpha, const double* A, const int* descA, 
                      double* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pdtrmm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

void Trmm
( char side, char uplo, char trans, char diag, int m, int n, 
  scomplex alpha, const scomplex* A, const int* descA, 
                        scomplex* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pctrmm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

void Trmm
( char side, char uplo, char trans, char diag, int m, int n, 
  dcomplex alpha, const dcomplex* A, const int* descA, 
                        dcomplex* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pztrmm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

// Trsm
// ----
void Trsm
( char side, char uplo, char trans, char diag, int m, int n, 
  float alpha, const float* A, const int* descA, 
                     float* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pstrsm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

void Trsm
( char side, char uplo, char trans, char diag, int m, int n, 
  double alpha, const double* A, const int* descA, 
                      double* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pdtrsm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

void Trsm
( char side, char uplo, char trans, char diag, int m, int n, 
  scomplex alpha, const scomplex* A, const int* descA, 
                        scomplex* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pctrsm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

void Trsm
( char side, char uplo, char trans, char diag, int m, int n, 
  dcomplex alpha, const dcomplex* A, const int* descA, 
                        dcomplex* B, const int* descB )
{ 
    int iA=1, jA=1, iB=1, jB=1;
    EL_SCALAPACK(pztrsm)
    ( &side, &uplo, &trans, &diag, &m, &n, 
      &alpha, A, &iA, &jA, descA, B, &iB, &jB, descB );
}

} // namespace pblas
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK
