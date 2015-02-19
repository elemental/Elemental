/*
   Copyright (c) 2009-2015, Jack Poulson
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

// Basic Linear Algebra Communication Subprograms
// ==============================================
int Csys2blacs_handle( MPI_Comm comm );
void Cblacs_gridinit
( int* context, const char* order, int gridHeight, int gridWidth );
void Cblacs_gridinfo
( int  context, int* gridHeight, int* gridWidth, int* gridRow, int* gridCol );
void Cfree_blacs_system_handle( int bhandle );
void Cblacs_gridexit( int context );
void Cblacs_exit( int notDone );

// Gemv
// ====
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
// ====
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
// ====
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
// ====
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
// ====
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

// Gemm
// ====
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

// TRMM
// ====
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

// TRSM
// ====
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

// Cholesky factorization
// ======================
void EL_SCALAPACK(pspotrf)
( const char* uplo, const int* n, float* A, const int* iA, const int* jA,
  const int* descA, int* info );
void EL_SCALAPACK(pdpotrf)
( const char* uplo, const int* n, double* A, const int* iA, const int* jA,
  const int* descA, int* info );
void EL_SCALAPACK(pcpotrf)
( const char* uplo, const int* n, scomplex* A, const int* iA, const int* jA,
  const int* descA, int* info );
void EL_SCALAPACK(pzpotrf)
( const char* uplo, const int* n, dcomplex* A, const int* iA, const int* jA,
  const int* descA, int* info );

// Two-sided TRSM/TRMM
// ===================
void EL_SCALAPACK(pssyngst)
( const int* typeB, const char* uplo, const int* n, 
        float* A, const int* iA, const int* jA, const int* descA,
  const float* B, const int* iB, const int* jB, const int* descB,
        float* scale, float* work, const int* workSize, int* info );
void EL_SCALAPACK(pdsyngst)
( const int* typeB, const char* uplo, const int* n, 
        double* A, const int* iA, const int* jA, const int* descA,
  const double* B, const int* iB, const int* jB, const int* descB,
        double* scale, double* work, const int* workSize, int* info );
void EL_SCALAPACK(pchengst)
( const int* typeB, const char* uplo, const int* n, 
        scomplex* A, const int* iA, const int* jA, const int* descA,
  const scomplex* B, const int* iB, const int* jB, const int* descB,
        float* scale, scomplex* work, const int* workSize, int* info );
void EL_SCALAPACK(pzhengst)
( const int* typeB, const char* uplo, const int* n, 
        dcomplex* A, const int* iA, const int* jA, const int* descA,
  const dcomplex* B, const int* iB, const int* jB, const int* descB,
        double* scale, dcomplex* work, const int* workSize, int* info );

// Hessenberg QR algorithm
// =======================

// Aggressive Early Deflation
// --------------------------
// NOTE: ScaLAPACK currently only supports AED for real matrices
void EL_SCALAPACK(pshseqr)
( const char* job, const char* compz, 
  const int* n, const int* ilo, const int* ihi, 
  float* H, const int* descH, float* wr, float* wi, 
  float* Q, const int* descQ, float* work, const int* workSize, 
  int* iWork, const int* iWorkSize, int* info );
void EL_SCALAPACK(pdhseqr)
( const char* job, const char* compz, 
  const int* n, const int* ilo, const int* ihi, 
  double* H, const int* descH, double* wr, double* wi, 
  double* Q, const int* descQ, double* work, const int* workSize, 
  int* iWork, const int* iWorkSize, int* info );

// Pipelined QR algorithm without AED
// ----------------------------------
void EL_SCALAPACK(pslahqr)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, float* H, const int* descH,
  float* wr, float* wi, const int* iloQ, const int* ihiQ, 
  float* Q, const int* descQ,
  float* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );
void EL_SCALAPACK(pdlahqr)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, double* H, const int* descH,
  double* wr, double* wi, const int* iloQ, const int* ihiQ, 
  double* Q, const int* descQ,
  double* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );
void EL_SCALAPACK(pclahqr)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, scomplex* H, const int* descH,
  scomplex* w, const int* iloQ, const int* ihiQ, scomplex* Q, const int* descQ,
  scomplex* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );
void EL_SCALAPACK(pzlahqr)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, dcomplex* H, const int* descH,
  dcomplex* w, const int* iloQ, const int* ihiQ, dcomplex* Q, const int* descQ,
  dcomplex* work, const int* workSize, int* iWork, const int* iWorkSize, 
  int* info );

// Pipelined QR algorithm with AED for big matrices
// ------------------------------------------------
void EL_SCALAPACK(pslaqr0)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, float* H, const int* descH, 
  float* wr, float* wi, const int* iloQ, const int* ihiQ, 
  float* Q, const int* descQ, float* work, const int* workSize, 
  int* iWork, const int* iWorkSize, int* info, const int* reclevel );
void EL_SCALAPACK(pdlaqr0)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, double* H, const int* descH, 
  double* wr, double* wi, const int* iloQ, const int* ihiQ,
  double* Q, const int* descQ, double* work, const int* workSize, 
  int* iWork, const int* iWorkSize, int* info, const int* reclevel );

// Pipelined QR algorithm with AED for small matrices
// --------------------------------------------------
void EL_SCALAPACK(pslaqr1)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, float* H, const int* descH, 
  float* wr, float* wi, const int* iloQ, const int* ihiQ, 
  float* Q, const int* descQ, float* work, const int* workSize, 
  int* iWork, const int* iWorkSize, int* info );
void EL_SCALAPACK(pdlaqr1)
( const EL_FORT_LOGICAL* wantt, const EL_FORT_LOGICAL* wantz, const int* n,
  const int* ilo, const int* ihi, double* H, const int* descH, 
  double* wr, double* wi, const int* iloQ, const int* ihiQ, 
  double* Q, const int* descQ, double* work, const int* workSize, 
  int* iWork, const int* iWorkSize, int* info );

} // extern "C"

namespace El {

namespace blacs {

int Handle( MPI_Comm comm )
{ return Csys2blacs_handle( comm ); }

int GridInit( int bhandle, bool colMajor, int gridHeight, int gridWidth )
{ 
    int context = bhandle;
    const char* order = ( colMajor ? "Col" : "Row" );
    Cblacs_gridinit( &context, order, gridHeight, gridWidth ); 
    return context;
}

int GridHeight( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridHeight;
}

int GridWidth( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridWidth;
}

int GridRow( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridRow;
}

int GridCol( int context )
{
    int gridHeight, gridWidth, gridRow, gridCol;
    Cblacs_gridinfo( context, &gridHeight, &gridWidth, &gridRow, &gridCol );
    return gridCol;
}

void FreeHandle( int bhandle )
{ Cfree_blacs_system_handle( bhandle ); }

void FreeGrid( int context )
{ Cblacs_gridexit( context ); }

void Exit( bool finished )
{ 
    int notDone = ( finished ? 0 : 1 );
    Cblacs_exit( notDone ); 
}

} // namespace blacs

namespace scalapack {

// Gemv
// ====
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
// ====
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
// ====
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
// ====
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
// ====
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


// Gemm
// ====
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
// ====
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
// ====
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

// Cholesky factorization
// ======================
void Cholesky( char uplo, int n, float* A, const int* descA )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::Cholesky"))
    int iA=1,jA=1,info;
    EL_SCALAPACK(pspotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pspotrf returned with info=",info);
}

void Cholesky( char uplo, int n, double* A, const int* descA )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::Cholesky"))
    int iA=1,jA=1,info;
    EL_SCALAPACK(pdpotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pdpotrf returned with info=",info);
}

void Cholesky( char uplo, int n, scomplex* A, const int* descA )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::Cholesky"))
    int iA=1,jA=1,info;
    EL_SCALAPACK(pcpotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pcpotrf returned with info=",info);
}

void Cholesky( char uplo, int n, dcomplex* A, const int* descA )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::Cholesky"))
    int iA=1,jA=1,info;
    EL_SCALAPACK(pzpotrf)( &uplo, &n, A, &iA, &jA, descA, &info );
    if( info != 0 )
        RuntimeError("pzpotrf returned with info=",info);
}

// Two-sided TRSM
// ==============
// NOTE: It is required that B have a positive diagonal
void TwoSidedTrsm
( char uplo, int n, float* A, const int* descA, 
  const float* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrsm"))
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale, dummyWork;
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<float> work( workSize );
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pssyngst exited with info=",info);
}

void TwoSidedTrsm
( char uplo, int n, double* A, const int* descA, 
  const double* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrsm"))
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale, dummyWork;
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<double> work( workSize );
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pdsyngst exited with info=",info);
}

void TwoSidedTrsm
( char uplo, int n, scomplex* A, const int* descA, 
  const scomplex* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrsm"))
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale;
    scomplex dummyWork;
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<scomplex> work( workSize );
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pchengst exited with info=",info);
}

void TwoSidedTrsm
( char uplo, int n, dcomplex* A, const int* descA, 
  const dcomplex* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrsm"))
    int typeB=1,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale;
    dcomplex dummyWork;
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<dcomplex> work( workSize );
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pzhengst exited with info=",info);
}

// Two-sided TRMM
// ==============
// NOTE: It is required that B have a positive diagonal
void TwoSidedTrmm
( char uplo, int n, float* A, const int* descA, 
  const float* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrmm"))
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale, dummyWork;
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<float> work( workSize );
    EL_SCALAPACK(pssyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pssyngst exited with info=",info);
}

void TwoSidedTrmm
( char uplo, int n, double* A, const int* descA, 
  const double* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrmm"))
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale, dummyWork;
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork;
    vector<double> work( workSize );
    EL_SCALAPACK(pdsyngst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pdsyngst exited with info=",info);
}

void TwoSidedTrmm
( char uplo, int n, scomplex* A, const int* descA, 
  const scomplex* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrmm"))
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    float scale;
    scomplex dummyWork;
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<scomplex> work( workSize );
    EL_SCALAPACK(pchengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pchengst exited with info=",info);
}

void TwoSidedTrmm
( char uplo, int n, dcomplex* A, const int* descA, 
  const dcomplex* B, const int* descB )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::TwoSidedTrmm"))
    int typeB=2,iA=1,jA=1,iB=1,jB=1,workSize=-1,info;
    double scale;
    dcomplex dummyWork;
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      &dummyWork, &workSize, &info );

    workSize = dummyWork.real();
    vector<dcomplex> work( workSize );
    EL_SCALAPACK(pzhengst)
    ( &typeB, &uplo, &n, A, &iA, &jA, descA, B, &iB, &jB, descB, &scale, 
      work.data(), &workSize, &info );
    if( info != 0 )
        RuntimeError("pzhengst exited with info=",info);
}

// Hessenberg Schur decomposition via the QR algorithm
// ===================================================
void HessenbergSchur
( int n, float* H, const int* descH, scomplex* w, bool fullTriangle, bool aed ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    const int ilo=1, ihi=n;
    vector<float> wr(n), wi(n);
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PSHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle?'S':'E'), compz='N';

        // Query the workspace sizes
        int workSize=-1, dummyIWork, iWorkSize=-1;
        float dummyWork;

        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pshseqr exited with info=",info);
    }
    else
    {
        EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                          wantz=EL_FORT_FALSE;

        // PSLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const Int context  = descH[1];
        const Int mb       = descH[4];
        const Int nb       = descH[5];
        const Int rowAlign = descH[7];
        const Int ldH      = descH[8];
        const Int ldQ      = descQ[8];

        const Int rowCut = 0;
        const Int rowRank = blacs::GridCol( context );
        const Int rowStride = blacs::GridWidth( context );

        const Int localWidth =
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pslahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, 0, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pslahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<float>(wr[j],wi[j]);
}

void HessenbergSchur
( int n, double* H, const int* descH, dcomplex* w, bool fullTriangle, bool aed )
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    const int ilo=1, ihi=n;
    vector<double> wr(n), wi(n);
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PDHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle?'S':'E'), compz='N';

        // Query the workspace sizes
        int workSize=-1, dummyIWork, iWorkSize=-1;
        double dummyWork;
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          0, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdhseqr exited with info=",info);
    }
    else
    {
        EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                          wantz=EL_FORT_FALSE;

        // PDLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const Int context  = descH[1];
        const Int mb       = descH[4];
        const Int nb       = descH[5];
        const Int rowAlign = descH[7];
        const Int ldH      = descH[8];
        const Int ldQ      = descQ[8];

        const Int rowCut = 0;
        const Int rowRank = blacs::GridCol( context );
        const Int rowStride = blacs::GridWidth( context );

        const Int localWidth = 
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pdlahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, 0, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdlahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<double>(wr[j],wi[j]);
}

void HessenbergSchur
( int n, scomplex* H, const int* descH, scomplex* w, bool fullTriangle, 
  bool aed ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                      wantz=EL_FORT_FALSE;
    if( aed )
        LogicError("AED is not supported for complex matrices");
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    scomplex dummyWork;
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<scomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pclahqr exited with info=",info);
}

void HessenbergSchur
( int n, dcomplex* H, const int* descH, dcomplex* w, bool fullTriangle,
  bool aed ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                      wantz=EL_FORT_FALSE;
    if( aed )
        LogicError("AED is not supported for complex matrices");
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    int descQ[9] = 
        { 1, descH[1], 0, 0, descH[4], descH[5], descH[6], descH[7], 1 };
    dcomplex dummyWork;
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<dcomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, 0, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pzlahqr exited with info=",info);
}

void HessenbergSchur
( int n, float* H, const int* descH, scomplex* w, float* Q, const int* descQ, 
  bool fullTriangle, bool multiplyQ, bool aed ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    const int ilo=1, ihi=n;
    vector<float> wr(n), wi(n);
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PSHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle?'S':'E'), compz=(multiplyQ?'V':'I');

        // Query the workspace sizes. Due to a bug in p{s,d}hseqr's workspace
        // querying, which is located in p{s,d}laqr1, 
        //    https://github.com/poulson/scalapack/commits/master 
        // we must be a bit more careful.
        int workSize=-1, dummyIWork=3, iWorkSize=-1;
        float dummyWork;
        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pshseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pshseqr exited with info=",info);
    }
    else
    {
        if( multiplyQ == false )
            LogicError("Forcing the matrix to identity is not yet supported");
        EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                          wantz=EL_FORT_TRUE;

        // PSLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const Int context  = descH[1];
        const Int mb       = descH[4];
        const Int nb       = descH[5];
        const Int rowAlign = descH[7];
        const Int ldH      = descH[8];
        const Int ldQ      = descQ[8];

        const Int rowCut = 0;
        const Int rowRank = blacs::GridCol( context );
        const Int rowStride = blacs::GridWidth( context );

        const Int localWidth =
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<float> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pslahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, Q, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pslahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<float>(wr[j],wi[j]);
}

void HessenbergSchur
( int n, double* H, const int* descH, dcomplex* w, 
  double* Q, const int* descQ, bool fullTriangle, bool multiplyQ, bool aed ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    const int ilo=1, ihi=n;
    vector<double> wr(n), wi(n);
    int info;
    if( aed )
    {
        cerr << 
          "WARNING: PDHSEQR seems to have a bug in its eigenvalue reordering" 
          << endl;
        const char job=(fullTriangle?'S':'E'), compz=(multiplyQ?'V':'I');

        // Query the workspace sizes. Due to a bug in p{s,d}hseqr's workspace
        // querying, which is located in p{s,d}laqr1, 
        //    https://github.com/poulson/scalapack/commits/master 
        // we must be a bit more careful.
        int workSize=-1, dummyIWork=3, iWorkSize=-1;
        double dummyWork;
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

        // Compute the eigenvalues in parallel
        workSize = dummyWork;
        iWorkSize = dummyIWork;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);
        EL_SCALAPACK(pdhseqr)
        ( &job, &compz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          Q, descQ, work.data(), &workSize, iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdhseqr exited with info=",info);
    }
    else
    {
        if( multiplyQ == false )
            LogicError("Forcing the matrix to identity is not yet supported");
        EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                          wantz=EL_FORT_TRUE;

        // PDLAHQR does not support a workspace query and instead assumes
        // that the workspace is at least 
        //     3*N + max(2*max(ldH,ldQ)+2*localWidth,numRowBlocks)

        const Int context  = descH[1];
        const Int mb       = descH[4];
        const Int nb       = descH[5];
        const Int rowAlign = descH[7];
        const Int ldH      = descH[8];
        const Int ldQ      = descQ[8];

        const Int rowCut = 0;
        const Int rowRank = blacs::GridCol( context );
        const Int rowStride = blacs::GridWidth( context );

        const Int localWidth =
            BlockedLength(n,rowRank,rowAlign,nb,rowCut,rowStride);
        const Int numRowBlocks = n/mb + 1;
        int workSize= 3*n + Max(2*Max(ldH,ldQ)+2*localWidth,numRowBlocks),
            iWorkSize=0;
        vector<double> work(workSize);
        vector<int> iWork(iWorkSize);

        // Compute the eigenvalues in parallel
        EL_SCALAPACK(pdlahqr)
        ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, wr.data(), wi.data(), 
          &ilo, &ihi, Q, descQ, work.data(), &workSize, 
          iWork.data(), &iWorkSize, &info );
        if( info != 0 )
            RuntimeError("pdlahqr exited with info=",info);
    }
    // Combine the real and imaginary components of the eigenvalues
    for( int j=0; j<n; ++j )
        w[j] = std::complex<double>(wr[j],wi[j]);
}

void HessenbergSchur
( int n, scomplex* H, const int* descH, scomplex* w, 
  scomplex* Q, const int* descQ, bool fullTriangle, bool multiplyQ, bool aed ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    if( !multiplyQ )
        LogicError("Forcing the matrix to identity is not yet supported");
    if( aed )
        LogicError("AED is not supported for complex matrices");
    EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                      wantz=EL_FORT_TRUE;
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    scomplex dummyWork;
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<scomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pclahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pclahqr exited with info=",info);
}

void HessenbergSchur
( int n, dcomplex* H, const int* descH, dcomplex* w, 
  dcomplex* Q, const int* descQ, bool fullTriangle, bool multiplyQ, bool aed ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergSchur"))
    if( !multiplyQ )
        LogicError("Forcing the matrix to identity is not yet supported");
    if( aed )
        LogicError("AED is not supported for complex matrices");
    EL_FORT_LOGICAL wantt=(fullTriangle?EL_FORT_TRUE:EL_FORT_FALSE), 
                      wantz=EL_FORT_TRUE;
    const int ilo=1, ihi=n;

    // Query the workspace sizes
    // NOTE: dummyIWork will currently be left unmodified (hence the
    //       zero initialization)!
    int workSize=-1, dummyIWork=0, iWorkSize=-1, info;
    dcomplex dummyWork;
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      &dummyWork, &workSize, &dummyIWork, &iWorkSize, &info );

    // Compute the eigenvalues in parallel
    workSize = dummyWork.real();
    iWorkSize = dummyIWork;
    vector<dcomplex> work(workSize);
    vector<int> iWork(iWorkSize);
    EL_SCALAPACK(pzlahqr)
    ( &wantt, &wantz, &n, &ilo, &ihi, H, descH, w, &ilo, &ihi, Q, descQ,
      work.data(), &workSize, iWork.data(), &iWorkSize, &info );
    if( info != 0 )
        RuntimeError("pzlahqr exited with info=",info);
}

// Hessenberg eigenvalues/pairs via the QR algorithm
// =================================================
void HessenbergEig( int n, float* H, const int* descH, scomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergEig"))
    HessenbergSchur( n, H, descH, w, false );
}

void HessenbergEig( int n, double* H, const int* descH, dcomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergEig"))
    HessenbergSchur( n, H, descH, w, false );
}

void HessenbergEig( int n, scomplex* H, const int* descH, scomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergEig"))
    HessenbergSchur( n, H, descH, w, false );
}

void HessenbergEig( int n, dcomplex* H, const int* descH, dcomplex* w ) 
{
    DEBUG_ONLY(CallStackEntry cse("scalapack::HessenbergEig"))
    HessenbergSchur( n, H, descH, w, false );
}

} // namespace scalapack
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK
