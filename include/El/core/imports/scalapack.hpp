/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_SCALAPACK_HPP
#define EL_IMPORTS_SCALAPACK_HPP

namespace El {

inline void AssertScaLAPACKSupport()
{
#ifndef EL_HAVE_SCALAPACK
    LogicError("Elemental was not compiled with ScaLAPACK support");
#endif
}

} // namespace El

#ifdef EL_HAVE_SCALAPACK

#if defined(EL_BUILT_SCALAPACK)
# define EL_SCALAPACK(name) FC_GLOBAL(name,name)
#else
# if defined(EL_HAVE_SCALAPACK_SUFFIX)
#  define EL_SCALAPACK(name) EL_CONCAT(name,EL_SCALAPACK_SUFFIX)
# else
#  define EL_SCALAPACK(name) name
# endif
#endif

// TODO: Decide which routines should be modified to use 64-bit integers if
//       BLAS/LAPACK were modified to do so...
#include "./scalapack/blacs.hpp"
#include "./scalapack/pblas.hpp"

namespace El {
namespace scalapack {

// Factorizations
// ==============

// Cholesky
// --------
void Cholesky( char uplo, int n, float* A, const int* descA );
void Cholesky( char uplo, int n, double* A, const int* descA );
void Cholesky( char uplo, int n, scomplex* A, const int* descA );
void Cholesky( char uplo, int n, dcomplex* A, const int* descA );

// QR
// --
void QR( int m, int n, float* A, const int* descA, float* tau );
void QR( int m, int n, double* A, const int* descA, double* tau );
void QR( int m, int n, scomplex* A, const int* descA, scomplex* tau );
void QR( int m, int n, dcomplex* A, const int* descA, dcomplex* tau );

// Solvers
// =======

// General linear solver
// ---------------------
void LinearSolve
( int n, int numRhs, float* A, const int* descA,
  int* ipiv, float* B, const int* descb );
void LinearSolve
( int n, int numRhs, double* A, const int* descA,
  int* ipiv, double* B, const int* descb );
void LinearSolve
( int n, int numRhs, scomplex* A, const int* descA,
  int* ipiv, scomplex* B, const int* descb );
void LinearSolve
( int n, int numRhs, dcomplex* A, const int* descA,
  int* ipiv, dcomplex* B, const int* descb );

// Spectral analysis
// =================

// SVD
// ---
void SingularValues
( int m, int n,
  float* A, const int* descA,
  float* s );
void SingularValues
( int m, int n,
  double* A, const int* descA,
  double* s );
void SingularValues
( int m, int n,
  scomplex* A, const int* descA,
  float* s );
void SingularValues
( int m, int n,
  dcomplex* A, const int* descA,
  double* s );

void SVD
( int m, int n,
  float* A, const int* descA,
  float* s,
  float* U, const int* descU,
  float* VH, const int* descVH );
void SVD
( int m, int n,
  double* A, const int* descA,
  double* s,
  double* U, const int* descU,
  double* VH, const int* descVH );
void SVD
( int m, int n,
  scomplex* A, const int* descA,
  float* s,
  scomplex* U, const int* descU,
  scomplex* VH, const int* descVH );
void SVD
( int m, int n,
  dcomplex* A, const int* descA,
  double* s,
  dcomplex* U, const int* descU,
  dcomplex* VH, const int* descVH );

// Hermitian eigenvalue problems
// -----------------------------

// Compute eigenvalues
// ^^^^^^^^^^^^^^^^^^^

// All eigenvalues
// """""""""""""""
void HermitianEig
( char uplo, int n,
  float* A, const int* descA,
  float* w );
void HermitianEig
( char uplo, int n,
  double* A, const int* descA,
  double* w );
void HermitianEig
( char uplo, int n,
  scomplex* A, const int* descA,
  float* w );
void HermitianEig
( char uplo, int n,
  dcomplex* A, const int* descA,
  double* w );

// Floating-point range
// """"""""""""""""""""
int HermitianEig
( char uplo, int n,
  float* A, const int* descA,
  float* w,
  float vl, float vu );
int HermitianEig
( char uplo, int n,
  double* A, const int* descA,
  double* w,
  double vl, double vu );
int HermitianEig
( char uplo, int n,
  scomplex* A, const int* descA,
  float* w,
  float vl, float vu );
int HermitianEig
( char uplo, int n,
  dcomplex* A, const int* descA,
  double* w,
  double vl, double vu );

// Index range
// """""""""""
void HermitianEig
( char uplo, int n,
  float* A, const int* descA,
  float* w,
  int il, int iu );
void HermitianEig
( char uplo, int n,
  double* A, const int* descA,
  double* w,
  int il, int iu );
void HermitianEig
( char uplo, int n,
  scomplex* A, const int* descA,
  float* w,
  int il, int iu );
void HermitianEig
( char uplo, int n,
  dcomplex* A, const int* descA,
  double* w,
  int il, int iu );

// Compute eigenpairs
// ^^^^^^^^^^^^^^^^^^

// All eigenpairs
// """"""""""""""
void HermitianEig
( char uplo, int n,
  float* A, const int* descA,
  float* w, 
  float* Z, const int* descZ );
void HermitianEig
( char uplo, int n,
  double* A, const int* descA,
  double* w, 
  double* Z, const int* descZ );
void HermitianEig
( char uplo, int n,
  scomplex* A, const int* descA,
  float* w, 
  scomplex* Z, const int* descZ );
void HermitianEig
( char uplo, int n,
  dcomplex* A, const int* descA,
  double* w, 
  dcomplex* Z, const int* descZ );

// Floating-point range
// """"""""""""""""""""
int HermitianEig
( char uplo, int n,
  float* A, const int* descA,
  float* w, 
  float* Z, const int* descZ,
  float vl, float vu );
int HermitianEig
( char uplo, int n,
  double* A, const int* descA,
  double* w, 
  double* Z, const int* descZ,
  double vl, double vu );
int HermitianEig
( char uplo, int n,
  scomplex* A, const int* descA,
  float* w, 
  scomplex* Z, const int* descZ,
  float vl, float vu );
int HermitianEig
( char uplo, int n,
  dcomplex* A, const int* descA,
  double* w, 
  dcomplex* Z, const int* descZ,
  double vl, double vu );

// Index range
// """""""""""
void HermitianEig
( char uplo, int n,
  float* A, const int* descA,
  float* w, 
  float* Z, const int* descZ,
  int il, int iu );
void HermitianEig
( char uplo, int n,
  double* A, const int* descA,
  double* w, 
  double* Z, const int* descZ,
  int il, int iu );
void HermitianEig
( char uplo, int n,
  scomplex* A, const int* descA,
  float* w, 
  scomplex* Z, const int* descZ,
  int il, int iu );
void HermitianEig
( char uplo, int n,
  dcomplex* A, const int* descA,
  double* w, 
  dcomplex* Z, const int* descZ,
  int il, int iu );

// Reduction of a generalized HPD EVP to standard form
// ---------------------------------------------------
// NOTE: It is required that B have a positive diagonal

// Two-sided Trsm
// ^^^^^^^^^^^^^^
void TwoSidedTrsm
( char uplo, int n, 
        float* A, const int* descA,
  const float* B, const int* descB );
void TwoSidedTrsm
( char uplo, int n, 
        double* A, const int* descA, 
  const double* B, const int* descB );
void TwoSidedTrsm
( char uplo, int n, 
        scomplex* A, const int* descA,
  const scomplex* B, const int* descB );
void TwoSidedTrsm
( char uplo, int n, 
        dcomplex* A, const int* descA,
  const dcomplex* B, const int* descB );

// Two-sided Trmm
// ^^^^^^^^^^^^^^
void TwoSidedTrmm
( char uplo, int n, 
        float* A, const int* descA,
  const float* B, const int* descB );
void TwoSidedTrmm
( char uplo, int n, 
        double* A, const int* descA,
  const double* B, const int* descB );
void TwoSidedTrmm
( char uplo, int n, 
        scomplex* A, const int* descA,
  const scomplex* B, const int* descB );
void TwoSidedTrmm
( char uplo, int n, 
        dcomplex* A, const int* descA,
  const dcomplex* B, const int* descB );

// Hessenberg Schur decomposition via the QR algorithm
// ---------------------------------------------------
// NOTE: In all of these routines, the matrix needs to be explicitly 
//       upper-Hessenberg before the call, otherwise behavior is unpredictable

void HessenbergSchur
( int n,
  float* H, const int* descH,
  scomplex* w,
  bool fullTriangle=false, bool aed=false );
void HessenbergSchur
( int n,
  double* H, const int* descH,
  dcomplex* w,
  bool fullTriangle=false, bool aed=false );
void HessenbergSchur
( int n,
  scomplex* H, const int* descH,
  scomplex* w,
  bool fullTriangle=false, bool aed=false );
void HessenbergSchur
( int n,
  dcomplex* H, const int* descH,
  dcomplex* w,
  bool fullTriangle=false, bool aed=false );

void HessenbergSchur
( int n,
  float* H, const int* descH,
  scomplex* w, 
  float* Q, const int* descQ,
  bool fullTriangle=true, bool multiplyQ=false, bool aed=false );
void HessenbergSchur
( int n,
  double* H, const int* descH,
  dcomplex* w, 
  double* Q, const int* descQ,
  bool fullTriangle=true, bool multiplyQ=false, bool aed=false );
void HessenbergSchur
( int n,
  scomplex* H, const int* descH,
  scomplex* w, 
  scomplex* Q, const int* descQ,
  bool fullTriangle=true, bool multiplyQ=false, bool aed=false );
void HessenbergSchur
( int n,
  dcomplex* H, const int* descH,
  dcomplex* w, 
  dcomplex* Q, const int* descQ,
  bool fullTriangle=true, bool multiplyQ=false, bool aed=false );

// Hessenberg eigenvalues/pairs
// ----------------------------
void HessenbergEig( int n, float* H, const int* descH, scomplex* w );
void HessenbergEig( int n, double* H, const int* descH, dcomplex* w );
void HessenbergEig( int n, scomplex* H, const int* descH, scomplex* w );
void HessenbergEig( int n, dcomplex* H, const int* descH, dcomplex* w );

// TODO: Compute the eigenvectors

} // namespace scalapack
} // namespace El

#endif // ifdef EL_HAVE_SCALAPACK
#endif // ifndef EL_IMPORTS_SCALAPACK_HPP
