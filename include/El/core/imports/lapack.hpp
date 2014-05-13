/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_LAPACK_HPP
#define EL_IMPORTS_LAPACK_HPP

namespace El {
namespace lapack {

// Machine constants
// =================

// Relative machine precision
template<typename R> R MachineEpsilon();
template<> float MachineEpsilon<float>();
template<> double MachineEpsilon<double>();

// Minimum number which can be inverted without overflow
template<typename R> R MachineSafeMin();
template<> float MachineSafeMin<float>();
template<> double MachineSafeMin<double>();

// Base of the machine, where the number is represented as 
//   (mantissa) x (base)^(exponent)
template<typename R> R MachineBase();
template<> float MachineBase<float>();
template<> double MachineBase<double>();

// Return the relative machine precision multiplied by the base
template<typename R> R MachinePrecision();
template<> float MachinePrecision<float>();
template<> double MachinePrecision<double>();

// Return the minimum exponent before (gradual) underflow occurs
template<typename R> R MachineUnderflowExponent();
template<> float MachineUnderflowExponent<float>();
template<> double MachineUnderflowExponent<double>();

// Return the underflow threshold: (base)^((underflow exponent)-1)
template<typename R> R MachineUnderflowThreshold();
template<> float MachineUnderflowThreshold<float>();
template<> double MachineUnderflowThreshold<double>();

// Return the largest exponent before overflow
template<typename R> R MachineOverflowExponent();
template<> float MachineOverflowExponent<float>();
template<> double MachineOverflowExponent<double>();

// Return the overflow threshold: (1-(rel. prec.)) * (base)^(overflow exponent)
template<typename R> R MachineOverflowThreshold();
template<> float MachineOverflowThreshold<float>();
template<> double MachineOverflowThreshold<double>();

// For safely computing norms without overflow/underflow
// =====================================================

float SafeNorm( float alpha, float beta );
double SafeNorm( double alpha, double beta );
float SafeNorm( float alpha, float beta, float gamma );
double SafeNorm( double alpha, double beta, double gamma );

float SafeNorm( Complex<float> alpha, float beta );
double SafeNorm( Complex<double> alpha, double beta );
float SafeNorm( float alpha, Complex<float> beta );
double SafeNorm( double alpha, Complex<double> beta );

// Givens rotations
// ================
//
// Given phi and gamma, compute a Givens rotation such that
//
//  |       c   s | |   phi |  = | rho |, where c^2 + |s|^2 = 1
//  | -conj(s)  c | | gamma |    |  0  |
//
// This routine uses the stable approach suggested by Kahan and Demmel and
// returns the value rho.
//

float    Givens( float    phi, float    gamma, float * c, float   * s );
double   Givens( double   phi, double   gamma, double* c, double  * s );
scomplex Givens( scomplex phi, scomplex gamma, float * c, scomplex* s );
dcomplex Givens( dcomplex phi, dcomplex gamma, double* c, dcomplex* s );

// Compute the eigen-values/pairs of a symmetric tridiagonal matrix
// ================================================================

// Compute eigenvalues
// -------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void SymmetricTridiagEig
( int n, float* d, float* e, float* w, float abstol=0 );
void SymmetricTridiagEig
( int n, double* d, double* e, double* w, double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int SymmetricTridiagEig
( int n, float* d, float* e, float* w, 
  float vl, float vu, float abstol=0 );
int SymmetricTridiagEig
( int n, double* d, double* e, double* w, 
  double vl, double vu, double abstol=0 );

// Index range
// ^^^^^^^^^^^
void SymmetricTridiagEig
( int n, float* d, float* e, float* w,
  int il, int iu, float abstol=0 );
void SymmetricTridiagEig
( int n, double* d, double* e, double* w,
  int il, int iu, double abstol=0 );

// Compute eigenpairs
// ------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void SymmetricTridiagEig
( int n, float* d, float* e, float* w, float* Z, int ldZ, float abstol=0 );
void SymmetricTridiagEig
( int n, double* d, double* e, double* w, double* Z, int ldZ, double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int SymmetricTridiagEig
( int n, float* d, float* e, float* w, float* Z, int ldZ,
  float vl, float vu, float abstol=0 );
int SymmetricTridiagEig
( int n, double* d, double* e, double* w, double* Z, int ldZ,
  double vl, double vu, double abstol=0 );

// Index range
// ^^^^^^^^^^^
void SymmetricTridiagEig
( int n, float* d, float* e, float* w, float* Z, int ldZ,
  int il, int iu, float abstol=0 );
void SymmetricTridiagEig
( int n, double* d, double* e, double* w, double* Z, int ldZ,
  int il, int iu, double abstol=0 );

// Compute the eigen-values/pairs of a Hermitian matrix
// ====================================================

// Compute eigenvalues
// -------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, int ldA, float* w, float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int HermitianEig
( char uplo, int n, float* A, int ldA, float* w,
  float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, double* A, int ldA, double* w,
  double vl, double vu, double abstol=0 );
int HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w,
  float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w,
  double vl, double vu, double abstol=0 );

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, int ldA, float* w,
  int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w,
  int il, int iu, double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w,
  int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w,
  int il, int iu, double abstol=0 );

// Compute eigenpairs
// ------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, int ldA, float* w, float* Z, int ldZ,
  float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double* Z, int ldZ,
  double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, scomplex* Z, int ldZ,
  float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, dcomplex* Z, int ldZ,
  double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
int HermitianEig
( char uplo, int n, float* A, int ldA, float* w, float* Z, int ldZ,
  float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double* Z, int ldZ,
  double vl, double vu, double abstol=0 );
int HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, scomplex* Z, int ldZ,
  float vl, float vu, float abstol=0 );
int HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, dcomplex* Z, int ldZ,
  double vl, double vu, double abstol=0 );

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, int n, float* A, int ldA, float* w, float* Z, int ldZ,
  int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, double* A, int ldA, double* w, double* Z, int ldZ,
  int il, int iu, double abstol=0 );
void HermitianEig
( char uplo, int n, scomplex* A, int ldA, float* w, scomplex* Z, int ldZ,
  int il, int iu, float abstol=0 );
void HermitianEig
( char uplo, int n, dcomplex* A, int ldA, double* w, dcomplex* Z, int ldZ,
  int il, int iu, double abstol=0 );

// Compute the SVD of a general matrix using a divide and conquer algorithm
// ========================================================================

void DivideAndConquerSVD
( int m, int n, float* A, int ldA, 
  float* s, float* U, int ldU, float* VT, int ldVT );
void DivideAndConquerSVD
( int m, int n, double* A, int ldA, 
  double* s, double* U, int ldU, double* VT, int ldVT );
void DivideAndConquerSVD
( int m, int n, scomplex* A, int ldA, 
  float* s, scomplex* U, int ldU, scomplex* VH, int ldVH );
void DivideAndConquerSVD
( int m, int n, dcomplex* A, int ldA, 
  double* s, dcomplex* U, int ldU, dcomplex* VH, int ldVH );

// Compute the SVD of a general matrix using the QR algorithm
// ==========================================================

void QRSVD
( int m, int n, float* A, int ldA, 
  float* s, float* U, int ldU, float* VT, int ldVT );
void QRSVD
( int m, int n, double* A, int ldA, 
  double* s, double* U, int ldU, double* VT, int ldVT );
void QRSVD
( int m, int n, scomplex* A, int ldA, 
  float* s, scomplex* U, int ldU, scomplex* VH, int ldVH );
void QRSVD
( int m, int n, dcomplex* A, int ldA, 
  double* s, dcomplex* U, int ldU, dcomplex* VH, int ldVH );

// Compute the singular values of a general matrix (using the QR algorithm)
// ========================================================================

void SVD( int m, int n, float* A, int ldA, float* s );
void SVD( int m, int n, double* A, int ldA, double* s );
void SVD( int m, int n, scomplex* A, int ldA, float* s );
void SVD( int m, int n, dcomplex* A, int ldA, double* s );

// Compute the singular values of a bidiagonal matrix via dqds
// ===========================================================

void BidiagDQDS( int n, float* d, float* e );
void BidiagDQDS( int n, double* d, double* e );

// Compute the SVD of a bidiagonal matrix using the QR algorithm
// =============================================================

void BidiagQRAlg
( char uplo, int n, int numColsVT, int numRowsU,
  float* d, float* e, float* VT, int ldVT, float* U, int ldU );
void BidiagQRAlg
( char uplo, int n, int numColsVT, int numRowsU, 
  double* d, double* e, double* VT, int ldVT, double* U, int ldU );
void BidiagQRAlg
( char uplo, int n, int numColsVH, int numRowsU,
  float* d, float* e, scomplex* VH, int ldVH, scomplex* U, int ldU );
void BidiagQRAlg
( char uplo, int n, int numColsVH, int numRowsU, 
  double* d, double* e, dcomplex* VH, int ldVH, dcomplex* U, int ldU );

// Compute the Schur decomposition of an upper Hessenberg matrix
// =============================================================

void HessenbergSchur
( int n, float* H, int ldH, scomplex* w, bool fullTriangle=false );
void HessenbergSchur
( int n, double* H, int ldH, dcomplex* w, bool fullTriangle=false );
void HessenbergSchur
( int n, scomplex* H, int ldH, scomplex* w, bool fullTriangle=false );
void HessenbergSchur
( int n, dcomplex* H, int ldH, dcomplex* w, bool fullTriangle=false );

void HessenbergSchur
( int n, float* H, int ldH, scomplex* w, float* Q, int ldQ, 
  bool fullTriangle=true, bool multiplyQ=false );
void HessenbergSchur
( int n, double* H, int ldH, dcomplex* w, double* Q, int ldQ, 
  bool fullTriangle=true, bool multiplyQ=false );
void HessenbergSchur
( int n, scomplex* H, int ldH, scomplex* w, scomplex* Q, int ldQ, 
  bool fullTriangle=false, bool multiplyQ=false );
void HessenbergSchur
( int n, dcomplex* H, int ldH, dcomplex* w, dcomplex* Q, int ldQ, 
  bool fullTriangle=false, bool multiplyQ=false );

// Compute the eigenvalues/pairs of an upper Hessenberg matrix
// ===========================================================

void HessenbergEig( int n, float* H, int ldH, scomplex* w );
void HessenbergEig( int n, double* H, int ldH, dcomplex* w );
void HessenbergEig( int n, scomplex* H, int ldH, scomplex* w );
void HessenbergEig( int n, dcomplex* H, int ldH, dcomplex* w );

// TODO: A version which computes eigenvectors

// Compute the Schur decomposition of a square matrix
// ==================================================

void Schur
( int n, float* A, int ldA, scomplex* w, bool fullTriangle=false );
void Schur
( int n, double* A, int ldA, dcomplex* w, bool fullTriangle=false );
void Schur
( int n, scomplex* A, int ldA, scomplex* w, bool fullTriangle=false );
void Schur
( int n, dcomplex* A, int ldA, dcomplex* w, bool fullTriangle=false );

void Schur
( int n, float* A, int ldA, scomplex* w, float* Q, int ldQ, 
  bool fullTriangle=true );
void Schur
( int n, double* A, int ldA, dcomplex* w, double* Q, int ldQ, 
  bool fullTriangle=true );
void Schur
( int n, scomplex* A, int ldA, scomplex* w, scomplex* Q, int ldQ, 
  bool fullTriangle=true );
void Schur
( int n, dcomplex* A, int ldA, dcomplex* w, dcomplex* Q, int ldQ, 
  bool fullTriangle=true );

// Compute the eigenvalues/pairs of a square matrix
// ================================================

void Eig( int n, float* A, int ldA, scomplex* w );
void Eig( int n, double* A, int ldA, dcomplex* w );
void Eig( int n, scomplex* A, int ldA, scomplex* w );
void Eig( int n, dcomplex* A, int ldA, dcomplex* w );

void Eig( int n, float* A, int ldA, scomplex* w, scomplex* X, int ldX );
void Eig( int n, float* A, int ldA, scomplex* w, float* XPacked, int ldX );
void Eig( int n, double* A, int ldA, dcomplex* w, dcomplex* X, int ldX );
void Eig( int n, double* A, int ldA, dcomplex* w, double* XPacked, int ldX );
void Eig( int n, scomplex* A, int ldA, scomplex* w, scomplex* X, int ldX );
void Eig( int n, dcomplex* A, int ldA, dcomplex* w, dcomplex* X, int ldX );

} // namespace lapack
} // namespace El

#endif // ifndef EL_IMPORTS_LAPACK_HPP
