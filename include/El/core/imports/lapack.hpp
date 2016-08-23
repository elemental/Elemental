/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_LAPACK_HPP
#define EL_IMPORTS_LAPACK_HPP

namespace El {

namespace lapack {

// For copying column-major matrices
// =================================
template<typename T>
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const T* A, BlasInt lda, T* B, BlasInt ldb );

void Copy
( char uplo, BlasInt m, BlasInt n, 
  const float* A, BlasInt lda, float* B, BlasInt ldb );
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const double* A, BlasInt lda, double* B, BlasInt ldb );
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb );
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb );
template<typename T>
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const T* A, BlasInt lda, T* B, BlasInt ldb );

// Generate a Householder reflector
// ================================
// NOTE: 
// Since LAPACK chooses to use the identity matrix, rather than a single
// coordinate negation, in cases where the mass is already entirely in the
// first entry, and the identity matrix cannot be represented as a Householder
// reflector, Elemental does not ever directly call LAPACK's Householder
// routines. Otherwise, the logic of routines such as ApplyPackedReflectors
// would need to be (unnecessarily) complicated.
//
// Furthermore, LAPACK defines H = I - tau [1; v] [1; v]' such that
// adjoint(H) [chi; x] = [beta; 0], but Elemental instead defines
// H = I - tau [1; v] [1; v]' such that H [chi; x] = [beta; 0].
//
template<typename F>
F Reflector( BlasInt n, F& chi, F* x, BlasInt incx );

// Apply a Householder reflector
// =============================
// Form either
//
//     C := (I - tau v v^H) C
//
// or,
//
//     C := C (I - tau v v^H)
//
// Note that 'work' must be of size n if onLeft and of size m if !onLeft.
//
template<typename F>
void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const F* v, BlasInt vInc, 
  const F& tau,
        F* C, BlasInt CLDim );
template<typename F>
void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const F* v, BlasInt vInc, 
  const F& tau,
        F* C, BlasInt CLDim,
        F* work );

// Compute the eigen-values/pairs of a symmetric tridiagonal matrix
// ================================================================

// Compute eigenvalues
// -------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float abstol=0 );
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, 
  float vl, float vu, float abstol=0 );
BlasInt SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, 
  double vl, double vu, double abstol=0 );

// Index range
// ^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w,
  BlasInt il, BlasInt iu, float abstol=0 );
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w,
  BlasInt il, BlasInt iu, double abstol=0 );

// Compute eigenpairs
// ------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, 
  float* d, float* e, float* w, float* Z, BlasInt ldZ, float abstol=0 );
void SymmetricTridiagEig
( BlasInt n, 
  double* d, double* e, double* w, double* Z, BlasInt ldZ, double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ,
  float vl, float vu, float abstol=0 );
BlasInt SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ,
  double vl, double vu, double abstol=0 );

// Index range
// ^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, float abstol=0 );
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, double abstol=0 );

// Compute the eigen-values/pairs of a Hermitian matrix
// ====================================================

// Compute eigenvalues
// -------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w, double abstol=0 );
void HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w, float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w, double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w,
  float vl, float vu, float abstol=0 );
BlasInt HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w,
  double vl, double vu, double abstol=0 );
BlasInt HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w,
  float vl, float vu, float abstol=0 );
BlasInt HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w,
  double vl, double vu, double abstol=0 );

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w,
  BlasInt il, BlasInt iu, float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w,
  BlasInt il, BlasInt iu, double abstol=0 );
void HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w,
  BlasInt il, BlasInt iu, float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w,
  BlasInt il, BlasInt iu, double abstol=0 );

// Compute eigenpairs
// ------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, 
  float* A, BlasInt ldA,
  float* w,
  float* Z, BlasInt ldZ,
  float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA,
  double* w,
  double* Z, BlasInt ldZ,
  double abstol=0 );
void HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA,
  float* w,
  scomplex* Z, BlasInt ldZ,
  float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA,
  double* w,
  dcomplex* Z, BlasInt ldZ,
  double abstol=0 );

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt HermitianEig
( char uplo, BlasInt n, 
  float* A, BlasInt ldA,
  float* w,
  float* Z, BlasInt ldZ,
  float vl, float vu,
  float abstol=0 );
BlasInt HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA,
  double* w,
  double* Z, BlasInt ldZ,
  double vl, double vu,
  double abstol=0 );
BlasInt HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA,
  float* w,
  scomplex* Z, BlasInt ldZ,
  float vl, float vu,
  float abstol=0 );
BlasInt HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA,
  double* w,
  dcomplex* Z, BlasInt ldZ,
  double vl, double vu,
  double abstol=0 );

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, 
  float* A, BlasInt ldA,
  float* w,
  float* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu,
  float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA,
  double* w,
  double* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu,
  double abstol=0 );
void HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA,
  float* w,
  scomplex* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu,
  float abstol=0 );
void HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA,
  double* w,
  dcomplex* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu,
  double abstol=0 );

// Compute the SVD of a general matrix using a divide and conquer algorithm
// ========================================================================

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  float* A, BlasInt ldA, 
  float* s,
  float* U, BlasInt ldU,
  float* VT, BlasInt ldVT,
  bool thin=true );
void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  double* A, BlasInt ldA, 
  double* s,
  double* U, BlasInt ldU,
  double* VT, BlasInt ldVT,
  bool thin=true );
void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  scomplex* A, BlasInt ldA, 
  float* s,
  scomplex* U, BlasInt ldU,
  scomplex* VH, BlasInt ldVH,
  bool thin=true );
void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  dcomplex* A, BlasInt ldA, 
  double* s,
  dcomplex* U, BlasInt ldU,
  dcomplex* VH, BlasInt ldVH,
  bool thin=true );

// Compute the SVD of a general matrix using the QR algorithm
// ==========================================================

void QRSVD
( BlasInt m, BlasInt n,
  float* A, BlasInt ldA, 
  float* s,
  float* U, BlasInt ldU,
  float* VT, BlasInt ldVT,
  bool thin=true, bool avoidU=false, bool avoidV=false );
void QRSVD
( BlasInt m, BlasInt n,
  double* A, BlasInt ldA, 
  double* s,
  double* U, BlasInt ldU,
  double* VT, BlasInt ldVT,
  bool thin=true, bool avoidU=false, bool avoidV=false );
void QRSVD
( BlasInt m, BlasInt n,
  scomplex* A, BlasInt ldA, 
  float* s,
  scomplex* U, BlasInt ldU,
  scomplex* VH, BlasInt ldVH,
  bool thin=true, bool avoidU=false, bool avoidV=false );
void QRSVD
( BlasInt m, BlasInt n,
  dcomplex* A, BlasInt ldA, 
  double* s,
  dcomplex* U, BlasInt ldU,
  dcomplex* VH, BlasInt ldVH,
  bool thin=true, bool avoidU=false, bool avoidV=false );

// Compute the singular values of a general matrix (using the QR algorithm)
// ========================================================================

void SVD( BlasInt m, BlasInt n, float* A, BlasInt ldA, float* s );
void SVD( BlasInt m, BlasInt n, double* A, BlasInt ldA, double* s );
void SVD( BlasInt m, BlasInt n, scomplex* A, BlasInt ldA, float* s );
void SVD( BlasInt m, BlasInt n, dcomplex* A, BlasInt ldA, double* s );

// Compute the singular values of a bidiagonal matrix via dqds
// ===========================================================

void BidiagDQDS( BlasInt n, float* d, float* e );
void BidiagDQDS( BlasInt n, double* d, double* e );

// Compute the SVD of a bidiagonal matrix using the QR algorithm
// =============================================================

void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVT, BlasInt numRowsU,
  float* d, float* e, float* VT, BlasInt ldVT, float* U, BlasInt ldU );
void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVT, BlasInt numRowsU, 
  double* d, double* e, double* VT, BlasInt ldVT, double* U, BlasInt ldU );
void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVH, BlasInt numRowsU,
  float* d, float* e, scomplex* VH, BlasInt ldVH, scomplex* U, BlasInt ldU );
void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVH, BlasInt numRowsU, 
  double* d, double* e, dcomplex* VH, BlasInt ldVH, dcomplex* U, BlasInt ldU );

// Reduce a general square matrix to upper Hessenberg form
// =======================================================

void Hessenberg
( BlasInt n,
  float* A, BlasInt ldA,
  float* tau );
void Hessenberg
( BlasInt n,
  double* A, BlasInt ldA,
  double* tau );
void Hessenberg
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* tau );
void Hessenberg
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* tau );

// Return the unitary matrix from the reduction to Hessenberg form
void HessenbergGenerateUnitary
( BlasInt n,
  float* A, BlasInt ldA,
  const float* tau );
void HessenbergGenerateUnitary
( BlasInt n,
  double* A, BlasInt ldA,
  const double* tau );
void HessenbergGenerateUnitary
( BlasInt n,
  scomplex* A, BlasInt ldA,
  const scomplex* tau );
void HessenbergGenerateUnitary
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  const dcomplex* tau );

// Solve a 2x2 linear system using LU with full pivoting, perturbing as
// necessary to ensure sufficiently large pivots
// ==============================================================
// (NOTE: This is primarily a helper function for SmallSylvester)
template<typename Real>
bool Solve2x2FullPiv
( const Real* A,
        Real* b,
        Real& scale,
  const Real& smallNum,
  const Real& minPiv );

// Solve a 4x4 linear system using LU with full pivoting, perturbing as
// necessary to ensure sufficiently large pivots
// ==============================================================
// (NOTE: This is primarily a helper function for SmallSylvester)
template<typename Real>
bool Solve4x4FullPiv
(       Real* A,
        Real* b,
        Real& scale,
  const Real& smallNum,
  const Real& minPiv );

// Small Sylvester
// ===============
// Solve a 1x1, 1x2, 2x1, or 2x2 real Sylvester equation, 
//
//   op_C(C) X +- X op_D(D) = scale*B,
//
// where op_C(C) is either C or C^T, op_D(D) is either D or D^T,
// and scale in (0,1] is determined by the subroutine.
//
// The fundamental technique is Gaussian Elimination with full pivoting,
// with pivots forced to be sufficiently large (and, if such a perturbation was
// performed, the routine returns 'true').
//
// The analogous LAPACK routines are {s,d}lasy2.
//
template<typename Real>
bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const Real* C, BlasInt CLDim,
  const Real* D, BlasInt DLDim,
  const Real* B, BlasInt BLDim,
        Real& scale,
        Real* X, BlasInt XLDim,
        Real& XInfNorm );

// Schur exchange
// ==============
// NOTE: In all cases below, 'work' is expected to be a preallocated buffer
//       of length n
template<typename Real>
void AdjacentSchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim, 
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  Real* work,
  bool testAccuracy=true );
template<typename Real>
void AdjacentSchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim, 
  Real* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  Real* work,
  bool testAccuracy=true );

template<typename Real>
void SchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim, 
  BlasInt j1,
  BlasInt j2,
  Real* work,
  bool testAccuracy=true );
template<typename Real>
void SchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim, 
  Real* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  Real* work,
  bool testAccuracy=true );

template<typename Real>
void SchurExchange
( BlasInt n,
  Complex<Real>* T, BlasInt TLDim, 
  BlasInt j1,
  BlasInt j2 );
template<typename Real>
void SchurExchange
( BlasInt n,
  Complex<Real>* T, BlasInt TLDim, 
  Complex<Real>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 );

// Compute the Schur decomposition of an upper Hessenberg matrix
// =============================================================

void HessenbergSchur
( BlasInt n,
  float* H, BlasInt ldH,
  scomplex* w,
  bool fullTriangle=false,
  bool useAED=true );
void HessenbergSchur
( BlasInt n,
  double* H, BlasInt ldH,
  dcomplex* w,
  bool fullTriangle=false,
  bool useAED=true );
void HessenbergSchur
( BlasInt n,
  scomplex* H, BlasInt ldH,
  scomplex* w,
  bool fullTriangle=false,
  bool useAED=true );
void HessenbergSchur
( BlasInt n,
  dcomplex* H, BlasInt ldH,
  dcomplex* w,
  bool fullTriangle=false,
  bool useAED=true );

void HessenbergSchur
( BlasInt n,
  float* H, BlasInt ldH,
  scomplex* w,
  float* Q, BlasInt ldQ, 
  bool fullTriangle=true,
  bool multiplyQ=false,
  bool useAED=true );
void HessenbergSchur
( BlasInt n,
  double* H, BlasInt ldH,
  dcomplex* w,
  double* Q, BlasInt ldQ, 
  bool fullTriangle=true,
  bool multiplyQ=false,
  bool useAED=true );
void HessenbergSchur
( BlasInt n,
  scomplex* H, BlasInt ldH,
  scomplex* w,
  scomplex* Q, BlasInt ldQ, 
  bool fullTriangle=false,
  bool multiplyQ=false,
  bool useAED=true );
void HessenbergSchur
( BlasInt n,
  dcomplex* H, BlasInt ldH,
  dcomplex* w,
  dcomplex* Q, BlasInt ldQ, 
  bool fullTriangle=false,
  bool multiplyQ=false,
  bool useAED=true );

// Compute the eigenvalues/pairs of an upper Hessenberg matrix
// ===========================================================

void HessenbergEig( BlasInt n, float* H, BlasInt ldH, scomplex* w );
void HessenbergEig( BlasInt n, double* H, BlasInt ldH, dcomplex* w );
void HessenbergEig( BlasInt n, scomplex* H, BlasInt ldH, scomplex* w );
void HessenbergEig( BlasInt n, dcomplex* H, BlasInt ldH, dcomplex* w );

// TODO: A version which computes eigenvectors

// Compute the Schur decomposition of a square matrix
// ==================================================

void Schur
( BlasInt n,
  float* A, BlasInt ldA,
  scomplex* w,
  bool fullTriangle=false,
  bool time=false );
void Schur
( BlasInt n,
  double* A, BlasInt ldA,
  dcomplex* w,
  bool fullTriangle=false,
  bool time=false );
void Schur
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* w,
  bool fullTriangle=false,
  bool time=false );
void Schur
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* w,
  bool fullTriangle=false,
  bool time=false );

void Schur
( BlasInt n,
  float* A, BlasInt ldA,
  scomplex* w,
  float* Q, BlasInt ldQ, 
  bool fullTriangle=true,
  bool time=false );
void Schur
( BlasInt n,
  double* A, BlasInt ldA,
  dcomplex* w,
  double* Q, BlasInt ldQ, 
  bool fullTriangle=true,
  bool time=false );
void Schur
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* w,
  scomplex* Q, BlasInt ldQ, 
  bool fullTriangle=true,
  bool time=false );
void Schur
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* w,
  dcomplex* Q, BlasInt ldQ, 
  bool fullTriangle=true,
  bool time=false );

// Compute eigenvectors of an upper (quasi-)triangular matrix
// ==========================================================

void QuasiTriangEig
( BlasInt n,
  float* T, BlasInt ldT,
  float* VR, BlasInt ldVR,
  bool accumulate=false );
void QuasiTriangEig
( BlasInt n,
  double* T, BlasInt ldT,
  double* VR, BlasInt ldVR,
  bool accumulate=false );

void TriangEig
( BlasInt n,
  scomplex* T, BlasInt ldT,
  scomplex* VR, BlasInt ldVR,
  bool accumulate=false );
void TriangEig
( BlasInt n,
  dcomplex* T, BlasInt ldT,
  dcomplex* VR, BlasInt ldVR,
  bool accumulate=false );

// Compute the eigenvalues/pairs of a square matrix
// ================================================

// TODO: Optional support for balancing

void Eig
( BlasInt n,
  float* A, BlasInt ldA,
  scomplex* w,
  bool time=false );
void Eig
( BlasInt n,
  double* A, BlasInt ldA,
  dcomplex* w,
  bool time=false );
void Eig
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* w,
  bool time=false );
void Eig
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* w,
  bool time=false );

void Eig
( BlasInt n, 
  float* A, BlasInt ldA,
  scomplex* w,
  scomplex* X, BlasInt ldX,
  bool time=false );
void Eig
( BlasInt n, 
  float* A, BlasInt ldA,
  scomplex* w,
  float* XPacked, BlasInt ldX,
  bool time=false );
void Eig
( BlasInt n, 
  double* A, BlasInt ldA,
  dcomplex* w,
  dcomplex* X, BlasInt ldX,
  bool time=false );
void Eig
( BlasInt n, 
  double* A, BlasInt ldA,
  dcomplex* w,
  double* XPacked, BlasInt ldX,
  bool time=false );
void Eig
( BlasInt n, 
  scomplex* A, BlasInt ldA,
  scomplex* w,
  scomplex* X, BlasInt ldX,
  bool time=false );
void Eig
( BlasInt n, 
  dcomplex* A, BlasInt ldA,
  dcomplex* w,
  dcomplex* X, BlasInt ldX,
  bool time=false );

} // namespace lapack

} // namespace El

#endif // ifndef EL_IMPORTS_LAPACK_HPP
