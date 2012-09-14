/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {
namespace lapack {

//
// Machine constants
//

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

//
// For safely computing norms without overflow/underflow
//

float SafeNorm( float alpha, float beta );
double SafeNorm( double alpha, double beta );
float SafeNorm( float alpha, float beta, float gamma );
double SafeNorm( double alpha, double beta, double gamma );

//
// Given phi and gamma, compute a Givens rotation such that
//
//  |       cs   sn | |   phi |  = | rho |, where cs^2 + |sn|^2 = 1
//  | -conj(sn)  cs | | gamma |    |  0  |
//
// This routine should use the stable approach suggested by Kahan and Demmel
//

void ComputeGivens
( float phi, float gamma,
  float* cs, float* sn, float* rho );

void ComputeGivens
( double phi, double gamma,
  double* cs, double* sn, double* rho );

void ComputeGivens
( scomplex phi, scomplex gamma,
  float* cs, scomplex* sn, scomplex* rho );

void ComputeGivens
( dcomplex phi, dcomplex gamma,
  double* cs, dcomplex* sn, dcomplex* rho );

//
// Cholesky factorization
//

void Cholesky( char uplo, int n, const float* A, int lda );
void Cholesky( char uplo, int n, const double* A, int lda );
void Cholesky( char uplo, int n, const scomplex* A, int lda );
void Cholesky( char uplo, int n, const dcomplex* A, int lda );

//
// LU factorization (with partial pivoting)
//

void LU( int m, int n, float* A, int lda, int* p );
void LU( int m, int n, double* A, int lda, int* p );
void LU( int m, int n, scomplex* A, int lda, int* p );
void LU( int m, int n, dcomplex* A, int lda, int* p );

//
// For reducing well-conditioned Hermitian generalized-definite EVP's
// to standard form.
//

void Hegst
( int itype, char uplo, 
  int n, float* A, int lda, const float* B, int ldb );
void Hegst
( int itype, char uplo,
  int n, double* A, int lda, const double* B, int ldb );
void Hegst
( int itype, char uplo,
  int n, scomplex* A, int lda, const scomplex* B, int ldb );
void Hegst
( int itype, char uplo,
  int n, dcomplex* A, int lda, const dcomplex* B, int ldb );

//
// For computing the inverse of a triangular matrix
//

void TriangularInverse
( char uplo, char diag, int n, const float* A, int lda );
void TriangularInverse
( char uplo, char diag, int n, const double* A, int lda );
void TriangularInverse
( char uplo, char diag, int n, const scomplex* A, int lda );
void TriangularInverse
( char uplo, char diag, int n, const dcomplex* A, int lda );

//
// Compute the SVD of a general matrix using a divide and conquer algorithm
//

void DivideAndConquerSVD
( int m, int n, float* A, int lda, 
  float* s, float* U, int ldu, float* VTrans, int ldvt );
void DivideAndConquerSVD
( int m, int n, double* A, int lda, 
  double* s, double* U, int ldu, double* VTrans, int ldvt );
void DivideAndConquerSVD
( int m, int n, scomplex* A, int lda, 
  float* s, scomplex* U, int ldu, scomplex* VAdj, int ldva );
void DivideAndConquerSVD
( int m, int n, dcomplex* A, int lda, 
  double* s, dcomplex* U, int ldu, dcomplex* VAdj, int ldva );

//
// Compute the SVD of a general matrix using the QR algorithm
//

void QRSVD
( int m, int n, float* A, int lda, 
  float* s, float* U, int ldu, float* VTrans, int ldvt );
void QRSVD
( int m, int n, double* A, int lda, 
  double* s, double* U, int ldu, double* VTrans, int ldvt );
void QRSVD
( int m, int n, scomplex* A, int lda, 
  float* s, scomplex* U, int ldu, scomplex* VAdj, int ldva );
void QRSVD
( int m, int n, dcomplex* A, int lda, 
  double* s, dcomplex* U, int ldu, dcomplex* VAdj, int ldva );

//
// Compute the singular values of a general matrix using the QR algorithm
//

void SingularValues( int m, int n, float* A, int lda, float* s );
void SingularValues( int m, int n, double* A, int lda, double* s );
void SingularValues( int m, int n, scomplex* A, int lda, float* s );
void SingularValues( int m, int n, dcomplex* A, int lda, double* s );

//
// Compute the SVD of a bidiagonal matrix using the QR algorithm
//

void BidiagQRAlg
( char uplo, int n, int numColsVTrans, int numRowsU,
  float* d, float* e, float* VTrans, int ldVTrans, float* U, int ldU );
void BidiagQRAlg
( char uplo, int n, int numColsVTrans, int numRowsU, 
  double* d, double* e, double* VTrans, int ldVTrans, double* U, int ldU );
void BidiagQRAlg
( char uplo, int n, int numColsVAdj, int numRowsU,
  float* d, float* e, scomplex* VAdj, int ldVAdj, scomplex* U, int ldU );
void BidiagQRAlg
( char uplo, int n, int numColsVAdj, int numRowsU, 
  double* d, double* e, dcomplex* VAdj, int ldVAdj, dcomplex* U, int ldU );

// 
// Compute the eigenvalues of an upper Hessenberg matrix
//

void HessenbergEig( int n, float* H, int ldh, scomplex* w );
void HessenbergEig( int n, double* H, int ldh, dcomplex* w );
void HessenbergEig( int n, scomplex* H, int ldh, scomplex* w );
void HessenbergEig( int n, dcomplex* H, int ldh, dcomplex* w );

} // namespace lapack
} // namespace elem
