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
#ifndef ELEMENTAL_IMPORTS_LAPACK_HPP
#define ELEMENTAL_IMPORTS_LAPACK_HPP 1

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

} // namespace lapack
} // namespace elem

extern "C" {

// Machine constants
float LAPACK(slamch)( const char* cmach );
double LAPACK(dlamch)( const char* cmach );

// Safe norms
float LAPACK(slapy2)
( const float* alpha, const float* beta );
double LAPACK(dlapy2)
( const double* alpha, const double* beta );
float LAPACK(slapy3)
( const float* alpha, const float* beta, const float* gamma );
double LAPACK(dlapy3)
( const double* alpha, const double* beta, const double* gamma );

// Safely compute a Givens rotation
void LAPACK(slartg)
( const float* phi, const float* gamma,
  float* c, float* s, float* rho );
void LAPACK(dlartg)
( const double* phi, const double* gamma,
  double* c, double* s, double* rho );
void LAPACK(clartg)
( const elem::scomplex* phi, const elem::scomplex* gamma,
  float* c, elem::scomplex* s, elem::scomplex* rho );
void LAPACK(zlartg)
( const elem::dcomplex* phi, const elem::dcomplex* gamma,
  double* c, elem::dcomplex* s, elem::dcomplex* rho );

// Cholesky factorization
void LAPACK(spotrf)
( const char* uplo, const int* n, const float* A, const int* lda,
  int* info );
void LAPACK(dpotrf)
( const char* uplo, const int* n, const double* A, const int* lda,
  int* info );
void LAPACK(cpotrf)
( const char* uplo, const int* n, const elem::scomplex* A, 
  const int* lda, int* info );
void LAPACK(zpotrf)
( const char* uplo, const int* n, const elem::dcomplex* A, 
  const int* lda, int* info );

// LU factorization (with partial pivoting)
void LAPACK(sgetrf)
( const int* m, const int* n, 
  float* A, const int* lda, int* p, int* info );
void LAPACK(dgetrf)
( const int* m, const int* n, 
  double* A, const int* lda, int* p, int* info );
void LAPACK(cgetrf)
( const int* m, const int* n, 
  elem::scomplex* A, const int* lda, int* p, int* info );
void LAPACK(zgetrf)
( const int* m, const int* n, 
  elem::dcomplex* A, const int* lda, int* p, int* info );

// For reducing well-conditioned Hermitian generalized EVP to Hermitian 
// standard form
void LAPACK(ssygst)
( const int* itype, const char* uplo, const int* n,
  float* A, int* lda, const float* B, int* ldb, int* info );
void LAPACK(dsygst)
( const int* itype, const char* uplo, const int* n,
  double* A, int* lda, const double* B, int* ldb, int* info );
void LAPACK(chegst)
( const int* itype, const char* uplo, const int* n,
        elem::scomplex* A, const int* lda, 
  const elem::scomplex* B, const int* ldb, int* info );
void LAPACK(zhegst)
( const int* itype, const char* uplo, const int* n,
        elem::dcomplex* A, const int* lda,
  const elem::dcomplex* B, const int* ldb, int* info );

// Triangular inversion
void LAPACK(strtri)
( const char* uplo, const char* diag, 
  const int* n, const float* A, const int* lda, int* info );
void LAPACK(dtrtri)
( const char* uplo, const char* diag, 
  const int* n, const double* A, const int* lda, int* info );
void LAPACK(ctrtri)
( const char* uplo, const char* diag,
  const int* n, const elem::scomplex* A, const int* lda, int* info );
void LAPACK(ztrtri)
( const char* uplo, const char* diag,
  const int* n, const elem::dcomplex* A, const int* lda, int* info );

} // extern "C"

#endif /* ELEMENTAL_IMPORTS_LAPACK_HPP */

