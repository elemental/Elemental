/*
   Copyright (c) 2009-2011, Jack Poulson
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

namespace elemental {
namespace lapack {

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

void Chol( char uplo, int n, const float* A, int lda );
void Chol( char uplo, int n, const double* A, int lda );
#ifndef WITHOUT_COMPLEX
void Chol( char uplo, int n, const scomplex* A, int lda );
void Chol( char uplo, int n, const dcomplex* A, int lda );
#endif

void Hegst
( int itype, char uplo, 
  int n, float* A, int lda, const float* B, int ldb );
void Hegst
( int itype, char uplo,
  int n, double* A, int lda, const double* B, int ldb );
#ifndef WITHOUT_COMPLEX
void Hegst
( int itype, char uplo,
  int n, scomplex* A, int lda, const scomplex* B, int ldb );
void Hegst
( int itype, char uplo,
  int n, dcomplex* A, int lda, const dcomplex* B, int ldb );
#endif

void LU( int m, int n, float* A, int lda, int* p );
void LU( int m, int n, double* A, int lda, int* p );
#ifndef WITHOUT_COMPLEX
void LU( int m, int n, scomplex* A, int lda, int* p );
void LU( int m, int n, dcomplex* A, int lda, int* p );
#endif

void LQ( int m, int n, float* A, int lda );
void LQ( int m, int n, double* A, int lda );
#ifndef WITHOUT_COMPLEX
void LQ( int m, int n, scomplex* A, int lda, scomplex* t );
void LQ( int m, int n, dcomplex* A, int lda, dcomplex* t );
#endif

void QR( int m, int n, float* A, int lda );
void QR( int m, int n, double* A, int lda );
#ifndef WITHOUT_COMPLEX
void QR( int m, int n, scomplex* A, int lda, scomplex* t );
void QR( int m, int n, dcomplex* A, int lda, dcomplex* t );
#endif

float SafeNorm( float alpha, float beta );
double SafeNorm( double alpha, double beta );
float SafeNorm( float alpha, float beta, float gamma );
double SafeNorm( double alpha, double beta, double gamma );

void SVD
( char UColumns, char VColumns, int m, int n, float* A, int lda,
  float* SigmaDiag, float* U, int ldu, float* VT, int ldvt );
void SVD
( char UColumns, char VColumns, int m, int n, double* A, int lda,
  double* SigmaDiag, double* U, int ldu, double* VT, int ldvt );
#ifndef WITHOUT_COMPLEX
void SVD
( char UColumns, char VColumns, int m, int n, scomplex* A, int lda,
  float* SigmaDiag, scomplex* U, int ldu, scomplex* VT, int ldvt );
void SVD
( char UColumns, char VColumns, int m, int n, dcomplex* A, int lda,
  double* SigmaDiag, dcomplex* U, int ldu, dcomplex* VT, int ldvt );
#endif

void HermitianTridiag( char uplo, int n, float* A, int lda );
void HermitianTridiag( char uplo, int n, double* A, int lda );
#ifndef WITHOUT_COMPLEX
void HermitianTridiag( char uplo, int n, scomplex* A, int lda, scomplex* t );
void HermitianTridiag( char uplo, int n, dcomplex* A, int lda, dcomplex* t );
#endif

void TriangularInverse
( char uplo, char diag, int n, const float* A, int lda );
void TriangularInverse
( char uplo, char diag, int n, const double* A, int lda );
#ifndef WITHOUT_COMPLEX
void TriangularInverse
( char uplo, char diag, int n, const scomplex* A, int lda );
void TriangularInverse
( char uplo, char diag, int n, const dcomplex* A, int lda );
#endif

} // lapack
} // elemental

extern "C" {

// Machine constants
float LAPACK(slamch)( const char* cmach );
double LAPACK(dlamch)( const char* cmach );

// LU factorization
void LAPACK(sgetrf)
( const int* m, const int* n, 
  float* A, const int* lda, int* p, int* info );

void LAPACK(dgetrf)
( const int* m, const int* n, 
  double* A, const int* lda, int* p, int* info );

#ifndef WITHOUT_COMPLEX
void LAPACK(cgetrf)
( const int* m, const int* n, 
  elemental::scomplex* A, const int* lda, int* p, int* info );

void LAPACK(zgetrf)
( const int* m, const int* n, 
  elemental::dcomplex* A, const int* lda, int* p, int* info );
#endif

// LQ factorization
void LAPACK(sgelqf)
( const int* m, const int* n, float* A, const int* lda, 
  float* t, float* work, const int* lwork, int* info );

void LAPACK(dgelqf)
( const int* m, const int* n, double* A, const int* lda,
  double* t, double* work, const int* lwork, int* info );

#ifndef WITHOUT_COMPLEX
void LAPACK(cgelqf)
( const int* m, const int* n, elemental::scomplex* A, const int* lda,
  elemental::scomplex* t, elemental::scomplex* work, const int* lwork, 
  int* info );

void LAPACK(zgelqf)
( const int* m, const int* n, elemental::dcomplex* A, const int* lda,
  elemental::dcomplex* t, elemental::dcomplex* work, const int* lwork, 
  int* info );
#endif

// QR factorization
void LAPACK(sgeqrf)
( const int* m, const int* n, float* A, const int* lda, 
  float* t, float* work, const int* lwork, int* info );

void LAPACK(dgeqrf)
( const int* m, const int* n, double* A, const int* lda,
  double* t, double* work, const int* lwork, int* info );

#ifndef WITHOUT_COMPLEX
void LAPACK(cgeqrf)
( const int* m, const int* n, elemental::scomplex* A, const int* lda,
  elemental::scomplex* t, elemental::scomplex* work, const int* lwork, 
  int* info );

void LAPACK(zgeqrf)
( const int* m, const int* n, elemental::dcomplex* A, const int* lda,
  elemental::dcomplex* t, elemental::dcomplex* work, const int* lwork, 
  int* info );
#endif

// Safe norms
float LAPACK(slapy2)
( const float* alpha, const float* beta );

double LAPACK(dlapy2)
( const double* alpha, const double* beta );

float LAPACK(slapy3)
( const float* alpha, const float* beta, const float* gamma );

double LAPACK(dlapy3)
( const double* alpha, const double* beta, const double* gamma );

// Cholesky factorization
void LAPACK(spotrf)
( const char* uplo, const int* n, const float* A, const int* lda,
  int* info );

void LAPACK(dpotrf)
( const char* uplo, const int* n, const double* A, const int* lda,
  int* info );
    
#ifndef WITHOUT_COMPLEX
void LAPACK(cpotrf)
( const char* uplo, const int* n, const elemental::scomplex* A, 
  const int* lda, int* info );
    
void LAPACK(zpotrf)
( const char* uplo, const int* n, const elemental::dcomplex* A, 
  const int* lda, int* info );
#endif

// Hermitian generalized EVP to hermitian standard EVP
void LAPACK(ssygst)
( const int* itype, const char* uplo, const int* n,
  float* A, int* lda, const float* B, int* ldb, int* info );

void LAPACK(dsygst)
( const int* itype, const char* uplo, const int* n,
  double* A, int* lda, const double* B, int* ldb, int* info );
 
#ifndef WITHOUT_COMPLEX
void LAPACK(chegst)
( const int* itype, const char* uplo, const int* n,
        elemental::scomplex* A, const int* lda, 
  const elemental::scomplex* B, const int* ldb, int* info );

void LAPACK(zhegst)
( const int* itype, const char* uplo, const int* n,
        elemental::dcomplex* A, const int* lda,
  const elemental::dcomplex* B, const int* ldb, int* info );
#endif

// SVD
void LAPACK(sgesvd)
( const char* UColumns, const char* VColumns, const int* m, const int* n,
  float* A, const int* lda, float* SigmaDiag, float* U, const int* ldu,
  float* VT, const int* ldvt, float* work, const int* lwork, int* info );

void LAPACK(dgesvd)
( const char* UColumns, const char* VColumns, const int* m, const int* n,
  double* A, const int* lda, double* SigmaDiag, double* U, const int* ldu,
  double* VT, const int* ldvt, double* work, const int* lwork, int* info );

#ifndef WITHOUT_COMPLEX
void LAPACK(cgesvd)
( const char* UColumns, const char* VColumns, const int* m, const int* n,
  elemental::scomplex* A, const int* lda, float* SigmaDiag, 
  elemental::scomplex* U, const int* ldu, elemental::scomplex* VT, 
  const int* ldvt, elemental::scomplex* work, const int* lwork, 
  float* realWork, int* info );

void LAPACK(zgesvd)
( const char* UColumns, const char* VColumns, const int* m, const int* n,
  elemental::dcomplex* A, const int* lda, double* SigmaDiag, 
  elemental::dcomplex* U, const int* ldu, elemental::dcomplex* VT, 
  const int* ldvt, elemental::dcomplex* work, const int* lwork, 
  double* realWork, int* info );
#endif

// Triangular inversion
void LAPACK(strtri)
( const char* uplo, const char* diag, 
  const int* n, const float* A, const int* lda, int* info );

void LAPACK(dtrtri)
( const char* uplo, const char* diag, 
  const int* n, const double* A, const int* lda, int* info );
    
#ifndef WITHOUT_COMPLEX
void LAPACK(ctrtri)
( const char* uplo, const char* diag,
  const int* n, const elemental::scomplex* A, const int* lda, int* info );
    
void LAPACK(ztrtri)
( const char* uplo, const char* diag,
  const int* n, const elemental::dcomplex* A, const int* lda, int* info );
#endif

// Reduction from hermitian to symmetric tridiagonal
void LAPACK(ssytrd)
( const char* uplo,
  const int* n, float* A, const int* lda,
  float* d, float* e, float* tau, float* work, int* lwork, int* info );

void LAPACK(dsytrd)
( const char* uplo,
  const int* n, double* A, const int* lda, 
  double* d, double* e, double* tau, double* work, int* lwork, int* info );

#ifndef WITHOUT_COMPLEX
void LAPACK(chetrd)
( const char* uplo,
  const int* n, elemental::scomplex* A, const int* lda,
  float* d, float* e, elemental::scomplex* tau, 
  elemental::scomplex* work, int* lwork, int* info );

void LAPACK(zhetrd)
( const char* uplo,
  const int* n, elemental::dcomplex* A, const int* lda,
  double* d, double* e, elemental::dcomplex* tau, 
  elemental::dcomplex* work, int* lwork, int* info );
#endif

} // extern "C"

#endif /* ELEMENTAL_IMPORTS_LAPACK_HPP */

