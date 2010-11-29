/*
   Copyright (c) 2009-2010, Jack Poulson
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
#ifndef ELEMENTAL_WRAPPERS_LAPACK_HPP
#define ELEMENTAL_WRAPPERS_LAPACK_HPP 1

#if defined(LAPACK_PRE) && defined(LAPACK_POST)
#define LAPACK(name) _ ## name ## _
#elif defined(LAPACK_PRE)
#define LAPACK(name) _ ## name
#elif defined(LAPACK_POST)
#define LAPACK(name) name ## _
#else
#define LAPACK(name) name
#endif

namespace elemental {
namespace wrappers {
namespace lapack {

void
Chol
( char uplo, int n, const float* A, int lda );

void
Chol
( char uplo, int n, const double* A, int lda );
 
#ifndef WITHOUT_COMPLEX
void
Chol
( char uplo, int n, const scomplex* A, int lda );

void
Chol
( char uplo, int n, const dcomplex* A, int lda );
#endif

void
Hegst
( int itype, char uplo, 
  int n, float* A, int lda, const float* B, int ldb );

void
Hegst
( int itype, char uplo,
  int n, double* A, int lda, const double* B, int ldb );

#ifndef WITHOUT_COMPLEX
void
Hegst
( int itype, char uplo,
  int n, scomplex* A, int lda, const scomplex* B, int ldb );

void
Hegst
( int itype, char uplo,
  int n, dcomplex* A, int lda, const dcomplex* B, int ldb );
#endif

void
LU
( int m, int n, float* A, int lda, int* p );

void
LU
( int m, int n, double* A, int lda, int* p );

#ifndef WITHOUT_COMPLEX
void
LU
( int m, int n, scomplex* A, int lda, int* p );

void
LU
( int m, int n, dcomplex* A, int lda, int* p );
#endif

void
QR
( int m, int n, float* A, int lda );

void
QR
( int m, int n, double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
QR
( int m, int n, scomplex* A, int lda, scomplex* t );

void
QR
( int m, int n, dcomplex* A, int lda, dcomplex* t );
#endif

float
SafeNorm
( float alpha, float beta );

double
SafeNorm
( double alpha, double beta );

float
SafeNorm
( float alpha, float beta, float gamma );

double
SafeNorm
( double alpha, double beta, double gamma );

void
SVD
( char UColumns, char VColumns, int m, int n, float* A, int lda,
  float* SigmaDiag, float* U, int ldu, float* VT, int ldvt );

void
SVD
( char UColumns, char VColumns, int m, int n, double* A, int lda,
  double* SigmaDiag, double* U, int ldu, double* VT, int ldvt );

#ifndef WITHOUT_COMPLEX
void
SVD
( char UColumns, char VColumns, int m, int n, scomplex* A, int lda,
  float* SigmaDiag, scomplex* U, int ldu, scomplex* VT, int ldvt );

void
SVD
( char UColumns, char VColumns, int m, int n, dcomplex* A, int lda,
  double* SigmaDiag, dcomplex* U, int ldu, dcomplex* VT, int ldvt );
#endif

void
Tridiag
( char uplo, int n, float* A, int lda );

void
Tridiag
( char uplo, int n, double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Tridiag
( char uplo, int n, scomplex* A, int lda, scomplex* t );

void
Tridiag
( char uplo, int n, dcomplex* A, int lda, dcomplex* t );
#endif

void
Trinv
( char uplo, char diag, int n, const float* A, int lda );

void
Trinv
( char uplo, char diag, int n, const double* A, int lda );

#ifndef WITHOUT_COMPLEX
void
Trinv
( char uplo, char diag, int n, const scomplex* A, int lda );

void
Trinv
( char uplo, char diag, int n, const dcomplex* A, int lda );
#endif

} // lapack
} // wrappers
} // elemental

extern "C" {

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

//----------------------------------------------------------------------------//
// Implementations begin here                                                 //
//----------------------------------------------------------------------------//

inline void
elemental::wrappers::lapack::Chol
( char uplo, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Chol");
#endif
    int info;
    LAPACK(spotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "spotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Chol
( char uplo, int n, const double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Chol");
#endif
    int info;
    LAPACK(dpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::Chol
( char uplo, int n, const scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Chol");
#endif
    int info;
    LAPACK(cpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Chol
( char uplo, int n, const dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Chol");
#endif
    int info;
    LAPACK(zpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::lapack::Hegst
( int itype, char uplo, int n,
  float* A, int lda, const float* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Hegst");
#endif
    int info;
    LAPACK(ssygst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ssygst returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Hegst
( int itype, char uplo, int n,
  double* A, int lda, const double* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Hegst");
#endif
    int info;
    LAPACK(dsygst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dsygst returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::Hegst
( int itype, char uplo, int n,
  scomplex* A, int lda, const scomplex* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Hegst");
#endif
    int info;
    LAPACK(chegst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "chegst returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Hegst
( int itype, char uplo, int n,
  dcomplex* A, int lda, const dcomplex* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Hegst");
#endif
    int info;
    LAPACK(zhegst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zhegst returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::lapack::LU
( int m, int n, float* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::LU");
#endif
    int info;
    LAPACK(sgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::LU
( int m, int n, double* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::LU");
#endif
    int info;
    LAPACK(dgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::LU
( int m, int n, scomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::LU");
#endif
    int info;
    LAPACK(cgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::LU
( int m, int n, dcomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::LU");
#endif
    int info;
    LAPACK(zgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::lapack::QR
( int m, int n, float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::QR");
#endif
    int info;
    int lwork;
    float workSize;

    std::vector<float> t(std::min(m,n));

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(sgeqrf)( &m, &n, A, &lda, &t[0], &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<float> work(lwork); 
    LAPACK(sgeqrf)( &m, &n, A, &lda, &t[0], &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::QR
( int m, int n, double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::QR");
#endif
    int info;
    int lwork;
    double workSize;

    std::vector<double> t(std::min(m,n));

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(dgeqrf)( &m, &n, A, &lda, &t[0], &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<double> work(lwork); 
    LAPACK(dgeqrf)( &m, &n, A, &lda, &t[0], &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::QR
( int m, int n, scomplex* A, int lda, scomplex* t )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::QR");
#endif
    int info;
    int lwork;
    scomplex workSize;

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(cgeqrf)( &m, &n, A, &lda, t, &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<scomplex> work(lwork); 
    LAPACK(cgeqrf)( &m, &n, A, &lda, t, &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::QR
( int m, int n, dcomplex* A, int lda, dcomplex* t )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::QR");
#endif
    int info;
    int lwork;
    dcomplex workSize;

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(zgeqrf)( &m, &n, A, &lda, t, &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<dcomplex> work(lwork); 
    LAPACK(zgeqrf)( &m, &n, A, &lda, t, &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}
#endif

inline float
elemental::wrappers::lapack::SafeNorm
( float alpha, float beta )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SafeNorm");
#endif
    float gamma = LAPACK(slapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

inline double
elemental::wrappers::lapack::SafeNorm
( double alpha, double beta )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SafeNorm");
#endif
    double gamma = LAPACK(dlapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

inline float
elemental::wrappers::lapack::SafeNorm
( float alpha, float beta, float gamma )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SafeNorm");
#endif
    float delta = LAPACK(slapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

inline double
elemental::wrappers::lapack::SafeNorm
( double alpha, double beta, double gamma )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SafeNorm");
#endif
    double delta = LAPACK(dlapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

inline void
elemental::wrappers::lapack::SVD
( char UColumns, char VColumns, int m, int n, float* A, int lda,
  float* SigmaDiag, float* U, int ldu, float* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SVD");
#endif
    int info;
    int lwork;
    float workSize;

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(sgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &workSize, &lwork, &info );

    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<float> work(lwork);
    LAPACK(sgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::SVD
( char UColumns, char VColumns, int m, int n, double* A, int lda,
  double* SigmaDiag, double* U, int ldu, double* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SVD");
#endif
    int info;
    int lwork;
    double workSize;

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(dgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &workSize, &lwork, &info );

    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<double> work(lwork);
    LAPACK(dgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::SVD
( char UColumns, char VColumns, int m, int n, scomplex* A, int lda,
  float* SigmaDiag, scomplex* U, int ldu, scomplex* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SVD");
#endif
    int info;
    int lwork;
    scomplex workSize;
    std::vector<float> realWork(5*std::min(m,n));

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(cgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &workSize, &lwork, &realWork[0], &info );

    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<scomplex> work(lwork);
    LAPACK(cgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &work[0], &lwork, &realWork[0], &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::SVD
( char UColumns, char VColumns, int m, int n, dcomplex* A, int lda,
  double* SigmaDiag, dcomplex* U, int ldu, dcomplex* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::SVD");
#endif
    int info;
    int lwork;
    dcomplex workSize;
    std::vector<double> realWork(5*std::min(m,n));

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(zgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &workSize, &lwork, &realWork[0], &info );

    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<dcomplex> work(lwork);
    LAPACK(zgesvd)
    ( &UColumns, &VColumns, &m, &n, A, &lda, SigmaDiag, U, &ldu, VT, &ldvt,
      &work[0], &lwork, &realWork[0], &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::lapack::Tridiag
( char uplo, int n, float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    int lwork;
    float workSize;

    std::vector<float> d(n);
    std::vector<float> e(n-1);
    std::vector<float> t(n-1);

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(ssytrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], &t[0], &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<float> work(lwork); 
    LAPACK(ssytrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], &t[0], &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ssytrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Tridiag
( char uplo, int n, double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    int lwork;
    double workSize;

    std::vector<double> d(n);
    std::vector<double> e(n-1);
    std::vector<double> t(n-1);

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(dsytrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], &t[0], &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<double> work(lwork); 
    LAPACK(dsytrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], &t[0], &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dsytrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::Tridiag
( char uplo, int n, scomplex* A, int lda, scomplex* t )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    int lwork;
    scomplex workSize;

    std::vector<float> d(n);
    std::vector<float> e(n-1);

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(chetrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], t, &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<scomplex> work(lwork); 
    LAPACK(chetrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], t, &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "chetrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Tridiag
( char uplo, int n, dcomplex* A, int lda, dcomplex* t )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    int lwork;
    dcomplex workSize;

    std::vector<double> d(n);
    std::vector<double> e(n-1);

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(zhetrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], t, &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<dcomplex> work(lwork); 
    LAPACK(zhetrd)
    ( &uplo, &n, A, &lda, &d[0], &e[0], t, &work[0], &lwork, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zhetrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::lapack::Trinv
( char uplo, char diag, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Trinv");
#endif
    int info;
    LAPACK(strtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "strtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Trinv
( char uplo, char diag, int n, const double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Trinv");
#endif
    int info;
    LAPACK(dtrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dtrtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::Trinv
( char uplo, char diag, int n, const scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Trinv");
#endif
    int info;
    LAPACK(ctrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ctrtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Trinv
( char uplo, char diag, int n, const dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Trinv");
#endif
    int info;
    LAPACK(ztrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ztrtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}
#endif

#endif /* ELEMENTAL_WRAPPERS_LAPACK_HPP */

