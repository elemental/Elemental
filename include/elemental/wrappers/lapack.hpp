/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_WRAPPERS_LAPACK_HPP
#define ELEMENTAL_WRAPPERS_LAPACK_HPP 1

#include "elemental/environment.hpp"

#ifdef LAPACK_UNDERSCORE
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
Tridiag
( char uplo, int n, float* A, int lda, float* d, float* e, float* t );

void
Tridiag
( char uplo, int n, double* A, int lda, double* d, double* e, double* t );

#ifndef WITHOUT_COMPLEX
void
Tridiag
( char uplo, int n, scomplex* A, int lda, float* d, float* e, scomplex* t );

void
Tridiag
( char uplo, int n, dcomplex* A, int lda, double* d, double* e, dcomplex* t );
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
void LAPACK(ssytd2)
( const char* uplo,
  const int* n, float* A, const int* lda,
  float* d, float* e, float* tau, int* info );

void LAPACK(dsytd2)
( const char* uplo,
  const int* n, double* A, const int* lda, 
  double* d, double* e, double* tau, int* info );

#ifndef WITHOUT_COMPLEX
void LAPACK(chetd2)
( const char* uplo,
  const int* n, elemental::scomplex* A, const int* lda,
  float* d, float* e, elemental::scomplex* tau, int* info );

void LAPACK(zhetd2)
( const char* uplo,
  const int* n, elemental::dcomplex* A, const int* lda,
  double* d, double* e, elemental::dcomplex* tau, int* info );
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
elemental::wrappers::lapack::Tridiag
( char uplo, int n, float* A, int lda, float* d, float* e, float* tau )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    LAPACK(ssytd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ssytd2 returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Tridiag
( char uplo, int n, double* A, int lda, double* d, double* e, double* tau )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    LAPACK(dsytd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dsytd2 returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::lapack::Tridiag
( char uplo, int n, scomplex* A, int lda, float* d, float* e, scomplex* tau )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    LAPACK(chetd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "chetd2 returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Tridiag
( char uplo, int n, dcomplex* A, int lda, double* d, double* e, dcomplex* tau )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Tridiag");
#endif
    int info;
    LAPACK(zhetd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zhetd2 returned with info = " << info;
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

