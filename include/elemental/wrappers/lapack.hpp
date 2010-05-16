/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_WRAPPERS_LAPACK_HPP
#define ELEMENTAL_WRAPPERS_LAPACK_HPP 1

#include "elemental/environment.hpp"

#ifdef FUNDERSCORE
#define C2F(name) name ## _
#else
#define C2F(name) name
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
( char uplo, int n, float* A, int lda, float* d, float* e, float* tau );

void
Tridiag
( char uplo, int n, double* A, int lda, double* d, double* e, double* tau );

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

void C2F(sgetrf)
( const int* m, const int* n, 
  float* A, const int* lda, int* p, int* info );

void C2F(dgetrf)
( const int* m, const int* n, 
  double* A, const int* lda, int* p, int* info );

#ifndef WITHOUT_COMPLEX
void C2F(cgetrf)
( const int* m, const int* n, 
  elemental::scomplex* A, const int* lda, int* p, int* info );

void C2F(zgetrf)
( const int* m, const int* n, 
  elemental::dcomplex* A, const int* lda, int* p, int* info );
#endif

float C2F(slapy2)
( const float* alpha, const float* beta );

double C2F(dlapy2)
( const double* alpha, const double* beta );

float C2F(slapy3)
( const float* alpha, const float* beta, const float* gamma );

double C2F(dlapy3)
( const double* alpha, const double* beta, const double* gamma );

void C2F(spotrf)
( const char* uplo, const int* n, const float* A, const int* lda,
  int* info );

void C2F(dpotrf)
( const char* uplo, const int* n, const double* A, const int* lda,
  int* info );
    
#ifndef WITHOUT_COMPLEX
void C2F(cpotrf)
( const char* uplo, const int* n, const elemental::scomplex* A, 
  const int* lda, int* info );
    
void C2F(zpotrf)
( const char* uplo, const int* n, const elemental::dcomplex* A, 
  const int* lda, int* info );
#endif
    
void C2F(strtri)
( const char* uplo, const char* diag, 
  const int* n, const float* A, const int* lda, int* info );

void C2F(dtrtri)
( const char* uplo, const char* diag, 
  const int* n, const double* A, const int* lda, int* info );
    
#ifndef WITHOUT_COMPLEX
void C2F(ctrtri)
( const char* uplo, const char* diag,
  const int* n, const elemental::scomplex* A, const int* lda, int* info );
    
void C2F(ztrtri)
( const char* uplo, const char* diag,
  const int* n, const elemental::dcomplex* A, const int* lda, int* info );
#endif

void C2F(ssytd2)
( const char* uplo,
  const int* n, float* A, const int* lda,
  float* d, float* e, float* tau, int* info );

void C2F(dsytd2)
( const char* uplo,
  const int* n, double* A, const int* lda, 
  double* d, double* e, double* tau, int* info );

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
    C2F(spotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "spotrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(dpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dpotrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(cpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cpotrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(zpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zpotrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(sgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgetrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(dgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgetrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(cgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgetrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(zgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgetrf returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    float gamma = C2F(slapy2)( &alpha, &beta );
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
    double gamma = C2F(dlapy2)( &alpha, &beta );
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
    float delta = C2F(slapy3)( &alpha, &beta, &gamma );
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
    double delta = C2F(dlapy3)( &alpha, &beta, &gamma );
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
    C2F(ssytd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ssytd2 returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(dsytd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dsytd2 returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
elemental::wrappers::lapack::Trinv
( char uplo, char diag, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::lapack::Trinv");
#endif
    int info;
    C2F(strtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "strtri returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(dtrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dtrtri returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(ctrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ctrtri returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
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
    C2F(ztrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ztrtri returned with info = " << info;
        const std::string& s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}
#endif

#endif /* ELEMENTAL_WRAPPERS_LAPACK_HPP */

