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

#include "Elemental/Environment.hpp"

#ifdef FUNDERSCORE
#define C2F(name) name ## _
#else
#define C2F(name) name
#endif

namespace Elemental {
namespace wrappers {
namespace LAPACK {

void
Chol
( const char uplo, const int n, const float* A, const int lda );

void
Chol
( const char uplo, const int n, const double* A, const int lda );
 
#ifndef WITHOUT_COMPLEX
void
Chol
( const char uplo, const int n, const scomplex* A, const int lda );

void
Chol
( const char uplo, const int n, const dcomplex* A, const int lda );
#endif

void
LU
( const int m, const int n, float* A, const int lda, int* p );

void
LU
( const int m, const int n, double* A, const int lda, int* p );

#ifndef WITHOUT_COMPLEX
void
LU
( const int m, const int n, scomplex* A, const int lda, int* p );

void
LU
( const int m, const int n, dcomplex* A, const int lda, int* p );
#endif

float
SafeNorm
( const float alpha, const float beta );

double
SafeNorm
( const double alpha, const double beta );

float
SafeNorm
( const float alpha, const float beta, const float gamma );

double
SafeNorm
( const double alpha, const double beta, const double gamma );

void
Tridiag
( const char uplo,
  const int n, float* A, const int lda, 
  float* d, float* e, float* tau );

void
Tridiag
( const char uplo,
  const int n, double* A, const int lda,
  double* d, double* e, double* tau );

void
Trinv
( const char uplo, const char diag,
  const int n, const float* A, const int lda );

void
Trinv
( const char uplo, const char diag,
  const int n, const double* A, const int lda );

#ifndef WITHOUT_COMPLEX
void
Trinv
( const char uplo, const char diag,
  const int n, const scomplex* A, const int lda );

void
Trinv
( const char uplo, const char diag,
  const int n, const dcomplex* A, const int lda );
#endif

} // LAPACK 
} // wrappers
} // Elemental

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
  Elemental::scomplex* A, const int* lda, int* p, int* info );

void C2F(zgetrf)
( const int* m, const int* n, 
  Elemental::dcomplex* A, const int* lda, int* p, int* info );
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
( const char* uplo, const int* n, const Elemental::scomplex* A, 
  const int* lda, int* info );
    
void C2F(zpotrf)
( const char* uplo, const int* n, const Elemental::dcomplex* A, 
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
  const int* n, const Elemental::scomplex* A, const int* lda, int* info );
    
void C2F(ztrtri)
( const char* uplo, const char* diag,
  const int* n, const Elemental::dcomplex* A, const int* lda, int* info );
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

// Implementations begin here

inline void
Elemental::wrappers::LAPACK::Chol
( const char uplo, const int n, const float* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Chol");
#endif
    int info;
    C2F(spotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "spotrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::Chol
( const char uplo, const int n, const double* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Chol");
#endif
    int info;
    C2F(dpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dpotrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::LAPACK::Chol
( const char uplo, const int n, const scomplex* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Chol");
#endif
    int info;
    C2F(cpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cpotrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::Chol
( const char uplo, const int n, const dcomplex* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Chol");
#endif
    int info;
    C2F(zpotrf)( &uplo, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zpotrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::LAPACK::LU
( const int m, const int n, float* A, const int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::LU");
#endif
    int info;
    C2F(sgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgetrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::LU
( const int m, const int n, double* A, const int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::LU");
#endif
    int info;
    C2F(dgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgetrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::LAPACK::LU
( const int m, const int n, scomplex* A, const int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::LU");
#endif
    int info;
    C2F(cgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgetrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::LU
( const int m, const int n, dcomplex* A, const int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::LU");
#endif
    int info;
    C2F(zgetrf)( &m, &n, A, &lda, p, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgetrf returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}
#endif

inline float
Elemental::wrappers::LAPACK::SafeNorm
( const float alpha, const float beta )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::SafeNorm");
#endif
    float gamma = C2F(slapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

inline double
Elemental::wrappers::LAPACK::SafeNorm
( const double alpha, const double beta )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::SafeNorm");
#endif
    double gamma = C2F(dlapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

inline float
Elemental::wrappers::LAPACK::SafeNorm
( const float alpha, const float beta, const float gamma )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::SafeNorm");
#endif
    float delta = C2F(slapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

inline double
Elemental::wrappers::LAPACK::SafeNorm
( const double alpha, const double beta, const double gamma )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::SafeNorm");
#endif
    double delta = C2F(dlapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

inline void
Elemental::wrappers::LAPACK::Tridiag
( const char uplo,
  const int n, float* A, const int lda, float* d, float* e, float* tau )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Tridiag");
#endif
    int info;
    C2F(ssytd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ssytd2 returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::Tridiag
( const char uplo,
  const int n, double* A, const int lda, double* d, double* e, double* tau )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Tridiag");
#endif
    int info;
    C2F(dsytd2)( &uplo, &n, A, &lda, d, e, tau, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dsytd2 returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::Trinv
( const char uplo, const char diag, const int n, const float* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Trinv");
#endif
    int info;
    C2F(strtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "strtri returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::Trinv
( const char uplo, const char diag, 
  const int n, const double* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Trinv");
#endif
    int info;
    C2F(dtrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dtrtri returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::LAPACK::Trinv
( const char uplo, const char diag,
  const int n, const scomplex* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Trinv");
#endif
    int info;
    C2F(ctrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ctrtri returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::LAPACK::Trinv
( const char uplo, const char diag, 
  const int n, const dcomplex* A, const int lda )
{
#ifndef RELEASE
    PushCallStack("wrappers::LAPACK::Trinv");
#endif
    int info;
    C2F(ztrtri)( &uplo, &diag, &n, A, &lda, &info );
#ifndef RELEASE
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ztrtri returned with info = " << info;
        const std::string s = msg.str();
        throw s.c_str();
    }
    PopCallStack();
#endif
}
#endif

#endif /* ELEMENTAL_WRAPPERS_LAPACK_HPP */

