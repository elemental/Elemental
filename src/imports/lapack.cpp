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
#include "elemental/core/environment.hpp"

//
// Machine constants
//

template<> float 
elemental::lapack::MachineEpsilon<float>()
{
    const char cmach = 'E';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachineEpsilon<double>()
{
    const char cmach = 'E';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
elemental::lapack::MachineSafeMin<float>()
{
    const char cmach = 'S';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachineSafeMin<double>()
{
    const char cmach = 'S';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
elemental::lapack::MachineBase<float>()
{
    const char cmach = 'B';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachineBase<double>()
{
    const char cmach = 'B';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
elemental::lapack::MachinePrecision<float>()
{
    const char cmach = 'P';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachinePrecision<double>()
{
    const char cmach = 'P';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
elemental::lapack::MachineUnderflowExponent<float>()
{
    const char cmach = 'M';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachineUnderflowExponent<double>()
{
    const char cmach = 'M';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
elemental::lapack::MachineUnderflowThreshold<float>()
{
    const char cmach = 'U';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachineUnderflowThreshold<double>()
{
    const char cmach = 'U';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
elemental::lapack::MachineOverflowExponent<float>()
{
    const char cmach = 'L';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachineOverflowExponent<double>()
{
    const char cmach = 'L';
    return LAPACK(dlamch)( &cmach );
}

template<> float 
elemental::lapack::MachineOverflowThreshold<float>()
{
    const char cmach = 'O';
    return LAPACK(slamch)( &cmach );
}

template<> double 
elemental::lapack::MachineOverflowThreshold<double>()
{
    const char cmach = 'O';
    return LAPACK(dlamch)( &cmach );
}

//
// Safely compute norms
//

float
elemental::lapack::SafeNorm
( float alpha, float beta )
{ return LAPACK(slapy2)( &alpha, &beta ); }

double
elemental::lapack::SafeNorm
( double alpha, double beta )
{ return LAPACK(dlapy2)( &alpha, &beta ); }

float
elemental::lapack::SafeNorm
( float alpha, float beta, float gamma )
{ return LAPACK(slapy3)( &alpha, &beta, &gamma ); }

double
elemental::lapack::SafeNorm
( double alpha, double beta, double gamma )
{ return LAPACK(dlapy3)( &alpha, &beta, &gamma ); }

//
// Safely compute Givens rotations (using Demmel and Kahan's algorithm)
//

void
elemental::lapack::ComputeGivens
( float phi, float gamma,
  float* c, float* s, float* rho )
{ LAPACK(slartg)( &phi, &gamma, c, s, rho ); }

void
elemental::lapack::ComputeGivens
( double phi, double gamma,
  double* c, double* s, double* rho )
{ LAPACK(dlartg)( &phi, &gamma, c, s, rho ); }

void
elemental::lapack::ComputeGivens
( scomplex phi, scomplex gamma,
  float* c, scomplex* s, scomplex* rho )
{ LAPACK(clartg)( &phi, &gamma, c, s, rho ); }

void
elemental::lapack::ComputeGivens
( dcomplex phi, dcomplex gamma,
  double* c, dcomplex* s, dcomplex* rho )
{ LAPACK(zlartg)( &phi, &gamma, c, s, rho ); }

//
// Cholesky factorization
//

void
elemental::lapack::Cholesky
( char uplo, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Cholesky");
#endif
    int info;
    LAPACK(spotrf)( &uplo, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "spotrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw NonHPDMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::Cholesky
( char uplo, int n, const double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Cholesky");
#endif
    int info;
    LAPACK(dpotrf)( &uplo, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "dpotrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw NonHPDMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::Cholesky
( char uplo, int n, const scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Cholesky");
#endif
    int info;
    LAPACK(cpotrf)( &uplo, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "cpotrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw NonHPDMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::Cholesky
( char uplo, int n, const dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Cholesky");
#endif
    int info;
    LAPACK(zpotrf)( &uplo, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "zpotrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw NonHPDMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

//
// LU factorization
//

void
elemental::lapack::LU
( int m, int n, float* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(sgetrf)( &m, &n, A, &lda, p, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "sgetrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::LU
( int m, int n, double* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(dgetrf)( &m, &n, A, &lda, p, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "dgetrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::LU
( int m, int n, scomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(cgetrf)( &m, &n, A, &lda, p, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "cgetrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::LU
( int m, int n, dcomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(zgetrf)( &m, &n, A, &lda, p, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "zgetrf returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

//
// Reduced a well-conditioned Hermitian generalized definite EVP to 
// standard form
//

void
elemental::lapack::Hegst
( int itype, char uplo, int n,
  float* A, int lda, const float* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("lapack::Hegst");
#endif
    int info;
    LAPACK(ssygst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ssygst returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::Hegst
( int itype, char uplo, int n,
  double* A, int lda, const double* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("lapack::Hegst");
#endif
    int info;
    LAPACK(dsygst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dsygst returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::Hegst
( int itype, char uplo, int n,
  scomplex* A, int lda, const scomplex* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("lapack::Hegst");
#endif
    int info;
    LAPACK(chegst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "chegst returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::Hegst
( int itype, char uplo, int n,
  dcomplex* A, int lda, const dcomplex* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("lapack::Hegst");
#endif
    int info;
    LAPACK(zhegst)( &itype, &uplo, &n, A, &lda, B, &ldb, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zhegst returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

//
// Triangular inversion
//

void
elemental::lapack::TriangularInverse
( char uplo, char diag, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::TriangularInverse");
#endif
    int info;
    LAPACK(strtri)( &uplo, &diag, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "strtri returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::TriangularInverse
( char uplo, char diag, int n, const double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::TriangularInverse");
#endif
    int info;
    LAPACK(dtrtri)( &uplo, &diag, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "dtrtri returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::TriangularInverse
( char uplo, char diag, int n, const scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::TriangularInverse");
#endif
    int info;
    LAPACK(ctrtri)( &uplo, &diag, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "ctrtri returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::TriangularInverse
( char uplo, char diag, int n, const dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::TriangularInverse");
#endif
    int info;
    LAPACK(ztrtri)( &uplo, &diag, &n, A, &lda, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "ztrtri returned with info = " << info;
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
        throw SingularMatrixException();
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

