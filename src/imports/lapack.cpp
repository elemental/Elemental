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
#include "elemental/core/environment.hpp"

namespace elem {
namespace lapack {

//
// Machine constants
//

template<>
float MachineEpsilon<float>()
{
    const char cmach = 'E';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachineEpsilon<double>()
{
    const char cmach = 'E';
    return LAPACK(dlamch)( &cmach );
}

template<> 
float MachineSafeMin<float>()
{
    const char cmach = 'S';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachineSafeMin<double>()
{
    const char cmach = 'S';
    return LAPACK(dlamch)( &cmach );
}

template<> 
float MachineBase<float>()
{
    const char cmach = 'B';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachineBase<double>()
{
    const char cmach = 'B';
    return LAPACK(dlamch)( &cmach );
}

template<>
float MachinePrecision<float>()
{
    const char cmach = 'P';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachinePrecision<double>()
{
    const char cmach = 'P';
    return LAPACK(dlamch)( &cmach );
}

template<> 
float MachineUnderflowExponent<float>()
{
    const char cmach = 'M';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachineUnderflowExponent<double>()
{
    const char cmach = 'M';
    return LAPACK(dlamch)( &cmach );
}

template<>
float MachineUnderflowThreshold<float>()
{
    const char cmach = 'U';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachineUnderflowThreshold<double>()
{
    const char cmach = 'U';
    return LAPACK(dlamch)( &cmach );
}

template<> 
float MachineOverflowExponent<float>()
{
    const char cmach = 'L';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachineOverflowExponent<double>()
{
    const char cmach = 'L';
    return LAPACK(dlamch)( &cmach );
}

template<> 
float MachineOverflowThreshold<float>()
{
    const char cmach = 'O';
    return LAPACK(slamch)( &cmach );
}

template<> 
double MachineOverflowThreshold<double>()
{
    const char cmach = 'O';
    return LAPACK(dlamch)( &cmach );
}

//
// Safely compute norms
//

float SafeNorm( float alpha, float beta )
{ return LAPACK(slapy2)( &alpha, &beta ); }

double SafeNorm( double alpha, double beta )
{ return LAPACK(dlapy2)( &alpha, &beta ); }

float SafeNorm( float alpha, float beta, float gamma )
{ return LAPACK(slapy3)( &alpha, &beta, &gamma ); }

double SafeNorm( double alpha, double beta, double gamma )
{ return LAPACK(dlapy3)( &alpha, &beta, &gamma ); }

//
// Safely compute Givens rotations (using Demmel and Kahan's algorithm)
//

void ComputeGivens
( float phi, float gamma,
  float* c, float* s, float* rho )
{ LAPACK(slartg)( &phi, &gamma, c, s, rho ); }

void ComputeGivens
( double phi, double gamma,
  double* c, double* s, double* rho )
{ LAPACK(dlartg)( &phi, &gamma, c, s, rho ); }

void ComputeGivens
( scomplex phi, scomplex gamma,
  float* c, scomplex* s, scomplex* rho )
{ LAPACK(clartg)( &phi, &gamma, c, s, rho ); }

void ComputeGivens
( dcomplex phi, dcomplex gamma,
  double* c, dcomplex* s, dcomplex* rho )
{ LAPACK(zlartg)( &phi, &gamma, c, s, rho ); }

//
// Cholesky factorization
//

void Cholesky( char uplo, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
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

void Cholesky( char uplo, int n, const double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
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

void Cholesky( char uplo, int n, const scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
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

void Cholesky( char uplo, int n, const dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
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

//
// LU factorization
//

void LU( int m, int n, float* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("LU");
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

void LU( int m, int n, double* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("LU");
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

void LU( int m, int n, scomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("LU");
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

void LU( int m, int n, dcomplex* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("LU");
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

//
// Reduced a well-conditioned Hermitian generalized definite EVP to 
// standard form
//

void Hegst
( int itype, char uplo, int n,
  float* A, int lda, const float* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("Hegst");
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

void Hegst
( int itype, char uplo, int n,
  double* A, int lda, const double* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("Hegst");
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

void Hegst
( int itype, char uplo, int n,
  scomplex* A, int lda, const scomplex* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("Hegst");
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

void Hegst
( int itype, char uplo, int n,
  dcomplex* A, int lda, const dcomplex* B, int ldb )
{
#ifndef RELEASE
    PushCallStack("Hegst");
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

//
// Triangular inversion
//

void TriangularInverse( char uplo, char diag, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
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

void TriangularInverse( char uplo, char diag, int n, const double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
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

void TriangularInverse
( char uplo, char diag, int n, const scomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
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

void TriangularInverse
( char uplo, char diag, int n, const dcomplex* A, int lda )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
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

} // namespace lapack
} // namespace elem
