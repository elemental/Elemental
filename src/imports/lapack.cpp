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

// Machine constants
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

void
elemental::lapack::Cholesky
( char uplo, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::Cholesky");
#endif
    int info;
    LAPACK(spotrf)( &uplo, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "spotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zpotrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

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
        throw std::logic_error( msg.str() );
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
        throw std::logic_error( msg.str() );
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
        throw std::logic_error( msg.str() );
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
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

void
elemental::lapack::LU
( int m, int n, float* A, int lda, int* p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
#endif
    int info;
    LAPACK(sgetrf)( &m, &n, A, &lda, p, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgetrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

void
elemental::lapack::LQ
( int m, int n, float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::LQ");
#endif
    int info;
    int lwork;
    float workSize;

    std::vector<float> t(std::min(m,n));

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(sgelqf)( &m, &n, A, &lda, &t[0], &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<float> work(lwork); 
    LAPACK(sgelqf)( &m, &n, A, &lda, &t[0], &work[0], &lwork, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgelqf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::LQ
( int m, int n, double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::LQ");
#endif
    int info;
    int lwork;
    double workSize;

    std::vector<double> t(std::min(m,n));

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(dgelqf)( &m, &n, A, &lda, &t[0], &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)workSize;
    std::vector<double> work(lwork); 
    LAPACK(dgelqf)( &m, &n, A, &lda, &t[0], &work[0], &lwork, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgelqf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::LQ
( int m, int n, scomplex* A, int lda, scomplex* t )
{
#ifndef RELEASE
    PushCallStack("lapack::LQ");
#endif
    int info;
    int lwork;
    scomplex workSize;

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(cgelqf)( &m, &n, A, &lda, t, &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<scomplex> work(lwork); 
    LAPACK(cgelqf)( &m, &n, A, &lda, t, &work[0], &lwork, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgelqf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::LQ
( int m, int n, dcomplex* A, int lda, dcomplex* t )
{
#ifndef RELEASE
    PushCallStack("lapack::LQ");
#endif
    int info;
    int lwork;
    dcomplex workSize;

    // Retrieve the optimal worksize
    lwork = -1;
    LAPACK(zgelqf)( &m, &n, A, &lda, t, &workSize, &lwork, &info );
    
    // Allocate the work buffer and make the actual call
    lwork = (int)std::real(workSize);
    std::vector<dcomplex> work(lwork); 
    LAPACK(zgelqf)( &m, &n, A, &lda, t, &work[0], &lwork, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgelqf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

void
elemental::lapack::QR
( int m, int n, float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::QR");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::QR
( int m, int n, double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::QR");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::QR
( int m, int n, scomplex* A, int lda, scomplex* t )
{
#ifndef RELEASE
    PushCallStack("lapack::QR");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::QR
( int m, int n, dcomplex* A, int lda, dcomplex* t )
{
#ifndef RELEASE
    PushCallStack("lapack::QR");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgeqrf returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

float
elemental::lapack::SafeNorm
( float alpha, float beta )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    float gamma = LAPACK(slapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

double
elemental::lapack::SafeNorm
( double alpha, double beta )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    double gamma = LAPACK(dlapy2)( &alpha, &beta );
#ifndef RELEASE
    PopCallStack();
#endif
    return gamma;
}

float
elemental::lapack::SafeNorm
( float alpha, float beta, float gamma )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    float delta = LAPACK(slapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

double
elemental::lapack::SafeNorm
( double alpha, double beta, double gamma )
{
#ifndef RELEASE
    PushCallStack("lapack::SafeNorm");
#endif
    double delta = LAPACK(dlapy3)( &alpha, &beta, &gamma );
#ifndef RELEASE
    PopCallStack();
#endif
    return delta;
}

void
elemental::lapack::SVD
( char UColumns, char VColumns, int m, int n, float* A, int lda,
  float* SigmaDiag, float* U, int ldu, float* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::SVD");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "sgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::SVD
( char UColumns, char VColumns, int m, int n, double* A, int lda,
  double* SigmaDiag, double* U, int ldu, double* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::SVD");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::SVD
( char UColumns, char VColumns, int m, int n, scomplex* A, int lda,
  float* SigmaDiag, scomplex* U, int ldu, scomplex* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::SVD");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "cgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::SVD
( char UColumns, char VColumns, int m, int n, dcomplex* A, int lda,
  double* SigmaDiag, dcomplex* U, int ldu, dcomplex* VT, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::SVD");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zgesvd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

void
elemental::lapack::HermitianTridiag
( char uplo, int n, float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianTridiag");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ssytrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::HermitianTridiag
( char uplo, int n, double* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianTridiag");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dsytrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
void
elemental::lapack::HermitianTridiag
( char uplo, int n, scomplex* A, int lda, scomplex* t )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianTridiag");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "chetrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::lapack::HermitianTridiag
( char uplo, int n, dcomplex* A, int lda, dcomplex* t )
{
#ifndef RELEASE
    PushCallStack("lapack::HermitianTridiag");
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "zhetrd returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

void
elemental::lapack::TriangularInverse
( char uplo, char diag, int n, const float* A, int lda )
{
#ifndef RELEASE
    PushCallStack("lapack::TriangularInverse");
#endif
    int info;
    LAPACK(strtri)( &uplo, &diag, &n, A, &lda, &info );
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "strtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "dtrtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ctrtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
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
    if( info != 0 )
    {
        std::ostringstream msg;
        msg << "ztrtri returned with info = " << info;
        throw std::logic_error( msg.str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

