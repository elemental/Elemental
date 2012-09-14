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
#include "elemental.hpp"

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

// Bidiagonal QR
void LAPACK(sbdsqr)
( const char* uplo, const int* n, const int* numColsVTrans, const int* numRowsU,
  const int* numColsC, float* d, float* e, float* VTrans, const int* ldVTrans,
  float* U, const int* ldU, float* C, const int* ldC, float* work, int* info );
void LAPACK(dbdsqr)
( const char* uplo, const int* n, const int* numColsVTrans, const int* numRowsU,
  const int* numColsC, double* d, double* e,
  double* VTrans, const int* ldVTrans, double* U, const int* ldU,
  double* C, const int* ldC, double* work, int* info );
void LAPACK(cbdsqr)
( const char* uplo, const int* n, const int* numColsVAdj, const int* numRowsU,
  const int* numColsC, float* d, float* e,
  elem::scomplex* VAdj, const int* ldVAdj, elem::scomplex* U, const int* ldU,
  elem::scomplex* C, const int* ldC, float* work, int* info );
void LAPACK(zbdsqr)
( const char* uplo, const int* n, const int* numColsVAdj, const int* numRowsU,
  const int* numColsC, double* d, double* e,
  elem::dcomplex* VAdj, const int* ldVAdj, elem::dcomplex* U, const int* ldU,
  elem::dcomplex* C, const int* ldC, double* work, int* info );

// Divide and Conquer SVD
void LAPACK(sgesdd)
( const char* jobz, const int* m, const int* n, float* A, const int* lda,
  float* s, float* U, const int* ldu, float* VTrans, const int* ldvt,
  float* work, const int* lwork, int* iwork, int* info );
void LAPACK(dgesdd)
( const char* jobz, const int* m, const int* n, double* A, const int* lda,
  double* s, double* U, const int* ldu, double* VTrans, const int* ldvt,
  double* work, const int* lwork, int* iwork, int* info );
void LAPACK(cgesdd)
( const char* jobz, const int* m, const int* n,
  elem::scomplex* A, const int* lda, float* s,
  elem::scomplex* U, const int* ldu, elem::scomplex* VTrans, const int* ldvt,
  elem::scomplex* work, const int* lwork, float* rwork,
  int* iwork, int* info );
void LAPACK(zgesdd)
( const char* jobz, const int* m, const int* n,
  elem::dcomplex* A, const int* lda, double* s,
  elem::dcomplex* U, const int* ldu, elem::dcomplex* VAdj, const int* ldva,
  elem::dcomplex* work, const int* lwork, double* rwork,
  int* iwork, int* info );

// QR-algorithm SVD
void LAPACK(sgesvd)
( const char* jobu, const char* jobvt, const int* m, const int* n,
  float* A, const int* lda,
  float* s, float* U, const int* ldu, float* VTrans, const int* ldvt,
  float* work, const int* lwork, int* info );
void LAPACK(dgesvd)
( const char* jobu, const char* jobvt, const int* m, const int* n,
  double* A, const int* lda,
  double* s, double* U, const int* ldu, double* VTrans, const int* ldvt,
  double* work, const int* lwork, int* info );
void LAPACK(cgesvd)
( const char* jobu, const char* jobva, const int* m, const int* n,
  elem::scomplex* A, const int* lda, float* s,
  elem::scomplex* U, const int* ldu, elem::scomplex* VTrans, const int* ldvt,
  elem::scomplex* work, const int* lwork, float* rwork, int* info );
void LAPACK(zgesvd)
( const char* jobu, const char* jobva, const int* m, const int* n,
  elem::dcomplex* A, const int* lda, double* s,
  elem::dcomplex* U, const int* ldu, elem::dcomplex* VAdj, const int* ldva,
  elem::dcomplex* work, const int* lwork, double* rwork, int* info );

// Hessenberg QR algorithm
void LAPACK(shseqr)
( const char* job, const char* compz, const int* n, 
  const int* ilo, const int* ihi, float* H, const int* ldh, 
  float* wr, float* wi, float* Z, const int* ldz, 
  float* work, const int* lwork, int* info );
void LAPACK(dhseqr)
( const char* job, const char* compz, const int* n, 
  const int* ilo, const int* ihi, double* H, const int* ldh, 
  double* wr, double* wi, double* Z, const int* ldz, 
  double* work, const int* lwork, int* info );
void LAPACK(chseqr)
( const char* job, const char* compz, const int* n,
  const int* ilo, const int* ihi, elem::scomplex* H, const int* ldh,
  elem::scomplex* w, elem::scomplex* Z, const int* ldz,
  elem::scomplex* work, const int* lwork, int* info );
void LAPACK(zhseqr)
( const char* job, const char* compz, const int* n,
  const int* ilo, const int* ihi, elem::dcomplex* H, const int* ldh,
  elem::dcomplex* w, elem::dcomplex* Z, const int* ldz,
  elem::dcomplex* work, const int* lwork, int* info );

} // extern "C"

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

void Cholesky( char uplo, int n, const double* A, int lda )
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

void Cholesky( char uplo, int n, const scomplex* A, int lda )
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

void Cholesky( char uplo, int n, const dcomplex* A, int lda )
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

//
// LU factorization
//

void LU( int m, int n, float* A, int lda, int* p )
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

void LU( int m, int n, double* A, int lda, int* p )
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

void LU( int m, int n, scomplex* A, int lda, int* p )
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

void LU( int m, int n, dcomplex* A, int lda, int* p )
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

//
// Reduced a well-conditioned Hermitian generalized definite EVP to 
// standard form
//

void Hegst
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

void Hegst
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

void Hegst
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

void Hegst
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

//
// Triangular inversion
//

void TriangularInverse( char uplo, char diag, int n, const float* A, int lda )
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

void TriangularInverse( char uplo, char diag, int n, const double* A, int lda )
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

void TriangularInverse
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

void TriangularInverse
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

//
// Bidiagonal QR algorithm for SVD
//

void BidiagQRAlg
( char uplo, int n, int numColsVTrans, int numRowsU,
  float* d, float* e, float* VTrans, int ldVTrans, float* U, int ldU )
{
#ifndef RELEASE
    PushCallStack("lapack::BidiagQRAlg");
#endif
    if( n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    int info;
    float* C=0;
    const int numColsC=0, ldC=1;
    std::vector<float> work( 4*n );
    LAPACK(sbdsqr)
    ( &uplo, &n, &numColsVTrans, &numRowsU, &numColsC, d, e, VTrans, &ldVTrans,
      U, &ldU, C, &ldC, &work[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        std::ostringstream msg;
        msg << "sbdsqr had " << info << " elements of e not converge";
        throw std::runtime_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void BidiagQRAlg
( char uplo, int n, int numColsVTrans, int numRowsU, 
  double* d, double* e, double* VTrans, int ldVTrans, double* U, int ldU )
{
#ifndef RELEASE
    PushCallStack("lapack::BidiagQRAlg");
#endif
    if( n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    int info;
    double* C=0;
    const int numColsC=0, ldC=1;
    std::vector<double> work( 4*n );
    LAPACK(dbdsqr)
    ( &uplo, &n, &numColsVTrans, &numRowsU, &numColsC, d, e, VTrans, &ldVTrans,
      U, &ldU, C, &ldC, &work[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        std::ostringstream msg;
        msg << "dbdsqr had " << info << " elements of e not converge";
        throw std::runtime_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void BidiagQRAlg
( char uplo, int n, int numColsVAdj, int numRowsU, 
  float* d, float* e, scomplex* VAdj, int ldVAdj, scomplex* U, int ldU )
{
#ifndef RELEASE
    PushCallStack("lapack::BidiagQRAlg");
#endif
    if( n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    int info;
    scomplex* C=0;
    const int numColsC=0, ldC=1;
    std::vector<float> work( 4*n );
    LAPACK(cbdsqr)
    ( &uplo, &n, &numColsVAdj, &numRowsU, &numColsC, d, e, VAdj, &ldVAdj,
      U, &ldU, C, &ldC, &work[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        std::ostringstream msg;
        msg << "cbdsqr had " << info << " elements of e not converge";
        throw std::runtime_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void BidiagQRAlg
( char uplo, int n, int numColsVAdj, int numRowsU, 
  double* d, double* e, dcomplex* VAdj, int ldVAdj, dcomplex* U, int ldU )
{
#ifndef RELEASE
    PushCallStack("lapack::BidiagQRAlg");
#endif
    if( n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    int info;
    dcomplex* C=0;
    const int numColsC=0, ldC=1;
    std::vector<double> work( 4*n );
    LAPACK(zbdsqr)
    ( &uplo, &n, &numColsVAdj, &numRowsU, &numColsC, d, e, VAdj, &ldVAdj,
      U, &ldU, C, &ldC, &work[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        std::ostringstream msg;
        msg << "zbdsqr had " << info << " elements of e not converge";
        throw std::runtime_error( msg.str().c_str() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Divide and Conquer SVD
//

void DivideAndConquerSVD
( int m, int n, float* A, int lda, 
  float* s, float* U, int ldu, float* VTrans, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::DivideAndConquerSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobz='S';
    int lwork=-1, info;
    float dummyWork;
    const int k = std::min(m,n);
    std::vector<int> iwork(8*k);

    LAPACK(sgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, &dummyWork, &lwork,
      &iwork[0], &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    LAPACK(sgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, &work[0], &lwork,
      &iwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("sgesdd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void DivideAndConquerSVD
( int m, int n, double* A, int lda, 
  double* s, double* U, int ldu, double* VTrans, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::DivideAndConquerSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobz='S';
    int lwork=-1, info;
    double dummyWork;
    const int k = std::min(m,n);
    std::vector<int> iwork(8*k);

    LAPACK(dgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, &dummyWork, &lwork,
      &iwork[0], &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    LAPACK(dgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, &work[0], &lwork,
      &iwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("dgesdd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void DivideAndConquerSVD
( int m, int n, scomplex* A, int lda, 
  float* s, scomplex* U, int ldu, scomplex* VAdj, int ldva )
{
#ifndef RELEASE
    PushCallStack("lapack::DivideAndConquerSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobz='S';
    int lwork=-1, info;
    const int k = std::min(m,n);
    const int K = std::max(m,n);
    const int lrwork = k*std::max(5*k+7,2*K+2*k+1);
    std::vector<float> rwork(lrwork);
    std::vector<int> iwork(8*k);

    scomplex dummyWork;
    LAPACK(cgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, &dummyWork, &lwork,
      &rwork[0], &iwork[0], &info );

    lwork = dummyWork.real;
    std::vector<scomplex> work(lwork);
    LAPACK(cgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, &work[0], &lwork,
      &rwork[0], &iwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("cgesdd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void DivideAndConquerSVD
( int m, int n, dcomplex* A, int lda, 
  double* s, dcomplex* U, int ldu, dcomplex* VAdj, int ldva )
{
#ifndef RELEASE
    PushCallStack("lapack::DivideAndConquerSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobz='S';
    int lwork=-1, info;
    dcomplex dummyWork;
    const int k = std::min(m,n);
    const int K = std::max(m,n);
    const int lrwork = k*std::max(5*k+7,2*K+2*k+1);
    std::vector<double> rwork(lrwork);
    std::vector<int> iwork(8*k);

    LAPACK(zgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, &dummyWork, &lwork,
      &rwork[0], &iwork[0], &info );

    lwork = dummyWork.real;
    std::vector<dcomplex> work(lwork);
    LAPACK(zgesdd)
    ( &jobz, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, &work[0], &lwork,
      &rwork[0], &iwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("zgesdd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// QR-algorithm SVD
//

void QRSVD
( int m, int n, float* A, int lda, 
  float* s, float* U, int ldu, float* VTrans, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::QRSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='S', jobvt='S';
    int lwork=-1, info;
    float dummyWork;

    LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("sgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void QRSVD
( int m, int n, double* A, int lda, 
  double* s, double* U, int ldu, double* VTrans, int ldvt )
{
#ifndef RELEASE
    PushCallStack("lapack::QRSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='S', jobvt='S';
    int lwork=-1, info;
    double dummyWork;

    LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VTrans, &ldvt, 
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("dgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void QRSVD
( int m, int n, scomplex* A, int lda, 
  float* s, scomplex* U, int ldu, scomplex* VAdj, int ldva )
{
#ifndef RELEASE
    PushCallStack("lapack::QRSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='S', jobva='S';
    int lwork=-1, info;
    const int k = std::min(m,n);
    std::vector<float> rwork(5*k);

    scomplex dummyWork;
    LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      &dummyWork, &lwork, &rwork[0], &info );

    lwork = dummyWork.real;
    std::vector<scomplex> work(lwork);
    LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      &work[0], &lwork, &rwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("cgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void QRSVD
( int m, int n, dcomplex* A, int lda, 
  double* s, dcomplex* U, int ldu, dcomplex* VAdj, int ldva )
{
#ifndef RELEASE
    PushCallStack("lapack::QRSVD");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='S', jobva='S';
    int lwork=-1, info;
    dcomplex dummyWork;
    const int k = std::min(m,n);
    std::vector<double> rwork(5*k);

    LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      &dummyWork, &lwork, &rwork[0], &info );

    lwork = dummyWork.real;
    std::vector<dcomplex> work(lwork);
    LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, U, &ldu, VAdj, &ldva, 
      &work[0], &lwork, &rwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("zgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// Compute singular values with QR algorithm
//

void SingularValues( int m, int n, float* A, int lda, float* s )
{
#ifndef RELEASE
    PushCallStack("lapack::SingularValues");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='N', jobvt='N';
    int fakeLDim=1, lwork=-1, info;
    float dummyWork;

    LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    LAPACK(sgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("sgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void SingularValues( int m, int n, double* A, int lda, double* s )
{
#ifndef RELEASE
    PushCallStack("lapack::SingularValues");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='N', jobvt='N';
    int fakeLDim=1, lwork=-1, info;
    double dummyWork;

    LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    LAPACK(dgesvd)
    ( &jobu, &jobvt, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("dgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void SingularValues( int m, int n, scomplex* A, int lda, float* s )
{
#ifndef RELEASE
    PushCallStack("lapack::SingularValues");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='N', jobva='N';
    int fakeLDim=1, lwork=-1, info;
    scomplex dummyWork;
    const int k = std::min(m,n);
    std::vector<float> rwork(5*k);

    LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, &rwork[0], &info );

    lwork = dummyWork.real;
    std::vector<scomplex> work(lwork);
    LAPACK(cgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &work[0], &lwork, &rwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("cgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void SingularValues( int m, int n, dcomplex* A, int lda, double* s )
{
#ifndef RELEASE
    PushCallStack("lapack::SingularValues");
#endif
    if( m==0 || n==0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char jobu='N', jobva='N';
    int fakeLDim=1, lwork=-1, info;
    dcomplex dummyWork;
    const int k = std::min(m,n);
    std::vector<double> rwork(5*k);

    LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &dummyWork, &lwork, &rwork[0], &info );

    lwork = dummyWork.real;
    std::vector<dcomplex> work(lwork);
    LAPACK(zgesvd)
    ( &jobu, &jobva, &m, &n, A, &lda, s, 0, &fakeLDim, 0, &fakeLDim, 
      &work[0], &lwork, &rwork[0], &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("zgesvd's updating process failed");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// 
// Compute the eigenvalues of an upper Hessenberg matrix
//

void HessenbergEig( int n, float* H, int ldh, scomplex* w )
{
#ifndef RELEASE
    PushCallStack("lapack::HessenbergEig");
#endif
    if( n == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    float dummyWork;
    std::vector<float> wr( n ), wi( n );
    LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, &wr[0], &wi[0], 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<float> work(lwork);
    LAPACK(shseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, &wr[0], &wi[0], 0, &fakeLDim, 
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("shseqr's failed to compute all eigenvalues");
    }

    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<float>(wr[i],wi[i]);
#ifndef RELEASE
    PopCallStack();
#endif
}

void HessenbergEig( int n, double* H, int ldh, dcomplex* w )
{
#ifndef RELEASE
    PushCallStack("lapack::HessenbergEig");
#endif
    if( n == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    double dummyWork;
    std::vector<double> wr( n ), wi( n );
    LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, &wr[0], &wi[0], 0, &fakeLDim,
      &dummyWork, &lwork, &info );

    lwork = dummyWork;
    std::vector<double> work(lwork);
    LAPACK(dhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, &wr[0], &wi[0], 0, &fakeLDim,
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("dhseqr's failed to compute all eigenvalues");
    }
    
    for( int i=0; i<n; ++i )
        w[i] = elem::Complex<double>(wr[i],wi[i]);
#ifndef RELEASE
    PopCallStack();
#endif
}

void HessenbergEig( int n, scomplex* H, int ldh, scomplex* w )
{
#ifndef RELEASE
    PushCallStack("lapack::HessenbergEig");
#endif
    if( n == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    scomplex dummyWork;
    LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork.real;
    std::vector<scomplex> work(lwork);
    LAPACK(chseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("chseqr's failed to compute all eigenvalues");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void HessenbergEig( int n, dcomplex* H, int ldh, dcomplex* w )
{
#ifndef RELEASE
    PushCallStack("lapack::HessenbergEig");
#endif
    if( n == 0 )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    const char job='E', compz='N';
    int ilo=1, ihi=n;
    int fakeLDim=1, lwork=-1, info;
    dcomplex dummyWork;
    LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      &dummyWork, &lwork, &info );

    lwork = dummyWork.real;
    std::vector<dcomplex> work(lwork);
    LAPACK(zhseqr)
    ( &job, &compz, &n, &ilo, &ihi, H, &ldh, w, 0, &fakeLDim, 
      &work[0], &lwork, &info );
    if( info < 0 )
    {
        std::ostringstream msg;
        msg << "Argument " << -info << " had illegal value";
        throw std::logic_error( msg.str().c_str() );
    }
    else if( info > 0 )
    {
        throw std::runtime_error("zhseqr's failed to compute all eigenvalues");
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace lapack
} // namespace elem
