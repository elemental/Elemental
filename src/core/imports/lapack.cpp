/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/lapack_like.hpp>

using El::FortranLogical;
using El::BlasInt;
using El::scomplex;
using El::dcomplex;

extern "C" {

// Copy matrices
void EL_LAPACK(slacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const float* A, const BlasInt* lda, float* B, const BlasInt* ldb );
void EL_LAPACK(dlacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const double* A, const BlasInt* lda, double* B, const BlasInt* ldb );
void EL_LAPACK(clacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const scomplex* A, const BlasInt* lda, scomplex* B, const BlasInt* ldb );
void EL_LAPACK(zlacpy)
( const char* uplo, const BlasInt* m, const BlasInt* n, 
  const dcomplex* A, const BlasInt* lda, dcomplex* B, const BlasInt* ldb );

// Symmetric tridiagonal eigensolvers (via MRRR)
void EL_LAPACK(sstevr)
( const char* job, const char* range, const BlasInt* n,
  float* d, float* e, const float* vl, const float* vu, 
  const BlasInt* il, const BlasInt* iu, const float* absTol, BlasInt * m, 
  float* w, float* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  float* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );
void EL_LAPACK(dstevr)
( const char* job, const char* range, const BlasInt* n,
  double* d, double* e, const double* vl, const double* vu, 
  const BlasInt* il, const BlasInt* iu, const double* absTol, BlasInt * m, 
  double* w, double* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  double* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );

// Hermitian eigensolvers (via MRRR)
void EL_LAPACK(ssyevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  float* A, const BlasInt* ldA, const float* vl, const float* vu, 
  const BlasInt* il, const BlasInt* iu, const float* absTol, BlasInt * m, 
  float* w, float* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  float* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );
void EL_LAPACK(dsyevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  double* A, const BlasInt* ldA, const double* vl, const double* vu, 
  const BlasInt* il, const BlasInt* iu, const double* absTol, BlasInt * m, 
  double* w, double* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  double* work, const BlasInt* workSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, 
  BlasInt* info );
void EL_LAPACK(cheevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  scomplex* A, const BlasInt* ldA, const float* vl, const float* vu, 
  const BlasInt* il, const BlasInt* iu, const float* absTol, BlasInt* m,
  float* w, scomplex* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  scomplex* work, const BlasInt* workSize, 
  float* rWork, const BlasInt* rWorkSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, BlasInt* info );
void EL_LAPACK(zheevr)
( const char* job, const char* range, const char* uplo, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA, const double* vl, const double* vu, 
  const BlasInt* il, const BlasInt* iu, const double* absTol, BlasInt* m,
  double* w, dcomplex* Z, const BlasInt* ldZ, BlasInt* isuppZ, 
  dcomplex* work, const BlasInt* workSize, 
  double* rWork, const BlasInt* rWorkSize, 
  BlasInt* iWork, const BlasInt* iWorkSize, BlasInt* info );

// Bidiagonal DQDS
void EL_LAPACK(slasq1)
( const BlasInt* n, float* d, float* e, float* work, BlasInt* info );
void EL_LAPACK(dlasq1)
( const BlasInt* n, double* d, double* e, double* work, BlasInt* info );

// Bidiagonal QR
void EL_LAPACK(sbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVT, const BlasInt* numRowsU,
  const BlasInt* numColsC, float* d, float* e, 
  float* VTrans, const BlasInt* ldVT,
  float* U, const BlasInt* ldU, float* C, const BlasInt* ldC, 
  float* work, BlasInt* info );
void EL_LAPACK(dbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVT, const BlasInt* numRowsU,
  const BlasInt* numColsC, double* d, double* e,
  double* VTrans, const BlasInt* ldVT, double* U, const BlasInt* ldU,
  double* C, const BlasInt* ldC, double* work, BlasInt* info );
void EL_LAPACK(cbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVH, const BlasInt* numRowsU,
  const BlasInt* numColsC, float* d, float* e,
  scomplex* VH, const BlasInt* ldVH, scomplex* U, const BlasInt* ldU,
  scomplex* C, const BlasInt* ldC, float* work, BlasInt* info );
void EL_LAPACK(zbdsqr)
( const char* uplo, const BlasInt* n, 
  const BlasInt* numColsVH, const BlasInt* numRowsU,
  const BlasInt* numColsC, double* d, double* e,
  dcomplex* VH, const BlasInt* ldVH, dcomplex* U, const BlasInt* ldU,
  dcomplex* C, const BlasInt* ldC, double* work, BlasInt* info );

// Divide and Conquer SVD
void EL_LAPACK(sgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n, 
  float* A, const BlasInt* ldA,
  float* s, float* U, const BlasInt* ldu, 
  float* VTrans, const BlasInt* ldvt,
  float* work, const BlasInt* workSize, BlasInt* iWork, BlasInt* info );
void EL_LAPACK(dgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n, 
  double* A, const BlasInt* ldA,
  double* s, double* U, const BlasInt* ldu, 
  double* VTrans, const BlasInt* ldvt,
  double* work, const BlasInt* workSize, BlasInt* iWork, BlasInt* info );
void EL_LAPACK(cgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n,
  scomplex* A, const BlasInt* ldA, float* s,
  scomplex* U, const BlasInt* ldu, scomplex* VTrans, const BlasInt* ldvt,
  scomplex* work, const BlasInt* workSize, float* rWork,
  BlasInt* iWork, BlasInt* info );
void EL_LAPACK(zgesdd)
( const char* jobz, const BlasInt* m, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA, double* s,
  dcomplex* U, const BlasInt* ldu, dcomplex* VH, const BlasInt* ldva,
  dcomplex* work, const BlasInt* workSize, double* rWork,
  BlasInt* iWork, BlasInt* info );

// QR-algorithm SVD [DQDS when no singular vectors desired]
void EL_LAPACK(sgesvd)
( const char* jobU, const char* jobVT, const BlasInt* m, const BlasInt* n,
  float* A, const BlasInt* ldA,
  float* s, float* U, const BlasInt* ldu, float* VTrans, const BlasInt* ldvt,
  float* work, const BlasInt* workSize, BlasInt* info );
void EL_LAPACK(dgesvd)
( const char* jobU, const char* jobVT, const BlasInt* m, const BlasInt* n,
  double* A, const BlasInt* ldA,
  double* s, double* U, const BlasInt* ldu, double* VTrans, const BlasInt* ldvt,
  double* work, const BlasInt* workSize, BlasInt* info );
void EL_LAPACK(cgesvd)
( const char* jobU, const char* jobVH, const BlasInt* m, const BlasInt* n,
  scomplex* A, const BlasInt* ldA, float* s,
  scomplex* U, const BlasInt* ldu, scomplex* VTrans, const BlasInt* ldvt,
  scomplex* work, const BlasInt* workSize, float* rWork, BlasInt* info );
void EL_LAPACK(zgesvd)
( const char* jobU, const char* jobVH, const BlasInt* m, const BlasInt* n,
  dcomplex* A, const BlasInt* ldA, double* s,
  dcomplex* U, const BlasInt* ldu, dcomplex* VH, const BlasInt* ldva,
  dcomplex* work, const BlasInt* workSize, double* rWork, BlasInt* info );

// Reduction to Hessenberg form
void EL_LAPACK(sgehrd)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  float* A, const BlasInt* ldA,
  float* tau,
  float* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(dgehrd)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  double* A, const BlasInt* ldA,
  double* tau,
  double* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(cgehrd)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  scomplex* A, const BlasInt* ldA,
  scomplex* tau,
  scomplex* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(zgehrd)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  dcomplex* A, const BlasInt* ldA,
  dcomplex* tau,
  dcomplex* work, const BlasInt* workSize,
  BlasInt* info );

// Generates a unitary matrix defined as the product of Householder reflectors
void EL_LAPACK(sorghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  float* A, const BlasInt* ldA,
  const float* tau,
  float* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(dorghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  double* A, const BlasInt* ldA,
  const double* tau,
  double* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(cunghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  scomplex* A, const BlasInt* ldA,
  const scomplex* tau,
  scomplex* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(zunghr)
( const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi, 
  dcomplex* A, const BlasInt* ldA,
  const dcomplex* tau,
  dcomplex* work, const BlasInt* workSize,
  BlasInt* info );

// Hessenberg QR algorithm (with AED)
void EL_LAPACK(shseqr)
( const char* job, const char* compZ, const BlasInt* n, 
  const BlasInt* ilo, const BlasInt* ihi,
  float* H, const BlasInt* ldH, 
  float* wr, float* wi,
  float* Z, const BlasInt* ldZ, 
  float* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(dhseqr)
( const char* job, const char* compZ, const BlasInt* n, 
  const BlasInt* ilo, const BlasInt* ihi,
  double* H, const BlasInt* ldH, 
  double* wr, double* wi,
  double* Z, const BlasInt* ldZ, 
  double* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(chseqr)
( const char* job, const char* compZ, const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi,
  scomplex* H, const BlasInt* ldH,
  scomplex* w,
  scomplex* Z, const BlasInt* ldZ,
  scomplex* work, const BlasInt* workSize,
  BlasInt* info );
void EL_LAPACK(zhseqr)
( const char* job, const char* compZ, const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi,
  dcomplex* H, const BlasInt* ldH,
  dcomplex* w,
  dcomplex* Z, const BlasInt* ldZ,
  dcomplex* work, const BlasInt* workSize,
  BlasInt* info );

// Hessenberg QR algorithm (without AED)
void EL_LAPACK(slahqr)
( const FortranLogical* wantT,
  const FortranLogical* wantZ,
  const BlasInt* n, 
  const BlasInt* ilo, const BlasInt* ihi,
  float* H, const BlasInt* ldH, 
  float* wr, float* wi,
  const BlasInt* iloZ, const BlasInt* ihiZ,
  float* Z, const BlasInt* ldZ, 
  BlasInt* info );
void EL_LAPACK(dlahqr)
( const FortranLogical* wantT,
  const FortranLogical* wantZ,
  const BlasInt* n, 
  const BlasInt* ilo, const BlasInt* ihi,
  double* H, const BlasInt* ldH, 
  double* wr, double* wi,
  const BlasInt* iloZ, const BlasInt* ihiZ,
  double* Z, const BlasInt* ldZ, 
  BlasInt* info );
void EL_LAPACK(clahqr)
( const FortranLogical* wantT,
  const FortranLogical* wantZ,
  const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi,
  scomplex* H, const BlasInt* ldH,
  scomplex* w,
  const BlasInt* iloZ, const BlasInt* ihiZ,
  scomplex* Z, const BlasInt* ldZ,
  BlasInt* info );
void EL_LAPACK(zlahqr)
( const FortranLogical* wantT,
  const FortranLogical* wantZ,
  const BlasInt* n,
  const BlasInt* ilo, const BlasInt* ihi,
  dcomplex* H, const BlasInt* ldH,
  dcomplex* w,
  const BlasInt* iloZ, const BlasInt* ihiZ,
  dcomplex* Z, const BlasInt* ldZ,
  BlasInt* info );

} // extern "C"

namespace El {
namespace lapack {

// Copy a matrix
// =============
template<typename T>
void Copy
( char uplo, BlasInt m, BlasInt n,
  const T* A, BlasInt lda,
        T* B, BlasInt ldb )
{
    if( std::toupper(uplo) == 'L' )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<m; ++i )
                B[i+j*ldb] = A[i+j*lda];
    }
    else if( std::toupper(uplo) == 'U' ) 
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                B[i+j*ldb] = A[i+j*lda];
    }
    else
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B[i+j*ldb] = A[i+j*lda];
    }
}
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Int* A, BlasInt lda,
        Int* B, BlasInt ldb );
#ifdef EL_HAVE_QD
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const DoubleDouble* A, BlasInt lda,
        DoubleDouble* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const QuadDouble* A, BlasInt lda,
        QuadDouble* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Complex<DoubleDouble>* A, BlasInt lda,
        Complex<DoubleDouble>* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Complex<QuadDouble>* A, BlasInt lda,
        Complex<QuadDouble>* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_QUAD
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Quad* A, BlasInt lda,
        Quad* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Complex<Quad>* A, BlasInt lda,
        Complex<Quad>* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_MPC
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const BigInt* A, BlasInt lda,
        BigInt* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const BigFloat* A, BlasInt lda,
        BigFloat* B, BlasInt ldb );
template void Copy
( char uplo, BlasInt m, BlasInt n, 
  const Complex<BigFloat>* A, BlasInt lda,
        Complex<BigFloat>* B, BlasInt ldb );
#endif

void Copy
( char uplo, BlasInt m, BlasInt n, 
  const float* A, BlasInt lda, float* B, BlasInt ldb )
{ EL_LAPACK(slacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const double* A, BlasInt lda, double* B, BlasInt ldb )
{ EL_LAPACK(dlacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb )
{ EL_LAPACK(clacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }
void Copy
( char uplo, BlasInt m, BlasInt n, 
  const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb )
{ EL_LAPACK(zlacpy)( &uplo, &m, &n, A, &lda, B, &ldb ); }

// Generate a Householder reflector
// ================================
// NOTE: 
// Since LAPACK chooses to use the identity matrix, rather than a single
// coordinate negation, in cases where the mass is already entirely in the
// first entry, and the identity matrix cannot be represented as a Householder
// reflector, Elemental does not ever directly call LAPACK's Householder
// routines. Otherwise, the logic of routines such as ApplyPackedReflectors
// would need to be (unnecessarily) complicated.
//
// Furthermore, LAPACK defines H = I - tau [1; v] [1; v]' such that
// adjoint(H) [chi; x] = [beta; 0], but Elemental instead defines
// H = I - tau [1; v] [1; v]' such that H [chi; x] = [beta; 0].
//
template<typename F>
F Reflector( BlasInt n, F& chi, F* x, BlasInt incx )
{
    DEBUG_CSE
    typedef Base<F> Real; 
    const Real zero(0);

    Real norm = blas::Nrm2( n-1, x, incx );
    F alpha = chi;
    if( norm == zero && ImagPart(alpha) == zero )
    {
        chi *= -1;
        return F(2); 
    }

    Real beta;
    if( RealPart(alpha) <= zero )
        beta = SafeNorm( alpha, norm );
    else
        beta = -SafeNorm( alpha, norm );

    // Rescale if the vector is too small
    const Real safeMin = limits::SafeMin<Real>();
    const Real epsilon = limits::Epsilon<Real>();
    const Real safeInv = safeMin/epsilon;
    Int count = 0;
    if( Abs(beta) < safeInv )
    {
        Real invOfSafeInv = Real(1)/safeInv;
        do
        {
            ++count;
            blas::Scal( n-1, invOfSafeInv, x, incx );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        norm = blas::Nrm2( n-1, x, incx );
        if( RealPart(alpha) <= 0 )
            beta = SafeNorm( alpha, norm );
        else
            beta = -SafeNorm( alpha, norm );
    }

    F tau = (beta-Conj(alpha)) / beta;
    blas::Scal( n-1, Real(1)/(alpha-beta), x, incx );

    // Undo the scaling
    for( Int j=0; j<count; ++j )
        beta *= safeInv;

    chi = beta;
    return tau;
}
template float Reflector
( BlasInt n, float& chi, float* x, BlasInt incx );
template Complex<float> Reflector
( BlasInt n, Complex<float>& chi, Complex<float>* x, BlasInt incx );
template double Reflector
( BlasInt n, double& chi, double* x, BlasInt incx );
template Complex<double> Reflector
( BlasInt n, Complex<double>& chi, Complex<double>* x, BlasInt incx );
#ifdef EL_HAVE_QD
template DoubleDouble Reflector
( BlasInt n, DoubleDouble& chi, DoubleDouble* x, BlasInt incx );
template QuadDouble Reflector
( BlasInt n, QuadDouble& chi, QuadDouble* x, BlasInt incx );
template Complex<DoubleDouble> Reflector
( BlasInt n,
  Complex<DoubleDouble>& chi,
  Complex<DoubleDouble>* x, BlasInt incx );
template Complex<QuadDouble> Reflector
( BlasInt n,
  Complex<QuadDouble>& chi,
  Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template Quad Reflector
( BlasInt n, Quad& chi, Quad* x, BlasInt incx );
template Complex<Quad> Reflector
( BlasInt n, Complex<Quad>& chi, Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template BigFloat Reflector
( BlasInt n, BigFloat& chi, BigFloat* x, BlasInt incx );
template Complex<BigFloat> Reflector
( BlasInt n, Complex<BigFloat>& chi, Complex<BigFloat>* x, BlasInt incx );
#endif

// Cf. LAPACK's ?LARF
//
// TODO(poulson):
// Provide wrappers for LAPACK for standard datatypes since, despite the
// differences in handling identity transformations, LAPACK's routine should
// always correctly apply a true Householder transformation.
//
template<typename F>
void ApplyReflector
( bool onLeft,
  BlasInt m,
  BlasInt n,
  const F* v, BlasInt vInc,
  const F& tau,
        F* C, BlasInt CLDim,
        F* work )
{
    // TODO(poulson): Add shortcuts for sparse v
    const F zero(0);
    if( onLeft )
    {
        for( BlasInt i=0; i<n; ++i )
            work[i] = zero;
        blas::Gemv( 'C', m, n, F(1), C, CLDim, v, vInc, F(0), work, 1 );
        blas::Ger( m, n, -tau, v, vInc, work, 1, C, CLDim );
    }
    else
    {
        for( BlasInt i=0; i<m; ++i )
            work[i] = zero;
        blas::Gemv( 'N', m, n, F(1), C, CLDim, v, vInc, F(0), work, 1 );
        blas::Ger( m, n, -tau, work, 1, v, vInc, C, CLDim );
    }
}

template<typename F>
void ApplyReflector
( bool onLeft,
  BlasInt m,
  BlasInt n,
  const F* v, BlasInt vInc,
  const F& tau,
        F* C, BlasInt CLDim )
{
    std::vector<F> work;
    if( onLeft )
        work.resize( n );
    else
        work.resize( m );
    ApplyReflector( onLeft, m, n, v, vInc, tau, C, CLDim, work.data() );
}

template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const float* v, BlasInt vInc,
  const float& tau,
        float* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<float>* v, BlasInt vInc,
  const Complex<float>& tau,
        Complex<float>* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const double* v, BlasInt vInc,
  const double& tau,
        double* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<double>* v, BlasInt vInc,
  const Complex<double>& tau,
        Complex<double>* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const DoubleDouble* v, BlasInt vInc,
  const DoubleDouble& tau,
        DoubleDouble* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<DoubleDouble>* v, BlasInt vInc,
  const Complex<DoubleDouble>& tau,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const QuadDouble* v, BlasInt vInc,
  const QuadDouble& tau,
        QuadDouble* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<QuadDouble>* v, BlasInt vInc,
  const Complex<QuadDouble>& tau,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Quad* v, BlasInt vInc,
  const Quad& tau,
        Quad* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<Quad>* v, BlasInt vInc,
  const Complex<Quad>& tau,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const BigFloat* v, BlasInt vInc,
  const BigFloat& tau,
        BigFloat* C, BlasInt CLDim );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<BigFloat>* v, BlasInt vInc,
  const Complex<BigFloat>& tau,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const float* v, BlasInt vInc,
  const float& tau,
        float* C, BlasInt CLDim,
        float* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<float>* v, BlasInt vInc,
  const Complex<float>& tau,
        Complex<float>* C, BlasInt CLDim,
        Complex<float>* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const double* v, BlasInt vInc,
  const double& tau,
        double* C, BlasInt CLDim,
        double* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<double>* v, BlasInt vInc,
  const Complex<double>& tau,
        Complex<double>* C, BlasInt CLDim,
        Complex<double>* work );
#ifdef EL_HAVE_QD
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const DoubleDouble* v, BlasInt vInc,
  const DoubleDouble& tau,
        DoubleDouble* C, BlasInt CLDim,
        DoubleDouble* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<DoubleDouble>* v, BlasInt vInc,
  const Complex<DoubleDouble>& tau,
        Complex<DoubleDouble>* C, BlasInt CLDim,
        Complex<DoubleDouble>* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const QuadDouble* v, BlasInt vInc,
  const QuadDouble& tau,
        QuadDouble* C, BlasInt CLDim,
        QuadDouble* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<QuadDouble>* v, BlasInt vInc,
  const Complex<QuadDouble>& tau,
        Complex<QuadDouble>* C, BlasInt CLDim,
        Complex<QuadDouble>* work );
#endif
#ifdef EL_HAVE_QUAD
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Quad* v, BlasInt vInc,
  const Quad& tau,
        Quad* C, BlasInt CLDim,
        Quad* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<Quad>* v, BlasInt vInc,
  const Complex<Quad>& tau,
        Complex<Quad>* C, BlasInt CLDim,
        Complex<Quad>* work );
#endif
#ifdef EL_HAVE_MPC
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const BigFloat* v, BlasInt vInc,
  const BigFloat& tau,
        BigFloat* C, BlasInt CLDim,
        BigFloat* work );
template void ApplyReflector
( bool onLeft, BlasInt m, BlasInt n,
  const Complex<BigFloat>* v, BlasInt vInc,
  const Complex<BigFloat>& tau,
        Complex<BigFloat>* C, BlasInt CLDim,
        Complex<BigFloat>* work );
#endif

// Compute the EVD of a symmetric tridiagonal matrix
// =================================================

BlasInt SymmetricTridiagEigWrapper
( char job,
  char range,
  BlasInt n,
  float* d,
  float* e,
  float vl, float vu,
  BlasInt il, BlasInt iu,
  float absTol,
  float* w,
  float* Z, BlasInt ldZ )
{
    DEBUG_CSE
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    float workDummy;
    EL_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<float> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(sstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sstevr's failed");
    return m;
}

BlasInt SymmetricTridiagEigWrapper
( char job, char range, BlasInt n, double* d, double* e, double vl, double vu,
  BlasInt il, BlasInt iu, double absTol, double* w, double* Z, BlasInt ldZ )
{
    DEBUG_CSE
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    double workDummy;
    EL_LAPACK(dstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<double> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(dstevr)
    ( &job, &range, &n, d, e, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dstevr's failed");
    return m;
}

// Compute eigenvalues
// -------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'N', 'A', n, d, e, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'N', 'A', n, d, e, 0, 0, 0, 0, absTol, w, 0, 1 );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float vl, float vu, 
  float absTol )
{
    DEBUG_CSE
    return SymmetricTridiagEigWrapper
    ( 'N', 'V', n, d, e, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double vl, double vu, 
  double absTol )
{
    DEBUG_CSE
    return SymmetricTridiagEigWrapper
    ( 'N', 'V', n, d, e, vl, vu, 0, 0, absTol, w, 0, 1 );
}

// Index range
// ^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, BlasInt il, BlasInt iu, 
  float absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'N', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, BlasInt il, BlasInt iu, 
  double absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'N', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}

// Compute eigenpairs
// ------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ, 
  float absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'V', 'A', n, d, e, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ, 
  double absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'V', 'A', n, d, e, 0, 0, 0, 0, absTol, w, Z, ldZ );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_CSE
    return SymmetricTridiagEigWrapper
    ( 'V', 'V', n, d, e, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_CSE
    return SymmetricTridiagEigWrapper
    ( 'V', 'V', n, d, e, vl, vu, 0, 0, absTol, w, Z, ldZ );
}

// Index range
// ^^^^^^^^^^^^^
void SymmetricTridiagEig
( BlasInt n, float* d, float* e, float* w, float* Z, BlasInt ldZ, 
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'V', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void SymmetricTridiagEig
( BlasInt n, double* d, double* e, double* w, double* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_CSE
    SymmetricTridiagEigWrapper
    ( 'V', 'I', n, d, e, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}

// Compute the EVD of a Hermitian matrix
// =====================================

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, float* A, BlasInt ldA, 
  float vl, float vu, BlasInt il, BlasInt iu, float absTol, 
  float* w, float* Z, BlasInt ldZ )
{
    DEBUG_CSE
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    float workDummy;
    EL_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &iWorkDummy, &iWorkSize,
      &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<float> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(ssyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("ssyevr's failed");
    return m;
}

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, double* A, BlasInt ldA, 
  double vl, double vu, BlasInt il, BlasInt iu, double absTol, 
  double* w, double* Z, BlasInt ldZ )
{
    DEBUG_CSE
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    double workDummy;
    EL_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, 
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy;
    iWorkSize = iWorkDummy;
    vector<double> work(workSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(dsyevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dsyevr's failed");
    return m;
}

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, scomplex* A, BlasInt ldA, 
  float vl, float vu, BlasInt il, BlasInt iu, float absTol, 
  float* w, scomplex* Z, BlasInt ldZ )
{
    DEBUG_CSE
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, rWorkSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    float rWorkDummy;
    scomplex workDummy;
    EL_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &rWorkDummy, &rWorkSize,
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy.real();
    rWorkSize = rWorkDummy;
    iWorkSize = iWorkDummy;
    vector<scomplex> work(workSize);
    vector<float> rWork(rWorkSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(cheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      rWork.data(), &rWorkSize, iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cheevr's failed");
    return m;
}

BlasInt HermitianEigWrapper
( char job, char range, char uplo, BlasInt n, dcomplex* A, BlasInt ldA, 
  double vl, double vu, BlasInt il, BlasInt iu, double absTol, 
  double* w, dcomplex* Z, BlasInt ldZ )
{
    DEBUG_CSE
    if( n == 0 )
        return 0;

    vector<BlasInt> isuppZ( 2*n );

    BlasInt workSize=-1, rWorkSize=-1, iWorkSize=-1, m, info;
    BlasInt iWorkDummy;
    double rWorkDummy;
    dcomplex workDummy;
    EL_LAPACK(zheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), &workDummy, &workSize, &rWorkDummy, &rWorkSize,
      &iWorkDummy, &iWorkSize, &info );

    workSize = workDummy.real();
    rWorkSize = rWorkDummy;
    iWorkSize = iWorkDummy;
    vector<dcomplex> work(workSize);
    vector<double> rWork(rWorkSize);
    vector<BlasInt> iWork(iWorkSize);
    EL_LAPACK(zheevr)
    ( &job, &range, &uplo, &n, A, &ldA, &vl, &vu, &il, &iu, &absTol, &m,
      w, Z, &ldZ, isuppZ.data(), work.data(), &workSize, 
      rWork.data(), &rWorkSize, iWork.data(), &iWorkSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zheevr's failed");
    return m;
}

// Compute the eigenvalues
// -----------------------

// All eigenvalues
// ^^^^^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w, double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w, float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w, double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, 0, 1 );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, 
  float vl, float vu, float absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w, 
  double vl, double vu, double absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w, 
  float vl, float vu, float absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}
BlasInt HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w, 
  double vl, double vu, double absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'N', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, 0, 1 );
}

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, 
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, double* A, BlasInt ldA, double* w, 
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, scomplex* A, BlasInt ldA, float* w, 
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}
void HermitianEig
( char uplo, BlasInt n, dcomplex* A, BlasInt ldA, double* w, 
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'N', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, 0, 1 );
}

// Compute the eigenpairs
// ----------------------

// All eigenpairs
// ^^^^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, float* A, BlasInt ldA, float* w, float* Z, BlasInt ldZ, 
  float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA, double* w, double* Z, BlasInt ldZ,
  double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA, float* w, scomplex* Z, BlasInt ldZ,
  float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA, double* w, dcomplex* Z, BlasInt ldZ,
  double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'A', uplo, n, A, ldA, 0, 0, 0, 0, absTol, w, Z, ldZ );
}

// Floating-point range
// ^^^^^^^^^^^^^^^^^^^^
BlasInt HermitianEig
( char uplo, BlasInt n, 
  float* A, BlasInt ldA, float* w, float* Z, BlasInt ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA, double* w, double* Z, BlasInt ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA, float* w, scomplex* Z, BlasInt ldZ,
  float vl, float vu, float absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}
BlasInt HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA, double* w, dcomplex* Z, BlasInt ldZ,
  double vl, double vu, double absTol )
{
    DEBUG_CSE
    return HermitianEigWrapper
    ( 'V', 'V', uplo, n, A, ldA, vl, vu, 0, 0, absTol, w, Z, ldZ );
}

// Index range
// ^^^^^^^^^^^
void HermitianEig
( char uplo, BlasInt n, 
  float* A, BlasInt ldA, float* w, float* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  double* A, BlasInt ldA, double* w, double* Z, BlasInt ldZ, 
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  scomplex* A, BlasInt ldA, float* w, scomplex* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, float absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}
void HermitianEig
( char uplo, BlasInt n, 
  dcomplex* A, BlasInt ldA, double* w, dcomplex* Z, BlasInt ldZ,
  BlasInt il, BlasInt iu, double absTol )
{
    DEBUG_CSE
    HermitianEigWrapper
    ( 'V', 'I', uplo, n, A, ldA, 0, 0, il+1, iu+1, absTol, w, Z, ldZ );
}

// Bidiagonal DQDS for singular values
// ===================================

void BidiagDQDS( BlasInt n, float* d, float* e )
{
    DEBUG_CSE
    BlasInt info;
    vector<float> work( 4*n );
    EL_LAPACK(slasq1)( &n, d, e, work.data(), &info );
    if( info != 0 )
    {
        ostringstream msg;
        if( info < 0 )
            msg << "Argument " << -info << " had an illegal value";
        else if( info == 1 )
            msg << "A split was marked in a positive value in E";
        else if( info == 2 )
            msg << "Current block of Z not bidiagonalized after 30*k its";
        else if( info == 3 )
            msg << "Termination criterion of outer while loop not met";
        RuntimeError( msg.str() );
    } 
}

void BidiagDQDS( BlasInt n, double* d, double* e )
{
    DEBUG_CSE
    BlasInt info;
    vector<double> work( 4*n );
    EL_LAPACK(dlasq1)( &n, d, e, work.data(), &info );
    if( info != 0 )
    {
        ostringstream msg;
        if( info < 0 )
            msg << "Argument " << -info << " had an illegal value";
        else if( info == 1 )
            msg << "A split was marked in a positive value in E";
        else if( info == 2 )
            msg << "Current block of Z not bidiagonalized after 30*k its";
        else if( info == 3 )
            msg << "Termination criterion of outer while loop not met";
        RuntimeError( msg.str() );
    } 
}

// Bidiagonal QR algorithm for SVD
// ===============================

void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVT, BlasInt numRowsU,
  float* d, float* e, float* VTrans, BlasInt ldVT, float* U, BlasInt ldU )
{
    DEBUG_CSE
    if( n==0 )
        return;

    BlasInt info;
    float* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<float> work( 4*n );
    EL_LAPACK(sbdsqr)
    ( &uplo, &n, &numColsVT, &numRowsU, &numColsC,
      d, e,
      VTrans, &ldVT,
      U, &ldU,
      C, &ldC,
      work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sbdsqr had ",info," elements of e not converge");
}

void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVT, BlasInt numRowsU, 
  double* d, double* e, double* VTrans, BlasInt ldVT, double* U, BlasInt ldU )
{
    DEBUG_CSE
    if( n==0 )
        return;

    BlasInt info;
    double* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<double> work( 4*n );
    EL_LAPACK(dbdsqr)
    ( &uplo, &n, &numColsVT, &numRowsU, &numColsC,
      d, e, 
      VTrans, &ldVT,
      U, &ldU,
      C, &ldC,
      work.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dbdsqr had ",info," elements of e not converge");
}

void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVH, BlasInt numRowsU, 
  float* d, float* e,
  scomplex* VH, BlasInt ldVH,
  scomplex* U,  BlasInt ldU )
{
    DEBUG_CSE
    if( n==0 )
        return;

    BlasInt info;
    scomplex* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<float> realWork;
    if( numColsVH == 0 && numRowsU == 0 && numColsC == 0 )
        realWork.resize( 2*n );
    else
        realWork.resize( Max(1,4*n-4) );
    EL_LAPACK(cbdsqr)
    ( &uplo, &n, &numColsVH, &numRowsU, &numColsC,
      d, e,
      VH, &ldVH,
      U, &ldU,
      C, &ldC,
      realWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cbdsqr had ",info," elements of e not converge");
}

void BidiagSVDQRAlg
( char uplo, BlasInt n, BlasInt numColsVH, BlasInt numRowsU, 
  double* d, double* e,
  dcomplex* VH, BlasInt ldVH,
  dcomplex* U,  BlasInt ldU )
{
    DEBUG_CSE
    if( n==0 )
        return;

    BlasInt info;
    dcomplex* C=0;
    const BlasInt numColsC=0, ldC=1;
    vector<double> realWork;
    if( numColsVH == 0 && numRowsU == 0 && numColsC == 0 )
        realWork.resize( 2*n );
    else
        realWork.resize( Max(1,4*n-4) );
    EL_LAPACK(zbdsqr)
    ( &uplo, &n, &numColsVH, &numRowsU, &numColsC,
      d, e,
      VH, &ldVH,
      U, &ldU,
      C, &ldC,
      realWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zbdsqr had ",info," elements of e not converge");
}

// Divide and Conquer SVD
// ======================

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  float* A, BlasInt ldA, 
  float* s,
  float* U, BlasInt ldu,
  float* VTrans, BlasInt ldvt,
  bool thin )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    float workDummy;
    const BlasInt k = Min(m,n);
    vector<BlasInt> iWork(8*k);

    EL_LAPACK(sgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      &workDummy, &workSize,
      iWork.data(),
      &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(sgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      work.data(), &workSize,
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesdd's updating process failed");
}

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  double* A, BlasInt ldA, 
  double* s,
  double* U, BlasInt ldu,
  double* VTrans, BlasInt ldvt,
  bool thin )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    double workDummy;
    const BlasInt k = Min(m,n);
    vector<BlasInt> iWork(8*k);

    EL_LAPACK(dgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      &workDummy, &workSize,
      iWork.data(),
      &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt,
      work.data(), &workSize,
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesdd's updating process failed");
}

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  scomplex* A, BlasInt ldA, 
  float* s,
  scomplex* U, BlasInt ldu,
  scomplex* VH, BlasInt ldva,
  bool thin )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    const BlasInt k = Min(m,n);
    const BlasInt K = Max(m,n);
    const BlasInt rWorkSize = k*Max(5*k+7,2*K+2*k+1);
    vector<float> rWork(rWorkSize);
    vector<BlasInt> iWork(8*k);

    scomplex workDummy;
    EL_LAPACK(cgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      &workDummy, &workSize,
      rWork.data(),
      iWork.data(),
      &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(cgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      work.data(), &workSize,
      rWork.data(),
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesdd's updating process failed");
}

void DivideAndConquerSVD
( BlasInt m, BlasInt n,
  dcomplex* A, BlasInt ldA, 
  double* s,
  dcomplex* U, BlasInt ldu,
  dcomplex* VH, BlasInt ldva,
  bool thin )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobz = ( thin ? 'S' : 'A' );
    BlasInt workSize=-1, info;
    dcomplex workDummy;
    const BlasInt k = Min(m,n);
    const BlasInt K = Max(m,n);
    const BlasInt rWorkSize = k*Max(5*k+7,2*K+2*k+1);
    vector<double> rWork(rWorkSize);
    vector<BlasInt> iWork(8*k);

    EL_LAPACK(zgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      &workDummy, &workSize,
      rWork.data(),
      iWork.data(),
      &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zgesdd)
    ( &jobz, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva,
      work.data(), &workSize,
      rWork.data(),
      iWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesdd's updating process failed");
}

// QR-algorithm SVD
// ================

void QRSVD
( BlasInt m, BlasInt n,
  float* A, BlasInt ldA, 
  float* s,
  float* U, BlasInt ldu,
  float* VTrans, BlasInt ldvt,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVT= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    float workDummy;

    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      &workDummy, &workSize,
      &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void QRSVD
( BlasInt m, BlasInt n,
  double* A, BlasInt ldA, 
  double* s,
  double* U, BlasInt ldu,
  double* VTrans, BlasInt ldvt,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVT= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    double workDummy;

    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      &workDummy, &workSize,
      &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VTrans, &ldvt, 
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void QRSVD
( BlasInt m, BlasInt n,
  scomplex* A, BlasInt ldA, 
  float* s,
  scomplex* U, BlasInt ldu,
  scomplex* VH, BlasInt ldva,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVH= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    const BlasInt k = Min(m,n);
    vector<float> rWork(5*k);

    scomplex workDummy;
    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      &workDummy, &workSize,
      rWork.data(),
      &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      work.data(), &workSize,
      rWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void QRSVD
( BlasInt m, BlasInt n,
  dcomplex* A, BlasInt ldA, 
  double* s,
  dcomplex* U, BlasInt ldu,
  dcomplex* VH, BlasInt ldva,
  bool thin, bool avoidU, bool avoidV )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU= ( avoidU ? 'N' : ( thin ? 'S' : 'A' ) ),
               jobVH= ( avoidV ? 'N' : ( thin ? 'S' : 'A' ) );
    BlasInt workSize=-1, info;
    dcomplex workDummy;
    const BlasInt k = Min(m,n);
    vector<double> rWork(5*k);

    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      &workDummy, &workSize,
      rWork.data(),
      &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n,
      A, &ldA,
      s,
      U, &ldu,
      VH, &ldva, 
      work.data(), &workSize,
      rWork.data(),
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Compute singular values (with DQDS)
// ===================================

void SVD( BlasInt m, BlasInt n, float* A, BlasInt ldA, float* s )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVT='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    float workDummy;

    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<float> work(workSize);
    EL_LAPACK(sgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("sgesvd's updating process failed");
}

void SVD( BlasInt m, BlasInt n, double* A, BlasInt ldA, double* s )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVT='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    double workDummy;

    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, &info );

    workSize = workDummy;
    vector<double> work(workSize);
    EL_LAPACK(dgesvd)
    ( &jobU, &jobVT, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("dgesvd's updating process failed");
}

void SVD( BlasInt m, BlasInt n, scomplex* A, BlasInt ldA, float* s )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVH='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    scomplex workDummy;
    const BlasInt k = Min(m,n);
    vector<float> rWork(5*k);

    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<scomplex> work(workSize);
    EL_LAPACK(cgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("cgesvd's updating process failed");
}

void SVD( BlasInt m, BlasInt n, dcomplex* A, BlasInt ldA, double* s )
{
    DEBUG_CSE
    if( m==0 || n==0 )
        return;

    const char jobU='N', jobVH='N';
    BlasInt fakeLDim=1, workSize=-1, info;
    dcomplex workDummy;
    const BlasInt k = Min(m,n);
    vector<double> rWork(5*k);

    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      &workDummy, &workSize, rWork.data(), &info );

    workSize = workDummy.real();
    vector<dcomplex> work(workSize);
    EL_LAPACK(zgesvd)
    ( &jobU, &jobVH, &m, &n, A, &ldA, s, 0, &fakeLDim, 0, &fakeLDim, 
      work.data(), &workSize, rWork.data(), &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," had an illegal value");
    else if( info > 0 )
        RuntimeError("zgesvd's updating process failed");
}

// Reduce a general square matrix to upper Hessenberg form
// =======================================================

void Hessenberg
( BlasInt n,
  float* A, BlasInt ldA,
  float* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy;

    // Reduce to Hessenberg form
    vector<float> work( workSize );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

void Hessenberg
( BlasInt n,
  double* A, BlasInt ldA,
  double* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy;

    // Reduce to Hessenberg form
    vector<double> work( workSize );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

void Hessenberg
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy.real();

    // Reduce to Hessenberg form
    vector<scomplex> work( workSize );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

void Hessenberg
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy.real();

    // Reduce to Hessenberg form
    vector<dcomplex> work( workSize );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

// Generate the unitary matrix used within the reduction to Hessenberg form
// ------------------------------------------------------------------------

void HessenbergGenerateUnitary
( BlasInt n,
  float* A, BlasInt ldA,
  const float* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    EL_LAPACK(sorghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy;

    // Generate the unitary matrix
    vector<float> work( workSize );
    EL_LAPACK(sorghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

void HessenbergGenerateUnitary
( BlasInt n,
  double* A, BlasInt ldA,
  const double* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    EL_LAPACK(dorghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy;

    // Generate the unitary matrix
    vector<double> work( workSize );
    EL_LAPACK(dorghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

void HessenbergGenerateUnitary
( BlasInt n,
  scomplex* A, BlasInt ldA,
  const scomplex* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    EL_LAPACK(cunghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy.real();

    // Generate the unitary matrix
    vector<scomplex> work( workSize );
    EL_LAPACK(cunghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

void HessenbergGenerateUnitary
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  const dcomplex* tau )
{
    // Query the workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    EL_LAPACK(zunghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      &workDummy, &workSize,
      &info );
    workSize = workDummy.real();

    // Generate the unitary matrix
    vector<dcomplex> work( workSize );
    EL_LAPACK(zunghr)
    ( &n, &ilo, &ihi,
      A, &ldA,
      tau,
      work.data(), &workSize,
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");
}

// Solve a 2x2 linear system using LU with full pivoting, perturbing the
// matrix as necessary to ensure sufficiently large pivots.
// NOTE: This is primarily a helper function for SmallSylvester
template<typename Real>
bool Solve2x2FullPiv
( const Real* A,
        Real* b,
        Real& scale,
  const Real& smallNum,
  const Real& minPiv )
{
    DEBUG_CSE
    const Real one(1), two(2);
    // Avoid tedious index calculations for the 2x2 LU with full pivoting
    // using cached tables in the same manner as LAPACK's {s,d}lasy2
    const BlasInt lambda21Ind[4] = { 1, 0, 3, 2 };
    const BlasInt ups12Ind[4] = { 2, 3, 0, 1 };
    const BlasInt ups22Ind[4] = { 3, 2, 1, 0 };
    const bool XSwapTable[4] = { false, false, true,  true };
    const bool BSwapTable[4] = { false, true,  false, true };

    bool perturbed = false;
    BlasInt iPiv = blas::MaxInd( 4, A, 1 );
    Real ups11 = A[iPiv];
    if( Abs(ups11) < minPiv )
    {
        ups11 = minPiv;
        perturbed = true;
    }
    Real ups12    = A[ups12Ind[iPiv]];
    Real lambda21 = A[lambda21Ind[iPiv]] / ups11;
    Real ups22    = A[ups22Ind[iPiv]] - ups12*lambda21;
    if( Abs(ups22) < minPiv )
    {
        ups22 = minPiv;
        perturbed = true;
    }
    if( BSwapTable[iPiv] )
    {
        Real tmp = b[1];
        b[1] = b[0] - lambda21*tmp;
        b[0] = tmp;
    }
    else
    {
        b[1] -= lambda21*b[0];
    }
    scale = one;
    if( (two*smallNum)*Abs(b[1]) > Abs(ups22) ||
        (two*smallNum)*Abs(b[0]) > Abs(ups11) )
    {
        b[0] *= scale;
        b[1] *= scale;
    }
    b[1] /= ups22;
    b[0] = b[0] / ups11 - (ups12/ups11)*b[1];
    if( XSwapTable[iPiv] )
    {
        Real tmp = b[1];
        b[1] = b[0];
        b[0] = tmp;
    }
    return perturbed;
}
template bool Solve2x2FullPiv
( const float* A,
        float* b,
        float& scale,
  const float& smallNum,
  const float& minPiv );
template bool Solve2x2FullPiv
( const double* A,
        double* b,
        double& scale,
  const double& smallNum,
  const double& minPiv );
#ifdef EL_HAVE_QUAD
template bool Solve2x2FullPiv
( const Quad* A,
        Quad* b,
        Quad& scale,
  const Quad& smallNum,
  const Quad& minPiv );
#endif
#ifdef EL_HAVE_QD
template bool Solve2x2FullPiv
( const DoubleDouble* A,
        DoubleDouble* b,
        DoubleDouble& scale,
  const DoubleDouble& smallNum,
  const DoubleDouble& minPiv );
template bool Solve2x2FullPiv
( const QuadDouble* A,
        QuadDouble* b,
        QuadDouble& scale,
  const QuadDouble& smallNum,
  const QuadDouble& minPiv );
#endif
#ifdef EL_HAVE_MPC
template bool Solve2x2FullPiv
( const BigFloat* A,
        BigFloat* b,
        BigFloat& scale,
  const BigFloat& smallNum,
  const BigFloat& minPiv );
#endif

// Solve a 4x4 linear system using LU with full pivoting, perturbing the
// matrix as necessary to ensure sufficiently large pivots.
// NOTE: This is primarily a helper function for SmallSylvester
template<typename Real>
bool Solve4x4FullPiv
(       Real* A,
        Real* b,
        Real& scale,
  const Real& smallNum,
  const Real& minPiv )
{
    DEBUG_CSE
    const Real zero(0), one(1), eight(8);
    bool perturbed = false;

    BlasInt iPiv, jPiv;
    BlasInt p[4];
    for( BlasInt i=0; i<3; ++i )
    {
        iPiv=jPiv=0;
        Real AMax = zero;
        for( BlasInt iMax=i; iMax<4; ++iMax )
        {
            for( BlasInt jMax=i; jMax<4; ++jMax )
            {
                if( Abs(A[iMax+jMax*4]) >= AMax )
                {
                    AMax = Abs(A[iMax+jMax*4]);
                    iPiv = iMax;
                    jPiv = jMax;
                }
            }
        }

        if( iPiv != i )
        {
            blas::Swap( 4, &A[iPiv], 4, &A[i], 4 );
            Real tmp = b[i];
            b[i] = b[iPiv];
            b[iPiv] = tmp;
        }
        if( jPiv != i )
        {
            blas::Swap( 4, &A[jPiv*4], 1, &A[i*4], 1 );
        }
        p[i] = jPiv;
        if( Abs(A[i+i*4]) < minPiv )
        {
            A[i+i*4] = minPiv;
            perturbed = true;
        }
        for( BlasInt j=i+1; j<4; ++j )
        {
            A[j+i*4] /= A[i+i*4];
            b[j] -= A[j+i*4]*b[i];
            for( BlasInt k=i+1; k<4; ++k )
            {
                A[j+k*4] -= A[j+i*4]*A[i+k*4];
            }
        }
    }
    if( Abs(A[3+3*4]) < minPiv )
    {
        A[3+3*4] = minPiv;
        perturbed = true;
    }
    scale = one;
    if( (eight*smallNum)*Abs(b[0]) > Abs(A[0+0*4]) ||
        (eight*smallNum)*Abs(b[1]) > Abs(A[1+1*4]) ||
        (eight*smallNum)*Abs(b[2]) > Abs(A[2+2*4]) ||
        (eight*smallNum)*Abs(b[3]) > Abs(A[3+3*4]) )
    {
        const Real xMax = blas::NrmInf( 4, b, 1 );
        scale = (one/eight) / xMax;
        b[0] *= scale;
        b[1] *= scale;
        b[2] *= scale;
        b[3] *= scale;
    }

    for( BlasInt i=0; i<4; ++i )
    {
        BlasInt k = 3-i;
        Real tmp = one / A[k+k*4];
        b[k] *= tmp;
        for( BlasInt j=k+1; j<4; ++j )
            b[k] -= (tmp*A[k+j*4])*b[j];
    }
    for( BlasInt i=0; i<3; ++i )
    {
        if( p[2-i] != 2-i )
        {
            Real tmp = b[2-i];
            b[2-i] = b[p[2-i]];
            b[p[2-i]] = tmp;
        }
    }
    return perturbed;
}
template bool Solve4x4FullPiv
( float* A,
  float* b,
  float& scale,
  const float& smallNum,
  const float& minPiv );
template bool Solve4x4FullPiv
( double* A,
  double* b,
  double& scale,
  const double& smallNum,
  const double& minPiv );
#ifdef EL_HAVE_QUAD
template bool Solve4x4FullPiv
( Quad* A,
  Quad* b,
  Quad& scale,
  const Quad& smallNum,
  const Quad& minPiv );
#endif
#ifdef EL_HAVE_QD
template bool Solve4x4FullPiv
( DoubleDouble* A,
  DoubleDouble* b,
  DoubleDouble& scale,
  const DoubleDouble& smallNum,
  const DoubleDouble& minPiv );
template bool Solve4x4FullPiv
( QuadDouble* A,
  QuadDouble* b,
  QuadDouble& scale,
  const QuadDouble& smallNum,
  const QuadDouble& minPiv );
#endif
#ifdef EL_HAVE_MPC
template bool Solve4x4FullPiv
( BigFloat* A,
  BigFloat* b,
  BigFloat& scale,
  const BigFloat& smallNum,
  const BigFloat& minPiv );
#endif

// Solve a 1x1, 1x2, 2x1, or 2x2 Sylvester equation, 
//
//   op_C(C) X +- X op_D(D) = scale*B,
//
// where op_C(C) is either C or C^T, op_D(D) is either D or D^T,
// and scale in (0,1] is determined by the subroutine.
//
// The fundamental technique is Gaussian Elimination with full pivoting,
// with pivots forced to be sufficiently large (and, if such a perturbation was
// performed, the routine returns 'true').
//
// The analogous LAPACK routines are {s,d}lasy2.
//
template<typename Real>
bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const Real* C, BlasInt CLDim,
  const Real* D, BlasInt DLDim,
  const Real* B, BlasInt BLDim,
        Real& scale,
        Real* X, BlasInt XLDim,
        Real& XInfNorm )
{
    DEBUG_CSE
    const Real one(1);
    const Real epsilon = limits::Epsilon<Real>();
    const Real smallNum = limits::SafeMin<Real>() / epsilon;
    const Real sgn = ( negate ? Real(-1) : Real(1) );

    bool perturbed = false;
    if( nC == 1 && nD == 1 )
    {
        // psi*chi + sgn*chi = beta
        const Real& psi   = C[0+0*CLDim];
        const Real& delta = D[0+0*DLDim];
        const Real& beta  = B[0+0*BLDim];
              Real& chi   = X[0+0*XLDim];

        Real tau = psi + sgn*delta;
        Real tauAbs = Abs(tau);
        if( tauAbs <= smallNum )
        {
            tau = tauAbs = smallNum;
            perturbed = true;
        }
        Real gamma = Abs(beta);
        scale = one;
        if( smallNum*gamma > tauAbs )
            scale = one/gamma;
        chi = (beta*scale) / tau;
        XInfNorm = Abs(chi);
    }
    else if( nC == 1 && nD == 2 )
    {
        // psi*x +- x*op(D) = b,
        //
        // where x = [chi0, chi1], D = |delta00, delta01|, b = [beta0, beta1].
        //                             |delta10, delta11|
        //
        const Real& psi     = C[0+0*CLDim];
        const Real& delta00 = D[0+0*DLDim];
        const Real& delta01 = D[0+1*DLDim];
        const Real& delta10 = D[1+0*DLDim];
        const Real& delta11 = D[1+1*DLDim];
        const Real& beta0   = B[0+0*BLDim];
        const Real& beta1   = B[0+1*BLDim];
              Real& chi0    = X[0+0*XLDim];
              Real& chi1    = X[0+1*XLDim];

        Real maxCD = Max( Abs(delta00), Abs(delta01) );
        maxCD = Max( maxCD, Abs(delta10) );
        maxCD = Max( maxCD, Abs(delta11) );
        maxCD = Max( maxCD, Abs(psi) );
        const Real minPiv = Max( epsilon*maxCD, smallNum );

        // 
        // In the case of no transpositions, solve
        //
        // | psi +- delta00,   +-delta01    | | chi0 | = | beta0 |
        // |    +-delta10,   psi +- delta11 | | chi1 |   | beta1 |
        //
        Real A[4], b[2];
        A[0] = psi + sgn*delta00;
        A[3] = psi + sgn*delta11;
        if( transD )
        {
            A[1] = sgn*delta10;
            A[2] = sgn*delta01;
        }
        else
        {
            A[1] = sgn*delta01;
            A[2] = sgn*delta10;
        }
        b[0] = beta0;
        b[1] = beta1;
        perturbed = Solve2x2FullPiv( A, b, scale, smallNum, minPiv );
        chi0 = b[0];
        chi1 = b[1];
        XInfNorm = Abs(chi0) + Abs(chi1);
    }
    else if( nC == 2 && nD == 1 )
    {
        // op(C)*x +- x*delta = b,
        //
        // where x = |chi0|, C = |psi00, psi01|, b = |beta0|.
        //           |chi1|      |psi10, psi11|      |beta1|
        //
        const Real& psi00 = C[0+0*CLDim];
        const Real& psi01 = C[0+1*CLDim];
        const Real& psi10 = C[1+0*CLDim];
        const Real& psi11 = C[1+1*CLDim];
        const Real& delta = D[0+0*DLDim];
        const Real& beta0 = B[0+0*BLDim];
        const Real& beta1 = B[1+0*BLDim];
              Real& chi0  = X[0+0*XLDim];
              Real& chi1  = X[1+0*XLDim];

        Real maxCD = Max( Abs(psi00), Abs(psi01) );
        maxCD = Max( maxCD, Abs(psi10) );
        maxCD = Max( maxCD, Abs(psi11) );
        maxCD = Max( maxCD, Abs(delta) );
        const Real minPiv = Max( epsilon*maxCD, smallNum );

        // 
        // In the case of no transpositions, solve
        //
        // | psi00 +- delta,      psi01      | | chi0 | = | beta0 |
        // |    psi10,        psi11 +- delta | | chi1 |   | beta1 |
        //
        Real A[4], b[2];
        A[0] = psi00 + sgn*delta;
        A[3] = psi11 + sgn*delta;
        if( transC )
        {
            A[1] = psi01;
            A[2] = psi10;
        }
        else
        {
            A[1] = psi10;
            A[2] = psi01;
        }
        b[0] = beta0;
        b[1] = beta1;
        perturbed = Solve2x2FullPiv( A, b, scale, smallNum, minPiv );
        chi0 = b[0];
        chi1 = b[1];
        XInfNorm = Max( Abs(chi0), Abs(chi1) );
    }
    else if( nC == 2 && nD == 2 )
    {
        // op(C)*X +- X*op(D) = B
        const Real& psi00   = C[0+0*CLDim];
        const Real& psi01   = C[0+1*CLDim];
        const Real& psi10   = C[1+0*CLDim];
        const Real& psi11   = C[1+1*CLDim];
        const Real& delta00 = D[0+0*DLDim];
        const Real& delta01 = D[0+1*DLDim];
        const Real& delta10 = D[1+0*DLDim];
        const Real& delta11 = D[1+1*DLDim];
        const Real& beta00  = B[0+0*BLDim];
        const Real& beta01  = B[0+1*BLDim];
        const Real& beta10  = B[1+0*BLDim];
        const Real& beta11  = B[1+1*BLDim];
              Real& chi00   = X[0+0*XLDim];
              Real& chi01   = X[0+1*XLDim];
              Real& chi10   = X[1+0*XLDim];
              Real& chi11   = X[1+1*XLDim];

        Real maxCD = Max( Abs(psi00), Abs(psi01) );
        maxCD = Max( maxCD, Abs(psi10) );
        maxCD = Max( maxCD, Abs(psi11) );
        maxCD = Max( maxCD, Abs(delta00) );
        maxCD = Max( maxCD, Abs(delta01) );
        maxCD = Max( maxCD, Abs(delta10) );
        maxCD = Max( maxCD, Abs(delta11) );
        const Real minPiv = Max( epsilon*maxCD, smallNum );

        Real A[16], b[4];
        A[0+0*4] = psi00 + sgn*delta00;
        A[1+1*4] = psi11 + sgn*delta00;
        A[2+2*4] = psi00 + sgn*delta11;
        A[3+3*4] = psi11 + sgn*delta11;
        if( transC )
        {
            A[0+1*4] = psi10;
            A[1+0*4] = psi01;
            A[2+3*4] = psi10;
            A[3+2*4] = psi01;
        }
        else
        {
            A[0+1*4] = psi01;
            A[1+0*4] = psi10;
            A[2+3*4] = psi01;
            A[3+2*4] = psi10;
        }
        if( transD )
        {
            A[0+2*4] = sgn*delta01;
            A[1+3*4] = sgn*delta01;
            A[2+0*4] = sgn*delta10;
            A[3+1*4] = sgn*delta10;
        }
        else
        {
            A[0+2*4] = sgn*delta10;
            A[1+3*4] = sgn*delta10;
            A[2+0*4] = sgn*delta01;
            A[3+1*4] = sgn*delta01;
        }
        A[0+3*4] = A[1+2*4] = A[2+1*4] = A[3+0*4] = 0;
        b[0] = beta00;
        b[1] = beta10;
        b[2] = beta01;
        b[3] = beta11;

        perturbed = Solve4x4FullPiv( A, b, scale, smallNum, minPiv );
        chi00 = b[0];
        chi10 = b[1];
        chi01 = b[2];
        chi11 = b[3];
        XInfNorm = Max( Abs(chi00)+Abs(chi01), Abs(chi10)+Abs(chi11) );
    }
    else
        LogicError("Invalid SmallSylvester sizes");

    return perturbed;
}
template bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const float* C, BlasInt CLDim,
  const float* D, BlasInt DLDim,
  const float* B, BlasInt BLDim,
        float& scale,
        float* X, BlasInt XLDim,
        float& XInfNorm );
template bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const double* C, BlasInt CLDim,
  const double* D, BlasInt DLDim,
  const double* B, BlasInt BLDim,
        double& scale,
        double* X, BlasInt XLDim,
        double& XInfNorm );
#ifdef EL_HAVE_QUAD
template bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const Quad* C, BlasInt CLDim,
  const Quad* D, BlasInt DLDim,
  const Quad* B, BlasInt BLDim,
        Quad& scale,
        Quad* X, BlasInt XLDim,
        Quad& XInfNorm );
#endif
#ifdef EL_HAVE_QD
template bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const DoubleDouble* C, BlasInt CLDim,
  const DoubleDouble* D, BlasInt DLDim,
  const DoubleDouble* B, BlasInt BLDim,
        DoubleDouble& scale,
        DoubleDouble* X, BlasInt XLDim,
        DoubleDouble& XInfNorm );
template bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const QuadDouble* C, BlasInt CLDim,
  const QuadDouble* D, BlasInt DLDim,
  const QuadDouble* B, BlasInt BLDim,
        QuadDouble& scale,
        QuadDouble* X, BlasInt XLDim,
        QuadDouble& XInfNorm );
#endif
#ifdef EL_HAVE_MPC
template bool SmallSylvester
( bool transC,
  bool transD,
  bool negate,
  BlasInt nC, BlasInt nD,
  const BigFloat* C, BlasInt CLDim,
  const BigFloat* D, BlasInt DLDim,
  const BigFloat* B, BlasInt BLDim,
        BigFloat& scale,
        BigFloat* X, BlasInt XLDim,
        BigFloat& XInfNorm );
#endif

namespace adjacent_schur {

// See the paper
//
//   Zhaojun Bai and James Demmel,
//   "On swapping diagonal blocks in real Schur form",
//   Linear Algebra and its Applications, 186:73--95, 1993,
//
// and the corresponding implementations in LAPACK's {s,d}laexc.
//
template<typename Real>
void Helper
( bool wantSchurVecs,
  BlasInt n,
  Real* T, BlasInt TLDim, 
  Real* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  Real* work,
  bool testAccuracy )
{
    DEBUG_CSE
    // Uphold LAPACK's conventions for allowing n1, n2 in {0,1,2}
    DEBUG_ONLY(
      if( n1 < 0 || n1 > 2 )
          LogicError("n1 must be in {0,1,2}");
      if( n2 < 0 || n2 > 2 )
          LogicError("n2 must be in {0,1,2}");
    )
    if( n == 0 || n1 == 0 || n2 == 0 || j1+n1 >= n )
        return;

    // As discussed by Bai/Demmel (building on a theorem of Ng/Peyton),
    // given the relationship
    //
    //   | A11, A12 | = | I, -X | | A11, 0   | | I, X |,
    //   |  0,  A22 |   | 0,  I | | 0,   A22 | | 0, I |
    //
    // where X is the solution of the Sylvester equation
    //
    //   A11 X - X A22 = A12,
    //
    // a Householder QR factorization of [-X; I] or [I, X], depending upon
    // which requires fewer reflections, can be used to swap A11 and A22
    // via a similarity transformation.
    //
    // For example,
    //
    //   Q' | -X | = | R | 
    //      |  I |   | 0 |    
    //
    // implies
    //
    //   | Q11', Q21' | | I, -X | = | Q11', R |
    //   | Q12', Q22' | | 0,  I |   | Q12', 0 |
    //
    // and therefore
    //
    //   Q' A Q
    //
    //    = | Q11', R | | A11, 0   | | Q11', R |^{-1}
    //      | Q12', 0 | | 0,   A22 | | Q12', 0 |
    //
    //    = | Q11', R | | A11, 0   | | 0,       inv(Q12)'             |,
    //      | Q12', 0 | | 0,   A22 | | inv(R), -inv(R) Q11' inv(Q12)' |
    //
    //    = | R A22 inv(R), Q11' A11 inv(Q12)' - R A22 inv(R) Q11' inv(Q12)' |.
    //      |       0,      Q12' A11 inv(Q12)'                               |
    //
    // On the other hand, should X have fewer rows than columns, it might
    // be preferable to instead pick a unitary Q such that
    //
    //   | I, X | | Q11, Q12 | = | 0,   R   |,
    //   | 0, I | | Q21, Q22 |   | Q21, Q22 |
    //
    // and therefore
    //
    //   Q' A Q 
    //
    //    = | 0,   R   |^{-1} | A11, 0   | | 0,   R   |
    //      | Q21, Q22 |      | 0,   A22 | | Q21, Q22 |
    //
    //    = | -inv(Q21)' Q22 inv(R), inv(Q21) | | A11,  0  | | 0,    R  |
    //      | inv(R),                    0    | |  0,  A22 | | Q21, Q22 |
    // 
    //    = | -inv(Q21)' Q22 inv(R) A11, inv(Q21) A22 | | 0,    R  |
    //      |         inv(R) A11,               0     | | Q21, Q22 |
    //
    //    = |inv(Q21) A22 Q21, -inv(Q21)' Q22 inv(R) A11 R + inv(Q21) A22 Q22|.
    //      |         0,                     inv(R) A11 R                    |
    //

    if( n1 == 1 && n2 == 1 )
    {
        Real tau11 = T[ j1   + j1   *TLDim];
        Real tau12 = T[ j1   +(j1+1)*TLDim];
        Real tau22 = T[(j1+1)+(j1+1)*TLDim];

        // Force the bottom-left entry of the similarity transformation
        //
        //   |  c, s | | tau11 tau12 | | c, -s |
        //   | -s  c | |   0   tau22 | | s,  c |
        //
        // to remain zero (and therefore swapping tau11 and tau22 and preserving
        // tau12)
        Real c, s;
        Givens( tau12, tau22-tau11, c, s );

        if( j1+2 < n )
        {
            blas::Rot
            ( n-(j1+2),
              &T[ j1   +(j1+2)*TLDim], TLDim,
              &T[(j1+1)+(j1+2)*TLDim], TLDim, c, s );
        }
        blas::Rot
        ( j1,
          &T[0+ j1   *TLDim], 1,
          &T[0+(j1+1)*TLDim], 1, c, s );

        T[ j1   + j1   *TLDim] = tau22;
        T[(j1+1)+(j1+1)*TLDim] = tau11;

        if( wantSchurVecs )
        {
            blas::Rot
            ( n,
              &Q[0+ j1   *QLDim], 1,
              &Q[0+(j1+1)*QLDim], 1, c, s );
        }
    }
    else
    {
        Real D[16], X[4];
        const BlasInt XLDim = 2;
        const BlasInt nSum = n1 + n2; 
        for( BlasInt j=0; j<nSum; ++j )
            for( BlasInt i=0; i<nSum; ++i )
                D[i+j*nSum] = T[(j1+i)+(j1+j)*TLDim];
        const Real DMax = blas::NrmInf( nSum*nSum, D, 1 );

        const Real epsilon = limits::Epsilon<Real>();
        const Real smallNum = limits::SafeMin<Real>() / epsilon;
        // NOTE: This value was increased from 10 because of slight failures
        const Real fudge = 20;
        const Real thresh = Max( fudge*epsilon*DMax, smallNum );

        // Solve the Sylvester equation T11*X - X*T22 = scale*T12 for X
        Real scale, XInfNorm;
        const bool transT11=false, transT22=false, negate=true;
        SmallSylvester
        ( transT11, transT22, negate, n1, n2,
          &D[0 +0 *nSum], nSum,
          &D[n1+n1*nSum], nSum,
          &D[0 +n1*nSum], nSum,
          scale, X, XLDim, XInfNorm );

        Real innerProd;
        const Real zero(0);
        if( n1 == 1 && n2 == 2 )
        {
            const Real tau11 = T[j1+j1*TLDim];

            // Compute the Householder reflection which satisfies
            //
            //   | 1, x | | gamma11, q12 | = | 0,    r  |,
            //   | 0, I | |   q21,   Q22 |   | q21, Q22 |
            //
            // where x is a row vector, gamma11 is scalar, q12 is a row vector,
            // q21 is a column vector, and r is a row vector.
            //
            // This can be accomplished via the Householder reflection
            //
            //   (I - phi [nu0; nu1; 1] [nu0, nu1, 1]) | 1    | = | 0   |
            //                                         | chi0 |   | 0   |
            //                                         | chi1 |   | rho |
            // // which is simply a permutation of the standard formulation
            // (and this fact is exploited by LAPACK).
            Real v[3];
            v[0] = scale;
            v[1] = X[0*XLDim];
            v[2] = X[1*XLDim];
            Real phi = Reflector( 3, v[2], v, 1 );
            v[2] = 1; 

            if( testAccuracy )
            {
                // Perform the candidate rotation out-of-place on D
                // ------------------------------------------------ 
                // Apply the rotation from the left,
                //   D := (I - phi v v') D = D - v (phi v' D)
                ApplyReflector( true, 3, 3, v, 1, phi, D, nSum, work );
                // Apply the rotation from the right,
                //   D := D (I - phi v v') = D - (phi D v) v'
                ApplyReflector( false, 3, 3, v, 1, phi, D, nSum, work );

                // Throw an exception if the rotation would be too inaccurate.
                // As in LAPACK, rather than simply measuring the size of the 
                // bottom-left block of the rotation of D (which should ideally
                // be zero), we also take into account our knowledge that the
                // bottom-right entry of D should by T(j1,j1)
                Real errMeasure = Max( Abs(D[2+0*3]), Abs(D[2+1*3]) );
                errMeasure = Max( errMeasure, Abs(D[2+2*3]-tau11) );
                if( errMeasure > thresh )
                    RuntimeError
                    ("Unacceptable Schur block swap: errMeasure, ",errMeasure,
                     " was greater than ",thresh);
            }

            // Perform the rotation on T
            // ------------------------- 
            // Apply the rotation from the left,
            //   T := (I - phi v v') T = T - v (phi v' T)
            ApplyReflector
            ( true, 3, n-j1, v, 1, phi, &T[j1+j1*TLDim], TLDim, work );
            // Apply the rotation from the right,
            //   T := T (I - phi v v') = T - (phi T v) v'
            ApplyReflector
            ( false, j1+2, 3, v, 1, phi, &T[j1*TLDim], TLDim, work );
            // Force our a priori knowledge that T is block upper-triangular
            T[(j1+2)+ j1   *TLDim] = zero;
            T[(j1+2)+(j1+1)*TLDim] = zero; 
            T[(j1+2)+(j1+2)*TLDim] = tau11;

            if( wantSchurVecs )
            {
                // Apply the rotation from the right,
                //   Q := Q (I - phi v v') = Q - (phi Q v) v'
                ApplyReflector
                ( false, n, 3, v, 1, phi, &Q[j1*QLDim], QLDim, work );
            }
        }
        else if( n1 == 2 && n2 == 1 )
        {
            const Real tau22 = T[(j1+2)+(j1+2)*TLDim];

            // Compute the Householder reflection which satisfies
            //
            //   | Q11,   q12   |^T | I, -x | = | Q11^T, r |,
            //   | q21, gamma22 |   | 0,  1 |   | q12^T, 0 |
            //
            // where x is a column vector, gamma22 is scalar, q21 is a row
            // vector, q12 is a column vector, and r is a column vector.
            //

            Real v[3];
            v[0] = -X[0];
            v[1] = -X[1];
            v[2] = scale;
            Real phi = Reflector( 3, v[0], &v[1], 1 );
            v[0] = 1; 

            if( testAccuracy )
            {
                // Perform the candidate rotation out-of-place on D
                // ------------------------------------------------ 
                // Apply the rotation from the left,
                //   D := (I - phi v v') D = D - v (phi v' D)
                ApplyReflector( true, 3, 3, v, 1, phi, D, nSum, work );
                // Apply the rotation from the right,
                //   D := D (I - phi v v') = D - (phi D v) v'
                ApplyReflector( false, 3, 3, v, 1, phi, D, nSum, work );

                // Throw an exception if the rotation would be too inaccurate.
                // As in LAPACK, rather than simply measuring the size of the 
                // bottom-left block of the rotation of D (which should ideally
                // be zero), we also take into account our knowledge that the
                // top-left entry of D should by tau22
                Real errMeasure = Max( Abs(D[1+0*3]), Abs(D[2+0*3]) );
                errMeasure = Max( errMeasure, Abs(D[0+0*3]-tau22) );
                if( errMeasure > thresh )
                    RuntimeError
                    ("Unacceptable Schur block swap: errMeasure, ",errMeasure,
                     " was greater than ",thresh);
            }

            // Perform the rotation on T
            // ------------------------- 
            // Apply the rotation from the right,
            //   T := T (I - phi v v') = T - (phi T v) v'
            ApplyReflector
            ( false, j1+3, 3, v, 1, phi, &T[j1*TLDim], TLDim, work );
            // Apply the rotation from the left,
            //   T := (I - phi v v') T = T - v (phi v' T)
            ApplyReflector
            ( true, 3, n-j1-1, v, 1, phi, &T[j1+(j1+1)*TLDim], TLDim, work );
            // Force our a priori knowledge that T is block upper-triangular
            T[ j1   +j1*TLDim] = tau22;
            T[(j1+1)+j1*TLDim] = zero;
            T[(j1+2)+j1*TLDim] = zero;

            if( wantSchurVecs )
            {
                // Apply the rotation from the right,
                //   Q := Q (I - phi v v') = Q - (phi Q v) v'
                ApplyReflector
                ( false, n, 3, v, 1, phi, &Q[j1*QLDim], QLDim, work );
            }
        }
        else
        {
            // Compute the Householder reflections which satisfy
            //
            //   Q1^T Q0^T | -X | = | R |
            //             |  I |   | 0 |
            //
            // using an inlined Householder QR factorization of [-X; I].
            //
            const Real& chi00 = X[0+0*XLDim];
            const Real& chi01 = X[0+1*XLDim];
            const Real& chi10 = X[1+0*XLDim];
            const Real& chi11 = X[1+1*XLDim];

            Real v0[3];
            v0[0] = -chi00;
            v0[1] = -chi10;
            v0[2] = scale;
            Real phi0 = Reflector( 3, v0[0], &v0[1], 1 );
            v0[0] = 1; 

            innerProd = -phi0*(chi01+v0[1]*chi11);
            Real v1[3];
            v1[0] = -innerProd*v0[1] - chi11;
            v1[1] = -innerProd*v0[2];
            v1[2] = scale;
            Real phi1 = Reflector( 3, v1[0], &v1[1], 1 ); 
            v1[0] = 1;

            if( testAccuracy )
            {
                // Perform the candidate rotation out-of-place on D
                // ------------------------------------------------ 
                // Apply the first rotation from the left,
                //   D := (I - phi0 v0 v0') D = D - v0 (phi0 v0' D)
                ApplyReflector( true, 3, 4, v0, 1, phi0, D, nSum, work );
                // Apply the first rotation from the right,
                //   D := D (I - phi0 v0 v0') = D - (phi0 D v0) v0'
                ApplyReflector( false, 4, 3, v0, 1, phi0, D, nSum, work );
                // Apply the second rotation from the left,
                //   D := (I - phi1 v1 v1') D = D - v1 (phi1 v1' D)
                ApplyReflector
                ( true, 3, 4, v1, 1, phi1, &D[1+0*nSum], nSum, work );
                // Apply the second rotation from the right,
                //   D := D (I - phi1 v1 v1') = D - (phi1 D v1) v1'
                ApplyReflector
                ( false, 4, 3, v1, 1, phi1, &D[0+1*nSum], nSum, work );

                // Throw an exception if the rotation would be too inaccurate.
                Real errMeasure = Max( Abs(D[2+0*4]), Abs(D[2+1*4]) );
                errMeasure = Max( errMeasure, Abs(D[3+0*4]) );
                errMeasure = Max( errMeasure, Abs(D[3+1*4]) );
                if( errMeasure > thresh )
                    RuntimeError
                    ("Unacceptable Schur block swap: errMeasure, ",errMeasure,
                     " was greater than ",thresh);
            }

            // Perform the rotation on T
            // ------------------------- 
            // Apply the first rotation from the left,
            //   T := (I - phi0 v0 v0') T = T - v0 (phi0 v0' T)
            ApplyReflector
            ( true, 3, n-j1, v0, 1, phi0, &T[j1+j1*TLDim], TLDim, work );
            // Apply the first rotation from the right,
            //   T := T (I - phi0 v0 v0') = T - (phi0 T v0) v0'
            ApplyReflector
            ( false, j1+4, 3, v0, 1, phi0, &T[j1*TLDim], TLDim, work );
            // Apply the second rotation from the left,
            //   T := (I - phi1 v1 v1') T = T - v1 (phi1 v1' T)
            ApplyReflector
            ( true, 3, n-j1, v1, 1, phi1, &T[(j1+1)+j1*TLDim], TLDim, work );
            // Apply the second rotation from the right,
            //   T := T (I - phi1 v1 v1') = T - (phi1 T v1) v1'
            ApplyReflector
            ( false, j1+4, 3, v1, 1, phi1, &T[(j1+1)*TLDim], TLDim, work );

            // Force our a priori knowledge that T is block upper-triangular
            T[(j1+2)+ j1   *TLDim] = zero;
            T[(j1+2)+(j1+1)*TLDim] = zero;
            T[(j1+3)+ j1   *TLDim] = zero;
            T[(j1+3)+(j1+1)*TLDim] = zero;

            if( wantSchurVecs )
            {
                // Apply the rotations from the right,
                //   Q := Q (I - phi v v') = Q - (phi Q v) v'
                ApplyReflector
                ( false, n, 3, v0, 1, phi0, &Q[j1*QLDim], QLDim, work );
                ApplyReflector
                ( false, n, 3, v1, 1, phi1, &Q[(j1+1)*QLDim], QLDim, work );
            }
        }

        // Clean up
        if( n1 == 2 )
        {
            // Force the rotated T11 into standard form
            Real c, s;
            Complex<Real> lambda0, lambda1;
            const Int j3 = j1+n2;
            const Int j4 = j3+1;
            schur::TwoByTwo
            ( T[j3+j3*TLDim], T[j3+j4*TLDim],
              T[j4+j3*TLDim], T[j4+j4*TLDim], lambda0, lambda1, c, s );
            if( j3+2 < n )
            {
                blas::Rot
                ( n-(j3+2), 
                  &T[j3+(j3+2)*TLDim], TLDim,
                  &T[j4+(j3+2)*TLDim], TLDim, c, s );
            }
            blas::Rot
            ( j3,
              &T[j3*TLDim], 1,
              &T[j4*TLDim], 1, c, s );
            if( wantSchurVecs )
            {
                blas::Rot
                ( n, 
                  &Q[j3*QLDim], 1,
                  &Q[j4*QLDim], 1, c, s );
            }
        }
        if( n2 == 2 )
        {
            // Force the rotated T22 into standard form
            Real c, s;
            Complex<Real> lambda0, lambda1;
            const Int j2 = j1+1;
            schur::TwoByTwo
            ( T[j1+j1*TLDim], T[j1+j2*TLDim],
              T[j2+j1*TLDim], T[j2+j2*TLDim], lambda0, lambda1, c, s );
            blas::Rot
            ( n-(j1+2), 
              &T[j1+(j1+2)*TLDim], TLDim,
              &T[j2+(j1+2)*TLDim], TLDim, c, s );
            blas::Rot
            ( j1,
              &T[j1*TLDim], 1,
              &T[j2*TLDim], 1, c, s );
            if( wantSchurVecs )
            {
                blas::Rot
                ( n, 
                  &Q[j1*QLDim], 1,
                  &Q[j2*QLDim], 1, c, s );
            }
        }
    }
}

} // namespace adjacent_schur

template<typename Real>
void AdjacentSchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  Real* work,
  bool testAccuracy )
{
    DEBUG_CSE
    bool wantSchurVecs = false;
    Real* Q=nullptr;
    BlasInt QLDim = 1;
    adjacent_schur::Helper
    ( wantSchurVecs, n, T, TLDim, Q, QLDim, j1, n1, n2, work, testAccuracy );
}

template<typename Real>
void AdjacentSchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim, 
  Real* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  Real* work,
  bool testAccuracy )
{
    DEBUG_CSE
    bool wantSchurVecs = true;
    adjacent_schur::Helper
    ( wantSchurVecs, n, T, TLDim, Q, QLDim, j1, n1, n2, work, testAccuracy );
}

template void AdjacentSchurExchange
( BlasInt n,
  float* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  float* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  float* T, BlasInt TLDim,
  float* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  float* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  double* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  double* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  double* T, BlasInt TLDim,
  double* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  double* work,
  bool testAccuracy );
#ifdef EL_HAVE_QUAD
template void AdjacentSchurExchange
( BlasInt n,
  Quad* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  Quad* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  Quad* T, BlasInt TLDim,
  Quad* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  Quad* work,
  bool testAccuracy );
#endif
#ifdef EL_HAVE_QD
template void AdjacentSchurExchange
( BlasInt n,
  DoubleDouble* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  DoubleDouble* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  DoubleDouble* T, BlasInt TLDim,
  DoubleDouble* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  DoubleDouble* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  QuadDouble* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  QuadDouble* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  QuadDouble* T, BlasInt TLDim,
  QuadDouble* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  QuadDouble* work,
  bool testAccuracy );
#endif
#ifdef EL_HAVE_MPC
template void AdjacentSchurExchange
( BlasInt n,
  BigFloat* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  BigFloat* work,
  bool testAccuracy );
template void AdjacentSchurExchange
( BlasInt n,
  BigFloat* T, BlasInt TLDim,
  BigFloat* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt n1,
  BlasInt n2,
  BigFloat* work,
  bool testAccuracy );
#endif

// Exchange two blocks of a real Schur decomposition
// =================================================
namespace schur_exchange {

// This following technique is an analogue of LAPACK's {s,d}trexc
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Helper
( bool wantSchurVecs,
  BlasInt n,
  Real* T, BlasInt TLDim, 
  Real* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  Real* work,
  bool testAccuracy )
{
    DEBUG_CSE
    const Real zero(0);
    if( n <= 1 )
        return;

    if( j1 > 0 && T[j1+(j1-1)*TLDim] != zero )
        --j1;
    const bool doubleFirst = ( j1 < n-1 && T[(j1+1)+j1*TLDim] != zero );
    const Int nb = ( doubleFirst ? 2 : 1 );

    if( j2 > 0 && T[j2+(j2-1)*TLDim] != zero )
        --j2;
    const bool doubleSecond = ( j1 < n-1 && T[(j1+1)+j1*TLDim] != zero );

    if( j1 == j2 )
        return; 

    bool splitDouble = false;
    if( j1 < j2 )
    {
        // We will translate down the j1 block, one swap at a time
        if( doubleFirst && !doubleSecond ) 
            --j2;
        if( !doubleFirst && doubleSecond )
            ++j2;

        for( Int j=j1; j<j2; )
        {
            if( !splitDouble )
            {
                bool doubleNext =
                  ( j+nb+1 < n && T[(j+nb+1)+(j+nb)*TLDim] != zero );
                Int nbNext = ( doubleNext ? 2 : 1 ); 
                adjacent_schur::Helper
                ( wantSchurVecs, n,
                  T, TLDim,
                  Q, QLDim,
                  j, nb, nbNext,
                  work, testAccuracy );

                j += nbNext;
                if( nb==2 && T[(j+1)+j*TLDim] == zero )
                {
                    // The 2x2 just split (this should be very rare)
                    splitDouble = true;
                }
            }
            else
            {
                bool doubleNext = ( j+3<n && T[(j+3)+(j+2)*TLDim] != zero );
                Int nbNext = ( doubleNext ? 2 : 1 );
                adjacent_schur::Helper
                ( wantSchurVecs, n,
                  T, TLDim,
                  Q, QLDim,
                  j+1, 1, nbNext,
                  work, testAccuracy );
                if( nbNext == 1 )
                {
                    adjacent_schur::Helper
                    ( wantSchurVecs, n,
                      T, TLDim,
                      Q, QLDim,
                      j, 1, nbNext,
                      work, testAccuracy );
                }
                else
                {
                    if( T[(j+2)+(j+1)*TLDim] != zero )
                    {
                        adjacent_schur::Helper
                        ( wantSchurVecs, n,
                          T, TLDim,
                          Q, QLDim,
                          j, 1, 2,
                          work, testAccuracy );
                    }
                    else
                    {
                        // The 2x2 just split (this should be very rare)
                        adjacent_schur::Helper
                        ( wantSchurVecs, n,
                          T, TLDim,
                          Q, QLDim,
                          j, 1, 1,
                          work, testAccuracy );
                        adjacent_schur::Helper
                        ( wantSchurVecs, n,
                          T, TLDim,
                          Q, QLDim,
                          j+1, 1, 1,
                          work, testAccuracy );
                    }
                }
                j += nbNext;
            }
        }
    }
    else
    {
        // We will translate up the j1 block, one swap at a time
        for( Int j=j1; j>j2; )
        {
            if( !splitDouble )
            {
                bool doubleNext = ( j>=2 && T[(j-1)+(j-2)*TLDim] != zero );
                Int nbNext = ( doubleNext ? 2 : 1 );
                adjacent_schur::Helper
                ( wantSchurVecs, n,
                  T, TLDim,
                  Q, QLDim,
                  j-nbNext, nbNext, nb,
                  work, testAccuracy );

                j -= nbNext;
                if( nb==2 && T[(j+1)+j*TLDim] == zero )
                    splitDouble = true;
            }
            else
            {
                // The 2x2 has split (this is very rare)
                bool doubleNext = ( j>=2 && T[(j-1)+(j-2)*TLDim] != zero );
                Int nbNext = ( doubleNext ? 2 : 1 );
                adjacent_schur::Helper
                ( wantSchurVecs, n,
                  T, TLDim,
                  Q, QLDim,
                  j-nbNext, nbNext, 1,
                  work, testAccuracy );
                if( nbNext == 1 )
                {
                    adjacent_schur::Helper
                    ( wantSchurVecs, n,
                      T, TLDim,
                      Q, QLDim,
                      j, nbNext, 1,
                      work, testAccuracy );
                }
                else
                {
                    if( T[j+(j-1)*TLDim] != zero )
                    {
                        adjacent_schur::Helper
                        ( wantSchurVecs, n,
                          T, TLDim,
                          Q, QLDim,
                          j-1, 2, 1,
                          work, testAccuracy );
                    }
                    else
                    {
                        // The 2x2 has just split (this is very rare)
                        adjacent_schur::Helper
                        ( wantSchurVecs, n,
                          T, TLDim,
                          Q, QLDim,
                          j, 1, 1,
                          work, testAccuracy );
                        adjacent_schur::Helper
                        ( wantSchurVecs, n,
                          T, TLDim,
                          Q, QLDim,
                          j-1, 1, 1,
                          work, testAccuracy );
                    }
                }
                j -= nbNext;
            }
        }
    }
}

template<typename Real>
void Helper
( bool wantSchurVecs,
  BlasInt n,
  Complex<Real>* T, BlasInt TLDim, 
  Complex<Real>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 )
{
    DEBUG_CSE
    if( n <= 1 || j1 == j2 )
        return;

    const Int jBeg   = ( j1<j2 ? j1 : j1-1 );
    const Int jEnd   = ( j1<j2 ? j2 : j2-1 );
    const Int stride = ( j1<j2 ? 1  :   -1 );

    for( Int j=jBeg; j!=jEnd; j+=stride )
    {
        // Save the j'th and (j+1)'th diagonal entries
        Complex<Real> tau00 = T[j+j*TLDim];
        Complex<Real> tau11 = T[(j+1)+(j+1)*TLDim];

        // Find the Givens rotation for swapping said diagonal entries
        Real c;
        Complex<Real> s;
        Givens( T[j+(j+1)*TLDim], tau11-tau00, c, s );

        // Apply the Givens rotation from the left
        if( j+2 < n )
            blas::Rot
            ( n-(j+2), &T[j    +(j+2)*TLDim], TLDim,
                       &T[(j+1)+(j+2)*TLDim], TLDim, c, s );
        // Apply the Givens rotation from the right
        if( j > 0 )
            blas::Rot( j, &T[j*TLDim], 1, &T[(j+1)*TLDim], 1, c, Conj(s) );
        if( wantSchurVecs )
            blas::Rot( n, &Q[j*QLDim], 1, &Q[(j+1)*QLDim], 1, c, Conj(s) );

        T[ j   + j   *TLDim] = tau11;
        T[(j+1)+(j+1)*TLDim] = tau00;
    }
}

} // namespace schur_exchange

template<typename Real>
void SchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2,
  Real* work,
  bool testAccuracy )
{
    DEBUG_CSE
    bool wantSchurVecs = false;
    Real* Q=nullptr;
    BlasInt QLDim = 1;
    schur_exchange::Helper
    ( wantSchurVecs, n, T, TLDim, Q, QLDim, j1, j2, work, testAccuracy );
}

template<typename Real>
void SchurExchange
( BlasInt n,
  Complex<Real>* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2 )
{
    DEBUG_CSE
    bool wantSchurVecs = false;
    Complex<Real>* Q=nullptr;
    BlasInt QLDim = 1;
    schur_exchange::Helper( wantSchurVecs, n, T, TLDim, Q, QLDim, j1, j2 );
}

template<typename Real>
void SchurExchange
( BlasInt n,
  Real* T, BlasInt TLDim, 
  Real* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  Real* work,
  bool testAccuracy )
{
    DEBUG_CSE
    bool wantSchurVecs = true;
    schur_exchange::Helper
    ( wantSchurVecs, n, T, TLDim, Q, QLDim, j1, j2, work, testAccuracy );
}

template<typename Real>
void SchurExchange
( BlasInt n,
  Complex<Real>* T, BlasInt TLDim, 
  Complex<Real>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 )
{
    DEBUG_CSE
    bool wantSchurVecs = true;
    schur_exchange::Helper
    ( wantSchurVecs, n, T, TLDim, Q, QLDim, j1, j2 );
}

template void SchurExchange
( BlasInt n,
  float* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2,
  float* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  float* T, BlasInt TLDim,
  float* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  float* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  Complex<float>* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  Complex<float>* T, BlasInt TLDim,
  Complex<float>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  double* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2,
  double* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  double* T, BlasInt TLDim,
  double* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  double* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  Complex<double>* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  Complex<double>* T, BlasInt TLDim,
  Complex<double>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 );
#ifdef EL_HAVE_QUAD
template void SchurExchange
( BlasInt n,
  Quad* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2,
  Quad* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  Quad* T, BlasInt TLDim,
  Quad* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  Quad* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  Complex<Quad>* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  Complex<Quad>* T, BlasInt TLDim,
  Complex<Quad>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 );
#endif
#ifdef EL_HAVE_QD
template void SchurExchange
( BlasInt n,
  DoubleDouble* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2,
  DoubleDouble* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  DoubleDouble* T, BlasInt TLDim,
  DoubleDouble* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  DoubleDouble* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  Complex<DoubleDouble>* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  Complex<DoubleDouble>* T, BlasInt TLDim,
  Complex<DoubleDouble>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  QuadDouble* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2,
  QuadDouble* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  QuadDouble* T, BlasInt TLDim,
  QuadDouble* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  QuadDouble* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  Complex<QuadDouble>* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  Complex<QuadDouble>* T, BlasInt TLDim,
  Complex<QuadDouble>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 );
#endif
#ifdef EL_HAVE_MPC
template void SchurExchange
( BlasInt n,
  BigFloat* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2,
  BigFloat* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  BigFloat* T, BlasInt TLDim,
  BigFloat* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2,
  BigFloat* work,
  bool testAccuracy );
template void SchurExchange
( BlasInt n,
  Complex<BigFloat>* T, BlasInt TLDim,
  BlasInt j1,
  BlasInt j2 );
template void SchurExchange
( BlasInt n,
  Complex<BigFloat>* T, BlasInt TLDim,
  Complex<BigFloat>* Q, BlasInt QLDim,
  BlasInt j1,
  BlasInt j2 );
#endif

// Compute the Schur decomposition of an upper Hessenberg matrix
// =============================================================

void HessenbergSchur
( BlasInt n,
  float* H, BlasInt ldH,
  scomplex* w,
  bool fullTriangle,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1;
    vector<float> wr( n ), wi( n );
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'), compZ='N';
        BlasInt workSize=-1;
        float workDummy;
        EL_LAPACK(shseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          0, &fakeLDim, &workDummy, &workSize, &info );

        workSize = workDummy;
        vector<float> work(workSize);
        EL_LAPACK(shseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          0, &fakeLDim, work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("shseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_FALSE;
        EL_LAPACK(slahqr)
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          &ilo, &ihi, 0, &fakeLDim, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("slahqr failed to compute all eigenvalues");
    }

    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<float>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n,
  double* H, BlasInt ldH,
  dcomplex* w,
  bool fullTriangle,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1;
    vector<double> wr( n ), wi( n );
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'), compZ='N';
        BlasInt workSize=-1;
        double workDummy;
        EL_LAPACK(dhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          0, &fakeLDim, &workDummy, &workSize, &info );

        workSize = workDummy;
        vector<double> work(workSize);
        EL_LAPACK(dhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          0, &fakeLDim, work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("dhseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_FALSE;
        EL_LAPACK(dlahqr) 
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          &ilo, &ihi, 0, &fakeLDim, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("dlahqr failed to compute all eigenvalues");
    }
    
    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<double>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n,
  scomplex* H, BlasInt ldH,
  scomplex* w,
  bool fullTriangle,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1;
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'), compZ='N';
        BlasInt workSize=-1;
        scomplex workDummy;
        EL_LAPACK(chseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
          &workDummy, &workSize, &info );

        workSize = workDummy.real();
        vector<scomplex> work(workSize);
        EL_LAPACK(chseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
          work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("chseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_FALSE;
        EL_LAPACK(clahqr) 
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, w,
          &ilo, &ihi, 0, &fakeLDim, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("clahqr failed to compute all eigenvalues");
    }
}

void HessenbergSchur
( BlasInt n,
  dcomplex* H, BlasInt ldH,
  dcomplex* w,
  bool fullTriangle,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    BlasInt fakeLDim=1;
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'), compZ='N';
        BlasInt workSize=-1;
        dcomplex workDummy;
        EL_LAPACK(zhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
          &workDummy, &workSize, &info );

        workSize = workDummy.real();
        vector<dcomplex> work(workSize);
        EL_LAPACK(zhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, 0, &fakeLDim, 
          work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("zhseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_FALSE;
        EL_LAPACK(zlahqr) 
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, w,
          &ilo, &ihi, 0, &fakeLDim, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("zlahqr failed to compute all eigenvalues");
    }
}

void HessenbergSchur
( BlasInt n,
  float* H, BlasInt ldH,
  scomplex* w,
  float* Q, BlasInt ldQ, 
  bool fullTriangle,
  bool multiplyQ,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    vector<float> wr( n ), wi( n );
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'),
                   compZ=(multiplyQ ? 'V' : 'I');
        BlasInt workSize=-1;
        float workDummy;
        EL_LAPACK(shseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
          &workDummy, &workSize, &info );

        workSize = workDummy;
        vector<float> work(workSize);
        EL_LAPACK(shseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
          work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("shseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_TRUE;
        EL_LAPACK(slahqr)
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          &ilo, &ihi, Q, &ldQ, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("slahqr failed to compute all eigenvalues");
    }

    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<float>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n,
  double* H, BlasInt ldH,
  dcomplex* w,
  double* Q, BlasInt ldQ, 
  bool fullTriangle,
  bool multiplyQ,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    vector<double> wr( n ), wi( n );
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'),
                   compZ=(multiplyQ ? 'V' : 'I');
        BlasInt workSize=-1;
        double workDummy;
        EL_LAPACK(dhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
          &workDummy, &workSize, &info );

        workSize = workDummy;
        vector<double> work(workSize);
        EL_LAPACK(dhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(), Q, &ldQ,
          work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("dhseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_TRUE;
        EL_LAPACK(dlahqr)
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, wr.data(), wi.data(),
          &ilo, &ihi, Q, &ldQ, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("dlahqr failed to compute all eigenvalues");
    }
    
    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<double>(wr[i],wi[i]);
}

void HessenbergSchur
( BlasInt n,
  scomplex* H, BlasInt ldH,
  scomplex* w,
  scomplex* Q, BlasInt ldQ,
  bool fullTriangle,
  bool multiplyQ,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'),
                   compZ=(multiplyQ ? 'V' : 'I');
        BlasInt workSize=-1;
        scomplex workDummy;
        EL_LAPACK(chseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ,
          &workDummy, &workSize, &info );

        workSize = workDummy.real();
        vector<scomplex> work(workSize);
        EL_LAPACK(chseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ,
          work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("chseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_TRUE;
        EL_LAPACK(clahqr)
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, w,
          &ilo, &ihi, Q, &ldQ, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("clahqr failed to compute all eigenvalues");
    }
}

void HessenbergSchur
( BlasInt n,
  dcomplex* H, BlasInt ldH,
  dcomplex* w,
  dcomplex* Q, BlasInt ldQ,
  bool fullTriangle,
  bool multiplyQ,
  bool useAED )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    BlasInt info;
    BlasInt ilo=1, ihi=n;
    if( useAED )
    {
        const char job=(fullTriangle ? 'S' : 'E'),
                   compZ=(multiplyQ ? 'V' : 'I');
        BlasInt workSize=-1;
        dcomplex workDummy;
        EL_LAPACK(zhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ,
          &workDummy, &workSize, &info );

        workSize = workDummy.real();
        vector<dcomplex> work(workSize);
        EL_LAPACK(zhseqr)
        ( &job, &compZ, &n, &ilo, &ihi, H, &ldH, w, Q, &ldQ,
          work.data(), &workSize, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("zhseqr failed to compute all eigenvalues");
    }
    else
    {
        FortranLogical wantT=(fullTriangle ? FORTRAN_TRUE : FORTRAN_FALSE),
                       wantZ=FORTRAN_TRUE;
        EL_LAPACK(zlahqr)
        ( &wantT, &wantZ, &n, &ilo, &ihi, H, &ldH, w,
          &ilo, &ihi, Q, &ldQ, &info );
        if( info < 0 )
            RuntimeError("Argument ",-info," had an illegal value");
        else if( info > 0 )
            RuntimeError("zlahqr failed to compute all eigenvalues");
    }
}

// Compute eigenvalues/pairs of an upper Hessenberg matrix
// =======================================================

void HessenbergEig( BlasInt n, float* H, BlasInt ldH, scomplex* w )
{
    DEBUG_CSE
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( BlasInt n, double* H, BlasInt ldH, dcomplex* w )
{
    DEBUG_CSE
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( BlasInt n, scomplex* H, BlasInt ldH, scomplex* w )
{
    DEBUG_CSE
    HessenbergSchur( n, H, ldH, w, false );
}

void HessenbergEig( BlasInt n, dcomplex* H, BlasInt ldH, dcomplex* w )
{
    DEBUG_CSE
    HessenbergSchur( n, H, ldH, w, false );
}

// TODO: Compute eigenpairs

// Compute the Schur decomposition of a square matrix
// ==================================================

// TODO: Add timing of components

void Schur
( BlasInt n,
  float* A, BlasInt ldA,
  scomplex* w,
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    vector<float> tau( n );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    vector<float> wr( n ), wi( n );
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<float> work( workSize );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<float>(wr[i],wi[i]);
}

void Schur
( BlasInt n,
  double* A, BlasInt ldA,
  dcomplex* w,
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    vector<double> tau( n );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    vector<double> wr( n ), wi( n );
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<double> work( workSize );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), 0, &fakeLDim,
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<double>(wr[i],wi[i]);
}

void Schur
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* w,
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    vector<scomplex> tau( n );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<scomplex> work( workSize );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Schur
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* w,
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    vector<dcomplex> tau( n );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='N';
    BlasInt fakeLDim=1, negOne=-1;
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<dcomplex> work( workSize );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Compute the eigenvalues
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, 0, &fakeLDim, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("zhseqr's failed to compute all eigenvalues");
}

void Schur
( BlasInt n,
  float* A, BlasInt ldA,
  scomplex* w,
  float* Q, BlasInt ldQ, 
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    float workDummy;
    vector<float> tau( n );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    vector<float> wr( n ), wi( n );
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<float> work( workSize );
    EL_LAPACK(sgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(sorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(shseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("shseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<float>(wr[i],wi[i]);
}

void Schur
( BlasInt n,
  double* A, BlasInt ldA,
  dcomplex* w,
  double* Q, BlasInt ldQ, 
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    double workDummy;
    vector<double> tau( n );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy;

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    vector<double> wr( n ), wi( n );
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy), workSize );

    // Reduce to Hessenberg form
    vector<double> work( workSize );
    EL_LAPACK(dgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(dorghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(dhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, wr.data(), wi.data(), Q, &ldQ, 
      work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("dhseqr's failed to compute all eigenvalues");

    // Return the complex eigenvalues
    for( BlasInt i=0; i<n; ++i )
        w[i] = Complex<double>(wr[i],wi[i]);
}

void Schur
( BlasInt n,
  scomplex* A, BlasInt ldA,
  scomplex* w,
  scomplex* Q, BlasInt ldQ, 
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    scomplex workDummy;
    vector<scomplex> tau( n );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<scomplex> work( workSize );
    EL_LAPACK(cgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(cunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(chseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, work.data(), &workSize, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

void Schur
( BlasInt n,
  dcomplex* A, BlasInt ldA,
  dcomplex* w,
  dcomplex* Q, BlasInt ldQ, 
  bool fullTriangle,
  bool time )
{
    DEBUG_CSE
    if( n == 0 )
        return;

    // Query the reduction to Hessenberg form workspace size
    BlasInt ilo=1, ihi=n, workSize=-1, info;
    dcomplex workDummy;
    vector<dcomplex> tau( n );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), &workDummy, &workSize, &info );
    workSize = workDummy.real();

    // Query the explicit Q formation workspace
    BlasInt negOne=-1; 
    EL_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Query the QR algorithm workspace size
    const char job = ( fullTriangle ? 'S' : 'E' ), compZ='V';
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, 
      &workDummy, &negOne, &info );
    workSize = Max( BlasInt(workDummy.real()), workSize );

    // Reduce to Hessenberg form
    vector<dcomplex> work( workSize );
    EL_LAPACK(zgehrd)
    ( &n, &ilo, &ihi, A, &ldA, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of reduction had an illegal value");

    // Copy the Householder vectors over
    for( BlasInt j=0; j<n; ++j )
        MemCopy( &Q[j*ldQ], &A[j*ldA], n );

    // Form the orthogonal matrix in place
    EL_LAPACK(zunghr)
    ( &n, &ilo, &ihi, Q, &ldQ, tau.data(), work.data(), &workSize, &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of formation had an illegal value");

    // Compute the Schur decomposition
    EL_LAPACK(zhseqr)
    ( &job, &compZ, &n, &ilo, &ihi, A, &ldA, w, Q, &ldQ, work.data(), &workSize, 
      &info );
    if( info < 0 )
        RuntimeError("Argument ",-info," of QR alg had an illegal value");
    else if( info > 0 )
        RuntimeError("chseqr's failed to compute all eigenvalues");
}

} // namespace lapack
} // namespace El

#include "./lapack/TriangEig.hpp"
#include "./lapack/Schur.hpp"
