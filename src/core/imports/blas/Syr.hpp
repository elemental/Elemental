/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(cher)
( const char* uplo, const BlasInt* m,
  const float* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* A, const BlasInt* ALDim );
void EL_BLAS(zher)
( const char* uplo, const BlasInt* m,
  const double* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* A, const BlasInt* ALDim );

void EL_BLAS(ssyr)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* x, const BlasInt* incx,
                            float* A, const BlasInt* ALDim  );
void EL_BLAS(dsyr)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* x, const BlasInt* incx,
                             double* A, const BlasInt* ALDim  );
// 'csyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csyr)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
        scomplex* A, const BlasInt* ALDim  );
// 'zsyr' is an auxilliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsyr)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
        dcomplex* A, const BlasInt* ALDim  );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Her
( char uplo, BlasInt m,
  const Base<T>& alpha,
  const T* x, BlasInt incx,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta;
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( x[j*incx], gamma );
            gamma *= alpha;
            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*Conj(x[j*incx]);
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( x[j*incx], gamma );
            gamma *= alpha;
            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*Conj(x[j*incx]);
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
}
template void Her
( char uplo, BlasInt m,
  const Int& alpha,
  const Int* x, BlasInt incx, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Her
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx,
        DoubleDouble* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx,
        QuadDouble* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const Complex<DoubleDouble>* x, BlasInt incx,
        Complex<DoubleDouble>* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const Complex<QuadDouble>* x, BlasInt incx,
        Complex<QuadDouble>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Her
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx,
        Quad* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Complex<Quad>* x, BlasInt incx, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Her
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx,
        BigInt* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx,
        BigFloat* A, BlasInt ALDim );
template void Her
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const Complex<BigFloat>* x, BlasInt incx,
        Complex<BigFloat>* A, BlasInt ALDim );
#endif

void Her
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
        float* A, BlasInt ALDim )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Her
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Her
( char uplo, BlasInt m,
  const float& alpha,
  const scomplex* x, BlasInt incx,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cher)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Her
( char uplo, BlasInt m,
  const double& alpha,
  const dcomplex* x, BlasInt incx,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zher)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

template<typename T>
void Syr
( char uplo, BlasInt m,
  const T& alpha,
  const T* x, BlasInt incx,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta;
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = x[j*incx];
            gamma *= alpha;
            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*x[j*incx];
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = x[j*incx];
            gamma *= alpha;
            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*x[j*incx];
                delta = x[i*incx];
                delta *= gamma;
                A[i+j*ALDim] += delta;
            }
        }
    }
}
template void Syr
( char uplo, BlasInt m,
  const Int& alpha,
  const Int* x, BlasInt incx, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Syr
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx,
        DoubleDouble* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx,
        QuadDouble* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* x, BlasInt incx,
        Complex<DoubleDouble>* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* x, BlasInt incx,
        Complex<QuadDouble>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syr
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx,
        Quad* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Syr
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx,
        BigInt* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx,
        BigFloat* A, BlasInt ALDim );
template void Syr
( char uplo, BlasInt m, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* x, BlasInt incx,
        Complex<BigFloat>* A, BlasInt ALDim );
#endif

void Syr
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
        float* A, BlasInt ALDim  )
{ EL_BLAS(ssyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Syr
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); }

void Syr
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx,
        scomplex* A, BlasInt ALDim )
{
    // Recall that 'csyr' is an LAPACK auxiliary routine
    EL_LAPACK(csyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); 
}

void Syr
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx,
        dcomplex* A, BlasInt ALDim )
{
    // Recall that 'zsyr' is an LAPACK auxiliary routine
    EL_LAPACK(zsyr)( &uplo, &m, &alpha, x, &incx, A, &ALDim ); 
}

} // namespace blas
} // namespace El
