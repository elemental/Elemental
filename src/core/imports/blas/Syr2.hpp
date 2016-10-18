/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(cher2)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* x, const BlasInt* incx,
  const scomplex* y, const BlasInt* incy,
        scomplex* A, const BlasInt* ALDim );
void EL_BLAS(zher2)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* y, const BlasInt* incy,
        dcomplex* A, const BlasInt* ALDim );

void EL_BLAS(ssyr2)
( const char* uplo, const BlasInt* m,
  const float* alpha,
  const float* x, const BlasInt* incx,
  const float* y, const BlasInt* incy,
        float* A, const BlasInt* ALDim  );
void EL_BLAS(dsyr2)
( const char* uplo, const BlasInt* m,
  const double* alpha,
  const double* x, const BlasInt* incx,
  const double* y, const BlasInt* incy,
        double* A, const BlasInt* ALDim  );

void EL_BLAS(csyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zsyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Her2
( char uplo, BlasInt m, 
  const T& alpha,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta, phi, psi, alphaConj;
    Conj( alpha, alphaConj );
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( y[j*incy], gamma );
            gamma *= alpha;
            Conj( x[j*incx], delta );
            delta *= alphaConj;
            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] += alpha*      x[i*incx]*Conj(y[j*incy]) + 
                //                 Conj(alpha)*y[i*incy]*Conj(x[j*incx]);
                phi = x[i*incx]; 
                phi *= gamma;
                psi = y[i*incy];
                psi *= delta;
                phi += psi;
                A[i+j*ALDim] += phi;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            Conj( y[j*incy], gamma );
            gamma *= alpha;
            Conj( x[j*incx], delta );
            delta *= alphaConj;
            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*      x[i*incx]*Conj(y[j*incy]) + 
                //                 Conj(alpha)*y[i*incy]*Conj(x[j*incx]);
                phi = x[i*incx]; 
                phi *= gamma;
                psi = y[i*incy];
                psi *= delta;
                phi += psi;
                A[i+j*ALDim] += phi;
            }
        }
    }
}
template void Her2
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* x, BlasInt incx, 
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Her2
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>* y, BlasInt incy, 
        Complex<DoubleDouble>* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>* y, BlasInt incy, 
        Complex<QuadDouble>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Her2
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Her2
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
template void Her2
( char uplo, BlasInt m, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>* y, BlasInt incy, 
        Complex<BigFloat>* A, BlasInt ALDim );
#endif

void Her2
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Her2
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Her2
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{ EL_BLAS(cher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Her2
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{ EL_BLAS(zher2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

template<typename T>
void Syr2
( char uplo, BlasInt m, 
  const T& alpha,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    T gamma, delta, phi, psi;
    if( std::toupper(uplo) == 'L' )
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = y[j*incy];
            gamma *= alpha;

            delta = x[j*incx];
            delta *= alpha;

            for( BlasInt i=j; i<m; ++i )
            {
                // A[i+j*ALDim] +=
                //   alpha*(x[i*incx]*y[j*incy]+y[i*incy]*x[j*incx])
                phi = x[i*incx]; 
                phi *= gamma;

                psi = y[i*incy];
                psi *= delta;

                phi += psi;

                A[i+j*ALDim] += phi;
            }
        }
    }
    else
    {
        for( BlasInt j=0; j<m; ++j )
        {
            gamma = y[j*incy];
            gamma *= alpha;

            delta = x[j*incx]; 
            delta *= alpha;

            for( BlasInt i=0; i<=j; ++i )
            {
                // A[i+j*ALDim] += alpha*x[i*incx]*y[j*incy] + 
                //                 alpha*y[i*incy]*x[j*incx];
                phi = x[i*incx]; 
                phi *= gamma;

                psi = y[i*incy];
                psi *= delta;

                phi += psi;

                A[i+j*ALDim] += phi;
            }
        }
    }
}
template void Syr2
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* x, BlasInt incx,
  const Int* y, BlasInt incy, 
        Int* A, BlasInt ALDim );
#ifdef EL_HAVE_QD
template void Syr2
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble* y, BlasInt incy, 
        DoubleDouble* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble* y, BlasInt incy, 
        QuadDouble* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>* y, BlasInt incy, 
        Complex<DoubleDouble>* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>* y, BlasInt incy, 
        Complex<QuadDouble>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syr2
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* x, BlasInt incx, 
  const Quad* y, BlasInt incy, 
        Quad* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>* y, BlasInt incy, 
        Complex<Quad>* A, BlasInt ALDim );
#endif
#ifdef EL_HAVE_MPC
template void Syr2
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* x, BlasInt incx, 
  const BigInt* y, BlasInt incy, 
        BigInt* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* x, BlasInt incx, 
  const BigFloat* y, BlasInt incy, 
        BigFloat* A, BlasInt ALDim );
template void Syr2
( char uplo, BlasInt m, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>* y, BlasInt incy, 
        Complex<BigFloat>* A, BlasInt ALDim );
#endif

void Syr2
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim )
{ EL_BLAS(ssyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Syr2
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim )
{ EL_BLAS(dsyr2)( &uplo, &m, &alpha, x, &incx, y, &incy, A, &ALDim ); }

void Syr2
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim )
{
    // csyr2 doesn't exist, so we route through csyr2k. However, csyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const BlasInt k = 1;
    const scomplex beta = 1.f;
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &ALDim );
}

void Syr2
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim )
{
    // zsyr2 doesn't exist, so we route through zsyr2k. However, zsyr2k expects 
    // contiguous access of 'x', so we treat x and y as a row vectors where 
    // their leading dimensions are 'incx' and 'incy'. Thus we must perform 
    // A += x' y + y' x
    const char trans = 'T';
    const BlasInt k = 1;
    const dcomplex beta = 1.;
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &m, &k, &alpha, x, &incx, y, &incy, &beta, A, &ALDim );
}

} // namespace blas
} // namespace El
