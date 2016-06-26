/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(sgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* x, const BlasInt* incx,
  const float* beta,        float* y, const BlasInt* incy );
void EL_BLAS(dgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* x, const BlasInt* incx,
  const double* beta,        double* y, const BlasInt* incy );
void EL_BLAS(cgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zgemv)
( const char* trans, const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Gemv
( char trans, BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    // TODO: Special-case alpha=0, alpha=1, and alpha=-1?
    T gamma, delta;
    if( std::toupper(trans) == 'N' )
    {
        if( m > 0 && n == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<m; ++i )
                y[i*incy] = 0;
            return;
        }

        Scal( m, beta, y, incy );
        for( BlasInt j=0; j<n; ++j )
        {
            gamma = x[j*incx];
            gamma *= alpha;
            for( BlasInt i=0; i<m; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                delta = A[i+j*ALDim];
                delta *= gamma;
                y[i*incy] += delta;
            }
        }
    }
    else if( std::toupper(trans) == 'T' )
    {
        if( n > 0 && m == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<n; ++i )
                y[i*incy] = 0;
            return;
        }
        Scal( n, beta, y, incy );

        // Prescale x to avoid a factor of two more work than necessary
        vector<T> xAlpha(m);
        for( Int j=0; j<m; ++j )
        {
            xAlpha[j] = x[j*incx];
            xAlpha[j] *= alpha;
        }

        for( BlasInt i=0; i<n; ++i )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                // y[i*incy] += alpha*A[j+i*ALDim]*x[j*incx];
                gamma = A[j+i*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }
    }
    else
    {
        if( n > 0 && m == 0 && beta == T(0) )
        {
            for( BlasInt i=0; i<n; ++i )
                y[i*incy] = 0;
            return;
        }
        Scal( n, beta, y, incy );

        // Prescale x to avoid a factor of two more work than necessary
        vector<T> xAlpha(m);
        for( Int j=0; j<m; ++j )
        {
            xAlpha[j] = x[j*incx];
            xAlpha[j] *= alpha;
        }

        for( BlasInt i=0; i<n; ++i )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                // y[i*incy] += alpha*Conj(A[j+i*ALDim])*x[j*incx];
                Conj( A[j+i*ALDim], gamma );
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }
    }
}
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Int& alpha,
  const Int* A, BlasInt ALDim,
  const Int* x, BlasInt incx, 
  const Int& beta,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble& beta,
        DoubleDouble* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble& beta,
        QuadDouble* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim,
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
  const Quad* x, BlasInt incx, 
  const Quad& beta,
        Quad* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>& beta,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
  const BigInt* x, BlasInt incx, 
  const BigInt& beta,
        BigInt* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
  const BigFloat* x, BlasInt incx, 
  const BigFloat& beta,
        BigFloat* y, BlasInt incy );
template void Gemv
( char trans, BlasInt m, BlasInt n, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim,
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* y, BlasInt incy );
#endif

void Gemv
( char trans, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(sgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(dgemv)
    ( &fixedTrans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

void Gemv
( char trans, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy )
{ EL_BLAS(cgemv)
  ( &trans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Gemv
( char trans, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zgemv)
  ( &trans, &m, &n, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

} // namespace blas
} // namespace El
