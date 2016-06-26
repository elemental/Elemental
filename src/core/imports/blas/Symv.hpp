/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(chemv)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
void EL_BLAS(zhemv)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

void EL_BLAS(ssymv)
( const char* uplo, const BlasInt* m,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* x, const BlasInt* incx,
  const float* beta,        float* y, const BlasInt* incy );
void EL_BLAS(dsymv)
( const char* uplo, const BlasInt* m,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* x, const BlasInt* incx,
  const double* beta,        double* y, const BlasInt* incy );
// 'csymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(csymv)
( const char* uplo, const BlasInt* m,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* x, const BlasInt* incx,
  const scomplex* beta,
        scomplex* y, const BlasInt* incy );
// 'zsymv' is an auxiliary LAPACK routine, but we will treat it as BLAS
void EL_LAPACK(zsymv)
( const char* uplo, const BlasInt* m,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* x, const BlasInt* incx,
  const dcomplex* beta,
        dcomplex* y, const BlasInt* incy );

} // extern "C"

namespace El {
namespace blas {

// TODO: Introduce some sort of blocking
template<typename T>
void Hemv
( char uplo, BlasInt m,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

    // y := beta*y
    if( beta == T(0) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] *= beta;
    }

    // Pre-scale x to avoid redundant computation
    vector<T> xAlpha(m);
    for( Int i=0; i<m; ++i )
    {
        xAlpha[i] = x[i*incx];
        xAlpha[i] *= alpha;
    }

    T gamma;
    if( std::toupper(uplo) == 'L' )
    {
        // Multiply with the lower triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=j; i<m; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the adjoint of the strictly lower triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=j+1; i<m; ++i )
            {
                // y[j*incy] += alpha*Conj(A[i+j*ALDim])*x[i*incx];
                Conj( A[i+j*ALDim], gamma );
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
    else
    {
        // Multiply with the upper triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=0; i<=j; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the adjoint of the strictly upper triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=0; i<j; ++i )
            {
                // y[j*incy] += alpha*Conj(A[i+j*ALDim])*x[i*incx];
                Conj( A[i+j*ALDim], gamma );
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
}

template void Hemv
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* x, BlasInt incx, 
  const Int& beta,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Hemv
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble& beta,
        DoubleDouble* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble& beta,
        QuadDouble* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim, 
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim, 
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Hemv
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* x, BlasInt incx, 
  const Quad& beta,
        Quad* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>& beta,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Hemv
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* x, BlasInt incx, 
  const BigInt& beta,
        BigInt* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* x, BlasInt incx, 
  const BigFloat& beta,
        BigFloat* y, BlasInt incy );
template void Hemv
( char uplo, BlasInt m, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim, 
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* y, BlasInt incy );
#endif

void Hemv
( char uplo, BlasInt m,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy )
{ EL_BLAS(chemv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Hemv
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy )
{ EL_BLAS(zhemv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

// TODO: Introduce some sort of blocking
template<typename T>
void Symv
( char uplo, BlasInt m,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

    // y := beta*y
    if( beta == T(0) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt i=0; i<m; ++i )
            y[i*incy] *= beta;
    }

    // Pre-scale x to avoid redundant computation
    vector<T> xAlpha(m);
    for( Int i=0; i<m; ++i )
    {
        xAlpha[i] = x[i*incx];
        xAlpha[i] *= alpha;
    }

    T gamma;
    if( std::toupper(uplo) == 'L' )
    {
        // Multiply with the lower triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=j; i<m; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the transpose of the strictly lower triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=j+1; i<m; ++i )
            {
                // y[j*incy] += alpha*A[i+j*ALDim]*x[i*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
    else
    {
        // Multiply with the upper triangle
        for( BlasInt j=0; j<m; ++j )
        {
            for( BlasInt i=0; i<=j; ++i )
            {
                // y[i*incy] += alpha*A[i+j*ALDim]*x[j*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[j];
                y[i*incy] += gamma;
            }
        }

        // Multiply with the transpose of the strictly upper triangle
        for( BlasInt j=0; j<m; ++j ) 
        {
            for( BlasInt i=0; i<j; ++i )
            {
                // y[j*incy] += alpha*A[i+j*ALDim]*x[i*incx];
                gamma = A[i+j*ALDim];
                gamma *= xAlpha[i];
                y[j*incy] += gamma;
            }
        }
    }
}

template void Symv
( char uplo, BlasInt m, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* x, BlasInt incx, 
  const Int& beta,
        Int* y, BlasInt incy );
#ifdef EL_HAVE_QD
template void Symv
( char uplo, BlasInt m, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* x, BlasInt incx, 
  const DoubleDouble& beta,
        DoubleDouble* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* x, BlasInt incx, 
  const QuadDouble& beta,
        QuadDouble* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim, 
  const Complex<DoubleDouble>* x, BlasInt incx, 
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim, 
  const Complex<QuadDouble>* x, BlasInt incx, 
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_QUAD
template void Symv
( char uplo, BlasInt m, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* x, BlasInt incx, 
  const Quad& beta,
        Quad* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* x, BlasInt incx, 
  const Complex<Quad>& beta,
        Complex<Quad>* y, BlasInt incy );
#endif
#ifdef EL_HAVE_MPC
template void Symv
( char uplo, BlasInt m, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* x, BlasInt incx, 
  const BigInt& beta,
        BigInt* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* x, BlasInt incx, 
  const BigFloat& beta,
        BigFloat* y, BlasInt incy );
template void Symv
( char uplo, BlasInt m, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim, 
  const Complex<BigFloat>* x, BlasInt incx, 
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* y, BlasInt incy );
#endif

void Symv
( char uplo, BlasInt m,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy )
{ EL_BLAS(ssymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, BlasInt m,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy )
{ EL_BLAS(dsymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy ); }

void Symv
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy )
{
    // Recall that 'csymv' is an LAPACK auxiliary routine
    EL_LAPACK(csymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

void Symv
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy )
{
    // Recall that 'zsymv' is an LAPACK auxiliary routine
    EL_LAPACK(zsymv)( &uplo, &m, &alpha, A, &ALDim, x, &incx, &beta, y, &incy );
}

} // namespace blas
} // namespace El
