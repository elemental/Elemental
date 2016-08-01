/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(strmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const float* A, const BlasInt* ALDim, float* x, const BlasInt* incx );
void EL_BLAS(dtrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const double* A, const BlasInt* ALDim, double* x, const BlasInt* incx );
void EL_BLAS(ctrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const scomplex* A, const BlasInt* ALDim,
        scomplex* x, const BlasInt* incx );
void EL_BLAS(ztrmv)
( const char* uplo, const char* trans, const char* diag, const BlasInt* m,
  const dcomplex* A, const BlasInt* ALDim,
        dcomplex* x, const BlasInt* incx );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Trmv
( char uplo, char trans, char diag,
  BlasInt m,
  const T* A, BlasInt ALDim,
        T* x, BlasInt incx )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    const bool lower = ( std::toupper(uplo) == 'L' );
    const bool conj = ( std::toupper(trans) == 'C' );
    const bool unitDiag = ( std::toupper(diag) == 'U' );
    const T zero(0);

    T gamma, delta;
    if( lower )
    {
        if( std::toupper(trans) == 'N' )
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                gamma = x[j*incx]; 
                if( gamma != zero )
                {
                    for( BlasInt i=m-1; i>j; --i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= gamma;
                        x[i*incx] += delta;
                    }
                    if( !unitDiag )
                        x[j*incx] *= A[j+j*ALDim];
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<m; ++j )
            {
                gamma = x[j*incx]; 
                if( conj )
                {
                    if( !unitDiag )
                    {
                        Conj( A[j+j*ALDim], delta );
                        gamma *= delta;
                    }
                    for( BlasInt i=j+1; i<m; ++i )
                    {
                        Conj( A[i+j*ALDim], delta );
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                else
                {
                    if( !unitDiag )
                        gamma *= A[j+j*ALDim];
                    for( BlasInt i=j+1; i<m; ++i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                x[j*incx] = gamma;
            }
        }
    }
    else
    {
        if( std::toupper(trans) == 'N' )
        {
            for( BlasInt j=0; j<m; ++j )
            {
                gamma = x[j*incx]; 
                if( gamma != zero )
                {
                    for( BlasInt i=0; i<j; ++i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= gamma;
                        x[i*incx] += delta;
                    }
                    if( !unitDiag )
                        x[j*incx] *= A[j+j*ALDim];
                }
            }
        }
        else
        {
            for( BlasInt j=m-1; j>=0; --j )
            {
                gamma = x[j*incx]; 
                if( conj )
                {
                    if( !unitDiag )
                    {
                        Conj( A[j+j*ALDim], delta );
                        gamma *= delta;
                    }
                    for( BlasInt i=j-1; i>=0; --i )
                    {
                        Conj( A[i+j*ALDim], delta );
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                else
                {
                    if( !unitDiag )
                        gamma *= A[j+j*ALDim];
                    for( BlasInt i=j-1; i>=0; --i )
                    {
                        delta = A[i+j*ALDim];
                        delta *= x[i*incx];
                        gamma += delta;
                    }
                }
                x[j*incx] = gamma;
            }
        }
    }
}
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Int* A, BlasInt ALDim,
        Int* x, BlasInt incx );
#ifdef EL_HAVE_QD
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const DoubleDouble* A, BlasInt ALDim,
        DoubleDouble* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const QuadDouble* A, BlasInt ALDim,
        QuadDouble* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
        Complex<DoubleDouble>* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<QuadDouble>* A, BlasInt ALDim,
        Complex<QuadDouble>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_QUAD
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Quad* A, BlasInt ALDim,
        Quad* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<Quad>* A, BlasInt ALDim,
        Complex<Quad>* x, BlasInt incx );
#endif
#ifdef EL_HAVE_MPC
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const BigInt* A, BlasInt ALDim,
        BigInt* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const BigFloat* A, BlasInt ALDim,
        BigFloat* x, BlasInt incx );
template void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const Complex<BigFloat>* A, BlasInt ALDim,
        Complex<BigFloat>* x, BlasInt incx );
#endif

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt ALDim,
        float* x, BlasInt incx )
{ EL_BLAS(strmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt ALDim,
        double* x, BlasInt incx )
{ EL_BLAS(dtrmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt ALDim,
        scomplex* x, BlasInt incx )
{ EL_BLAS(ctrmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* x, BlasInt incx )
{ EL_BLAS(ztrmv)( &uplo, &trans, &diag, &m, A, &ALDim, x, &incx ); }

} // namespace blas
} // namespace El
