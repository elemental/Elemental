/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(cherk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const float* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zherk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const double* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(ssyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* ALDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dsyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* ALDim,
  const double* beta,        double* C, const BlasInt* CLDim );
void EL_BLAS(csyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zsyrk)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Base<T>& alpha,
  const T* A, BlasInt ALDim,
  const Base<T>& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( beta == Base<T>(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] = 0;
    }
    else if( beta != Base<T>(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] *= beta;
    }

    const bool normal = ( std::toupper(trans) == 'N' );
    const bool lower = ( std::toupper(uplo) == 'L' );

    T gamma, delta;
    if( normal )
    {
        // C := alpha A A^H + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[j+l*ALDim], delta );
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[j+l*ALDim], delta );
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^H A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Int& alpha,
  const Int* A, BlasInt ALDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const DoubleDouble& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
  const DoubleDouble& beta,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const QuadDouble& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim,
  const QuadDouble& beta,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Quad& alpha,
  const Complex<Quad>* A, BlasInt ALDim,
  const Quad& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
template void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const BigFloat& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim,
  const BigFloat& beta,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    const char transFixed = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(ssyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    const char transFixed = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(dsyrk)
    ( &uplo, &transFixed, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const scomplex* A, BlasInt ALDim,
  const float& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(cherk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const dcomplex* A, BlasInt ALDim,
  const double& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zherk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

template<typename T>
void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] = 0;
    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<n; ++i )
                C[i+j*CLDim] *= beta;
    }

    const bool normal = ( std::toupper(trans) == 'N' );
    const bool lower = ( std::toupper(uplo) == 'L' );

    T gamma, delta;
    if( normal )
    {
        // C := alpha A A^T + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[j+l*ALDim];
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[j+l*ALDim];
                        delta *= A[i+l*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^T A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[l+i*ALDim];
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<=j; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[l+i*ALDim];
                        delta *= A[l+j*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim, 
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim, 
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
template void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim, 
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(csyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

void Syrk
( char uplo, char trans, BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zsyrk)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, &beta, C, &CLDim );
}

} // namespace blas
} // namespace El
