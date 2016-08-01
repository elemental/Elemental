/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(cher2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const float* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zher2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const double* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(ssyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* B, const BlasInt* BLDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dsyr2k)
( const char* uplo, const char* trans,
  const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* B, const BlasInt* BLDim,
  const double* beta,        double* C, const BlasInt* CLDim );
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
void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* B, BlasInt BLDim,
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

    const T alphaConj = Conj(alpha);
    const bool normal = ( std::toupper(trans) == 'N' );
    const bool lower = ( std::toupper(uplo) == 'L' );

    T gamma, delta, phi;
    if( normal )
    {
        // C := alpha A B^H + Conj(alpha) B A^H + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( B[j+l*BLDim], phi );
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        Conj( A[j+l*ALDim], phi );
                        phi *= B[i+l*BLDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
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
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( B[j+l*BLDim], phi );
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        Conj( A[j+l*ALDim], phi );
                        phi *= B[i+l*BLDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^H B + Conj(alpha) B^H A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], phi );
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        Conj( B[l+i*BLDim], phi );
                        phi *= A[l+j*ALDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
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
                    delta = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], phi );
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        Conj( B[l+i*BLDim], phi );
                        phi *= A[l+j*ALDim];
                        delta += phi;
                    }
                    gamma *= alpha;
                    delta *= alphaConj;
                    gamma += delta;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Int& alpha,
  const Int* A, BlasInt ALDim,
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
  const Complex<DoubleDouble>* B, BlasInt BLDim,
  const DoubleDouble& beta,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim,
  const Complex<QuadDouble>* B, BlasInt BLDim,
  const QuadDouble& beta,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim,
  const Complex<Quad>* B, BlasInt BLDim,
  const Quad& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
template void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim,
  const Complex<BigFloat>* B, BlasInt BLDim,
  const BigFloat& beta,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(ssyr2k)
    ( &uplo, &transFixed, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    const char transFixed = ( trans == 'C' ? 'T' : trans );
    EL_BLAS(dsyr2k)
    ( &uplo, &transFixed, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex* B, BlasInt BLDim,
  const float& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(cher2k)
    ( &uplo, &trans, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
  const dcomplex* B, BlasInt BLDim,
  const double& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zher2k)
    ( &uplo, &trans, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

template<typename T>
void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
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

    T gamma, phi;
    if( normal )
    {
        // C := alpha A B^T + alpha B A^T + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        phi = B[j+l*BLDim];
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        phi = A[j+l*ALDim];
                        phi *= B[i+l*BLDim];
                        gamma += phi;
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
                        phi = B[j+l*BLDim];
                        phi *= A[i+l*ALDim];
                        gamma += phi;

                        phi = A[j+l*ALDim];
                        phi *= B[i+l*BLDim];
                        gamma += phi;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
    else
    {
        // C := alpha A^T B + alpha B^T A + beta C
        if( lower )
        {
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=j; i<n; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        phi = A[l+i*ALDim];
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        phi = B[l+i*BLDim];
                        phi *= A[l+j*ALDim];
                        gamma += phi;
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
                        phi = A[l+i*ALDim];
                        phi *= B[l+j*BLDim];
                        gamma += phi;

                        phi = B[l+i*BLDim];
                        phi *= A[l+j*ALDim];
                        gamma += phi;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim, 
  const Complex<DoubleDouble>* B, BlasInt BLDim,
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim, 
  const Complex<QuadDouble>* B, BlasInt BLDim,
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
template void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim, 
  const Complex<BigFloat>* B, BlasInt BLDim,
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(csyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zsyr2k)
    ( &uplo, &trans, &n, &k, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

} // namespace blas
} // namespace El
