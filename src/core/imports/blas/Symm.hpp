/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(chemm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zhemm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

void EL_BLAS(ssymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* B, const BlasInt* BLDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dsymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* B, const BlasInt* BLDim,
  const double* beta,        double* C, const BlasInt* CLDim );
void EL_BLAS(csymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zsymm)
( const char* side, const char* uplo,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Hemm
( char side, char uplo,
  BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

    // Scale C
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] = 0;

    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] *= beta;
    }

    // Naive implementation
    T gamma, delta;
    if( std::toupper(side) == 'L' && std::toupper(uplo) == 'L' )
    {
        // C := alpha tril(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=i; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha tril(A,-1)^H B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i+1; l<m; ++l )
                {
                    Conj( A[l+i*ALDim], delta );
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'L' && std::toupper(uplo) == 'U' )
    {
        // C := alpha triu(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i; l<m; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha triu(A,-)^H B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<i; ++l )
                {
                    Conj( A[l+i*ALDim], delta );
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'L' )
    {
        // C := alpha B tril(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j; l<n; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B tril(A,-1)^H + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<j; ++l )
                {
                    Conj( A[j+l*ALDim], delta );
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'U' )
    {
        // C := alpha B triu(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=j; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B triu(A,1)^H + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j+1; l<n; ++l )
                {
                    Conj( A[j+l*ALDim], delta );
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    EL_DEBUG_ONLY(
      else
          LogicError("Unsuported Hemm option");
    )
}
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Int& alpha,
  const Int* A, BlasInt ALDim,
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
  const Complex<DoubleDouble>* B, BlasInt BLDim,
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim,
  const Complex<QuadDouble>* B, BlasInt BLDim,
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim,
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
template void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim,
  const Complex<BigFloat>* B, BlasInt BLDim,
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(chemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zhemm)
    ( &side, &uplo, &m, &n,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

template<typename T>
void Symm
( char side, char uplo,
  BlasInt m, BlasInt n, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation

    // Scale C
    if( beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] = 0;

    }
    else if( beta != T(1) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] *= beta;
    }

    // Naive implementation
    T gamma, delta;
    if( std::toupper(side) == 'L' && std::toupper(uplo) == 'L' )
    {
        // C := alpha tril(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=i; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha tril(A,-1)^T B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i+1; l<m; ++l )
                {
                    delta = A[l+i*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'L' && std::toupper(uplo) == 'U' )
    {
        // C := alpha triu(A) B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=i; l<m; ++l )
                {
                    delta = A[i+l*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha triu(A,-)^T B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<i; ++l )
                {
                    delta = A[l+i*ALDim];
                    delta *= B[l+j*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'L' )
    {
        // C := alpha B tril(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j; l<n; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B tril(A,-1)^T + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<j; ++l )
                {
                    delta = A[j+l*ALDim];
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    else if( std::toupper(side) == 'R' && std::toupper(uplo) == 'U' )
    {
        // C := alpha B triu(A) + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=0; l<=j; ++l )
                {
                    delta = B[i+l*BLDim];
                    delta *= A[l+j*ALDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }

        // C := alpha B triu(A,1)^T + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt i=0; i<m; ++i )
            {
                gamma = 0;
                for( BlasInt l=j+1; l<n; ++l )
                {
                    delta = A[j+l*ALDim];
                    delta *= B[i+l*BLDim];
                    gamma += delta;
                }
                gamma *= alpha;
                C[i+j*CLDim] += gamma;
            }
        }
    }
    EL_DEBUG_ONLY(
      else
          LogicError("Unsuported Symm option");
    )
}
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Int& alpha,
  const Int* A, BlasInt ALDim, 
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim, 
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim, 
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim, 
  const Complex<DoubleDouble>* B, BlasInt BLDim,
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim, 
  const Complex<QuadDouble>* B, BlasInt BLDim,
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim, 
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim, 
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim, 
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
template void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim, 
  const Complex<BigFloat>* B, BlasInt BLDim,
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(ssymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(csymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Symm
( char side, char uplo, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zsymm)
    ( &side, &uplo, &m, &n, &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

} // namespace blas
} // namespace El
