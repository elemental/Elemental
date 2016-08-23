/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {


void EL_BLAS(sgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const float* alpha, const float* A, const BlasInt* ALDim,
                      const float* B, const BlasInt* BLDim,
  const float* beta,        float* C, const BlasInt* CLDim );
void EL_BLAS(dgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const double* alpha, const double* A, const BlasInt* ALDim,
                       const double* B, const BlasInt* BLDim,
  const double* beta,        double* C, const BlasInt* CLDim );
void EL_BLAS(cgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zgemm)
( const char* transA, const char* transB,
  const BlasInt* m, const BlasInt* n, const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k,
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    if( m > 0 && n > 0 && k == 0 && beta == T(0) )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                C[i+j*CLDim] = 0;
        return;
    }

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
    if( std::toupper(transA) == 'N' && std::toupper(transB) == 'N' )
    {
        // C := alpha A B + C
        for( BlasInt j=0; j<n; ++j )
        {
            for( BlasInt l=0; l<k; ++l )
            {
                gamma = alpha;
                gamma *= B[l+j*BLDim];
                for( BlasInt i=0; i<m; ++i )
                {
                    delta = A[i+l*ALDim];
                    delta *= gamma;
                    C[i+j*CLDim] += delta;
                }
            }
        }
    }
    else if( std::toupper(transA) == 'N' )
    {
        if( std::toupper(transB) == 'T' )
        {
            // C := alpha A B^T + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt l=0; l<k; ++l )
                {
                    gamma = alpha;
                    gamma *= B[j+l*BLDim];
                    for( BlasInt i=0; i<m; ++i )
                    {
                        delta = A[i+l*ALDim];
                        delta *= gamma;
                        C[i+j*CLDim] += delta;
                    }
                }
            }
        }
        else
        {
            // C := alpha A B^H + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt l=0; l<k; ++l )
                {
                    Conj( B[j+l*BLDim], gamma );
                    gamma *= alpha;
                    for( BlasInt i=0; i<m; ++i )
                    {
                        delta = A[i+l*ALDim];
                        delta *= gamma;
                        C[i+j*CLDim] += delta;
                    }
                }
            }
        }
    }
    else if( std::toupper(transB) == 'N' )
    {
        if( std::toupper(transA) == 'T' )
        {
            // C := alpha A^T B + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
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
        else
        {
            // C := alpha A^H B + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
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
    }
    else
    {
        if( std::toupper(transA) == 'T' && std::toupper(transB) == 'T' )
        {
            // C := alpha A^T B^T + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        delta = A[l+i*ALDim];
                        delta *= B[j+l*BLDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else if( std::toupper(transA) == 'T' )
        {
            // C := alpha A^T B^H + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( B[j+l*BLDim], delta );
                        delta *= A[l+i*ALDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else if( std::toupper(transB) == 'T' )
        {
            // C := alpha A^H B^T + C
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );
                        delta *= B[j+l*BLDim];
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
        else
        {
            // C := alpha A^H B^H + C
            T phi;
            for( BlasInt j=0; j<n; ++j )
            {
                for( BlasInt i=0; i<m; ++i )
                {
                    gamma = 0;
                    for( BlasInt l=0; l<k; ++l )
                    {
                        Conj( A[l+i*ALDim], delta );       
                        Conj( B[j+l*BLDim], phi );
                        delta *= phi;
                        gamma += delta;
                    }
                    gamma *= alpha;
                    C[i+j*CLDim] += gamma;
                }
            }
        }
    }
}
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Int& alpha,
  const Int* A, BlasInt ALDim,
  const Int* B, BlasInt BLDim,
  const Int& beta,
        Int* C, BlasInt CLDim );
#ifdef EL_HAVE_QD
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
  const DoubleDouble* B, BlasInt BLDim,
  const DoubleDouble& beta,
        DoubleDouble* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
  const QuadDouble* B, BlasInt BLDim,
  const QuadDouble& beta,
        QuadDouble* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
  const Complex<DoubleDouble>* B, BlasInt BLDim,
  const Complex<DoubleDouble>& beta,
        Complex<DoubleDouble>* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim,
  const Complex<QuadDouble>* B, BlasInt BLDim,
  const Complex<QuadDouble>& beta,
        Complex<QuadDouble>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
  const Quad* B, BlasInt BLDim,
  const Quad& beta,
        Quad* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim, 
  const Complex<Quad>* B, BlasInt BLDim,
  const Complex<Quad>& beta,
        Complex<Quad>* C, BlasInt CLDim );
#endif
#ifdef EL_HAVE_MPC
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
  const BigInt* B, BlasInt BLDim,
  const BigInt& beta,
        BigInt* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
  const BigFloat* B, BlasInt BLDim,
  const BigFloat& beta,
        BigFloat* C, BlasInt CLDim );
template void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim,
  const Complex<BigFloat>* B, BlasInt BLDim,
  const Complex<BigFloat>& beta,
        Complex<BigFloat>* C, BlasInt CLDim );
#endif

void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    const char fixedTransA = ( std::toupper(transA) == 'C' ? 'T' : transA );
    const char fixedTransB = ( std::toupper(transB) == 'C' ? 'T' : transB );
    EL_BLAS(sgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Gemm
( char transA, char transB,
  BlasInt m, BlasInt n, BlasInt k, 
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }      

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    const char fixedTransA = ( std::toupper(transA) == 'C' ? 'T' : transA );
    const char fixedTransB = ( std::toupper(transB) == 'C' ? 'T' : transB );
    EL_BLAS(dgemm)
    ( &fixedTransA, &fixedTransB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }      

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    EL_BLAS(cgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k, 
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( std::toupper(transA) == 'N' )
      {
          if( ALDim < Max(m,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",m=",m);
      }
      else
      {
          if( ALDim < Max(k,1) )
              LogicError("ALDim was too small: ALDim=",ALDim,",k=",k);
      }      

      if( std::toupper(transB) == 'N' )
      {
          if( BLDim < Max(k,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",k=",k);
      }
      else
      {
          if( BLDim < Max(n,1) )
              LogicError("BLDim was too small: BLDim=",BLDim,",n=",n);
      }

      if( CLDim < Max(m,1) )
          LogicError("CLDim was too small: CLDim=",CLDim,",m=",m);
    )
    EL_BLAS(zgemm)
    ( &transA, &transB, &m, &n, &k,
      &alpha, A, &ALDim, B, &BLDim, &beta, C, &CLDim );
}

} // namespace blas
} // namespace El
