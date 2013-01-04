/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace blas {

//
// NOTE: templated routines are custom and not wrappers
//

//----------------------------------------------------------------//
// Level 1 BLAS                                                   //
//----------------------------------------------------------------//
void Axpy
( int n, float alpha, const float* x, int incx, float* y, int incy );
void Axpy
( int n, double alpha, const double* x, int incx, double* y, int incy );
void Axpy
( int n, scomplex alpha, const scomplex* x, int incx, scomplex* y, int incy );
void Axpy
( int n, dcomplex alpha, const dcomplex* x, int incx, dcomplex* y, int incy );
template<typename T>
void Axpy( int n, T alpha, const T* x, int incx, T* y, int incy );

void Copy( int n, const float* x, int incx, float* y, int incy );
void Copy( int n, const double* x, int incx, double* y, int incy );
void Copy( int n, const scomplex* x, int incx, scomplex* y, int incy );
void Copy( int n, const dcomplex* x, int incx, dcomplex* y, int incy );
template<typename T>
void Copy( int n, const T* x, int incx, T* y, int incy );

float Dot( int n, const float* x, int incx, const float* y, int incy );
double Dot( int n, const double* x, int incx, const double* y, int incy );
scomplex Dot( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dot( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dot( int n, const T* x, int incx, const T* y, int incy );

float Dotc
( int n, const float* x, int incx, const float* y, int incy );
double Dotc
( int n, const double* x, int incx, const double* y, int incy );
scomplex Dotc
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotc
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dotc( int n, const T* x, int incx, const T* y, int incy );

float Dotu
( int n, const float* x, int incx, const float* y, int incy );
double Dotu
( int n, const double* x, int incx, const double* y, int incy );
scomplex Dotu
( int n, const scomplex* x, int incx, const scomplex* y, int incy );
dcomplex Dotu
( int n, const dcomplex* x, int incx, const dcomplex* y, int incy );
template<typename T>
T Dotu( int n, const T* x, int incx, const T* y, int incy );

float Nrm2( int n, const float* x, int incx );
double Nrm2( int n, const double* x, int incx );
float Nrm2( int n, const scomplex* x, int incx );
double Nrm2( int n, const dcomplex* x, int incx );
template<typename F> F Nrm2( int n, const F* x, int incx );

void Scal( int n, float alpha, float* x, int incx );
void Scal( int n, double alpha, double* x, int incx );
void Scal( int n, scomplex alpha, scomplex* x, int incx );
void Scal( int n, dcomplex alpha, dcomplex* x, int incx );
template<typename F> void Scal( int n, F alpha, F* x, int incx );
            
//----------------------------------------------------------------//
// Level 2 BLAS                                                   //
//----------------------------------------------------------------//
void Gemv
( char trans, int m, int n,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );
void Gemv
( char trans, int m, int n,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );
void Gemv
( char trans, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );
void Gemv
( char trans, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );
template<typename T>
void Gemv
( char trans, int m, int n,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy );

void Ger
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Ger
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Ger
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Ger
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
template<typename T>
void Ger
( char trans, int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
                 T* A, int lda );

void Gerc
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Gerc
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Gerc
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Gerc
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
template<typename T>
void Gerc
( char trans, int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
                 T* A, int lda );

void Geru
( int m, int n,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Geru
( int m, int n,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Geru
( int m, int n,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Geru
( int m, int n,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
template<typename T>
void Geru
( char trans, int m, int n,
  T alpha, const T* x, int incx, const T* y, int incy,
                 T* A, int lda );

void Hemv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );
void Hemv
( char uplo, int m,
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );
void Hemv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );
void Hemv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );
template<typename T>
void Hemv
( char uplo, int m,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy );

void Her
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda );
void Her
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda );
void Her
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda );
void Her
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda );
template<typename T>
void Hemv( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda );

void Her2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Her2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Her2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Her2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
template<typename T>
void Her2
( char uplo, int m,
  T alpha, const T* x, int incx, const T* y, int incy, 
                 T* A, int lda );

void Symv
( char uplo, int m,
  float alpha, const float* A, int lda, const float* x, int incx,
  float beta,        float* y, int incy );
void Symv
( char uplo, int m, 
  double alpha, const double* A, int lda, const double* x, int incx,
  double beta,        double* y, int incy );
void Symv
( char uplo, int m,
  scomplex alpha, const scomplex* A, int lda, const scomplex* x, int incx,
  scomplex beta,        scomplex* y, int incy );
void Symv
( char uplo, int m,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* x, int incx,
  dcomplex beta,        dcomplex* y, int incy );
template<typename T>
void Symv
( char uplo, int m,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy );

void Syr
( char uplo, int m,
  float alpha, const float* x, int incx, float* A, int lda );
void Syr
( char uplo, int m,
  double alpha, const double* x, int incx, double* A, int lda );
void Syr
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, scomplex* A, int lda ); 
void Syr
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, dcomplex* A, int lda );
template<typename T>
void Syr( char uplo, int m, T alpha, const T* x, int incx, T* A, int lda );

void Syr2
( char uplo, int m,
  float alpha, const float* x, int incx, const float* y, int incy,
                     float* A, int lda );
void Syr2
( char uplo, int m,
  double alpha, const double* x, int incx, const double* y, int incy,
                      double* A, int lda );
void Syr2
( char uplo, int m,
  scomplex alpha, const scomplex* x, int incx, const scomplex* y, int incy,
                        scomplex* A, int lda );
void Syr2
( char uplo, int m,
  dcomplex alpha, const dcomplex* x, int incx, const dcomplex* y, int incy,
                        dcomplex* A, int lda );
template<typename T>
void Syr2
( char uplo, int m,
  T alpha, const T* x, int incx, const T* y, int incy,
                 T* A, int lda );

void Trmv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx );
void Trmv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx );
void Trmv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx );
void Trmv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx );
template<typename T>
void Trmv
( char uplo, char trans, char diag, int m,
  const T* A, int lda, T* x, int incx );

void Trsv
( char uplo, char trans, char diag, int m,
  const float* A, int lda, float* x, int incx );
void Trsv
( char uplo, char trans, char diag, int m,
  const double* A, int lda, double* x, int incx );
void Trsv
( char uplo, char trans, char diag, int m,
  const scomplex* A, int lda, scomplex* x, int incx );
void Trsv
( char uplo, char trans, char diag, int m,
  const dcomplex* A, int lda, dcomplex* x, int incx );
template<typename T>
void Trsv
( char uplo, char trans, char diag, int m,
  const T* A, int lda, T* x, int incx );

//----------------------------------------------------------------//
// Level 3 BLAS                                                   //
//----------------------------------------------------------------//
void Gemm
( char transA, char transB, int m, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Gemm
( char transA, char transB, int m, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Gemm
( char transA, char transB, int m, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Gemm
( char transA, char transB, int m, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Gemm
( char transA, char transB, int m, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Hemm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Hemm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Hemm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Hemm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Hemm
( char side, char uplo, int m, int n,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Her2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Her2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Her2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Her2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Her2k
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Herk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, float beta, float* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, double beta, double* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc );
void Herk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Herk
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda,
  T beta,        T* C, int ldc );

void Symm
( char side, char uplo, int m, int n,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Symm
( char side, char uplo, int m, int n,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Symm
( char side, char uplo, int m, int n,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Symm
( char side, char uplo, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Symm
( char side, char uplo, int m, int n,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Syr2k
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda, const float* B, int ldb,
  float beta,        float* C, int ldc );
void Syr2k
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda, const double* B, int ldb,
  double beta,        double* C, int ldc );
void Syr2k
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda, const scomplex* B, int ldb,
  scomplex beta,        scomplex* C, int ldc );
void Syr2k
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda, const dcomplex* B, int ldb,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Syr2k
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc );

void Syrk
( char uplo, char trans, int n, int k,
  float alpha, const float* A, int lda,
  float beta,        float* C, int ldc );
void Syrk
( char uplo, char trans, int n, int k,
  double alpha, const double* A, int lda,
  double beta,        double* C, int ldc );
void Syrk
( char uplo, char trans, int n, int k,
  scomplex alpha, const scomplex* A, int lda,
  scomplex beta,        scomplex* C, int ldc );
void Syrk
( char uplo, char trans, int n, int k,
  dcomplex alpha, const dcomplex* A, int lda,
  dcomplex beta,        dcomplex* C, int ldc );
template<typename T>
void Syrk
( char uplo, char trans, int n, int k,
  T alpha, const T* A, int lda,
  T beta,        T* C, int ldc );

void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb );
void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb );
void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb );
void Trmm
( char side,  char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb );
template<typename T>
void Trmm
( char side, char uplo, char trans, char unit, int m, int n,
  T alpha, const T* A, int lda, T* B, int ldb );

void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  float alpha, const float* A, int lda, float* B, int ldb );
void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  double alpha, const double* A, int lda, double* B, int ldb );
void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  scomplex alpha, const scomplex* A, int lda, scomplex* B, int ldb );
void Trsm
( char side,  char uplo, char trans, char unit, int m, int n,
  dcomplex alpha, const dcomplex* A, int lda, dcomplex* B, int ldb );
template<typename T>
void Trsm
( char side, char uplo, char trans, char unit, int m, int n,
  T alpha, const T* A, int lda, T* B, int ldb );

} // namespace blas
} // namespace elem

//
// Templated wrappers
//

namespace elem {
namespace blas {

//
// Level 1 BLAS
//

template<typename T>
inline void Axpy
( int n, T alpha, const T* x, int incx, T* y, int incy )
{
    for( int i=0; i<n; ++i )
        y[i*incy] += alpha*x[i*incx];
}

template<typename T>
inline void Copy
( int n, const T* x, int incx, T* y, int incy )
{
    for( int i=0; i<n; ++i )
        y[i*incy] = x[i*incx];
}

template<typename T>
inline T Dot( int n, const T* x, int incx, const T* y, int incy )
{ Dotc( n, x, incx, y, incy ); }

template<typename T>
inline T Dotc( int n, const T* x, int incx, const T* y, int incy )
{
    T alpha = 0;
    for( int i=0; i<n; ++i )
        alpha += Conj(x[i*incx])*y[i*incy];
    return alpha;
}

template<typename T>
inline T Dotu( int n, const T* x, int incx, const T* y, int incy )
{
    T alpha = 0;
    for( int i=0; i<n; ++i )
        alpha += x[i*incx]*y[i*incy];
    return alpha;
}

// TODO: templated Nrm2

template<typename T>
void Scal( int n, T alpha, T* x, int incx )
{
    for( int i=0; i<n; ++i )
        x[i*incx] *= alpha;
}

// 
// Level 2 BLAS
//

template<typename T>
void Gemv
( char trans, int m, int n,
  T alpha, const T* A, int lda, const T* x, int incx,
  T beta,        T* y, int incy )
{
    if( trans == 'N' )
    {
        if( m > 0 && n == 0 && beta == 0 )
        {
            for( int i=0; i<m; ++i )
                y[i*incy] = 0;   
            return;
        }
        Scal( m, beta, y, incy );
        for( int i=0; i<m; ++i ) 
            for( int j=0; j<n; ++j )
                y[i*incy] += alpha*A[i+j*lda]*x[j*incx];
    }
    else if( trans == 'T' ) 
    {
        if( n > 0 && m == 0 && beta == 0 )
        {
            for( int i=0; i<n; ++i )
                y[i*incy] = 0;   
            return;
        }
        Scal( n, beta, y, incy );
        for( int i=0; i<n; ++i ) 
            for( int j=0; j<m; ++j )
                y[i*incy] += alpha*A[j+i*lda]*x[j*incx];
    }
    else
    {
        if( n > 0 && m == 0 && beta == 0 )
        {
            for( int i=0; i<n; ++i )
                y[i*incy] = 0;   
            return;
        }
        Scal( n, beta, y, incy );
        for( int i=0; i<n; ++i ) 
            for( int j=0; j<m; ++j )
                y[i*incy] += alpha*Conj(A[j+i*lda])*x[j*incx];
    }
}

// TODO: templated Ger
// TODO: templated Gerc
// TODO: templated Geru
// TODO: templated Hemv
// TODO: templated Her
// TODO: templated Her2
// TODO: templated Symv
// TODO: templated Syr
// TODO: templated Syr2
// TODO: templated Trmv
// TODO: templated Trsv

//
// Level 3 BLAS
//

template<typename T>
void Gemm
( char transA, char transB, int m, int n, int k,
  T alpha, const T* A, int lda, const T* B, int ldb,
  T beta,        T* C, int ldc )
{
    if( m > 0 && n > 0 && k == 0 && beta == 0 )
    {
        for( int j=0; j<n; ++j )
            for( int i=0; i<m; ++i )
                C[i+j*ldc] = 0;
        return;
    }

    // Scale C
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            C[i+j*ldc] *= beta;

    // Naive implementation
    if( transA == 'N' && transB == 'N' )
    {
        // C := alpha A B + C
        for( int j=0; j<n; ++j )
            for( int i=0; i<m; ++i )
                for( int l=0; l<k; ++l )
                    C[i+j*ldc] += alpha*A[i+l*lda]*B[l+j*ldb];
    }
    else if( transA == 'N' )
    {
        if( transB == 'T' )
        {
            // C := alpha A B^T + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[i+l*lda]*B[j+l*ldb];
        }
        else
        {
            // C := alpha A B^H + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[i+l*lda]*Conj(B[j+l*ldb]);
        }
    }
    else if( transB == 'N' )
    {
        if( transA == 'T' )
        {
            // C := alpha A^T B + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*B[l+j*ldb];
        }
        else
        {
            // C := alpha A^H B + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*B[l+j*ldb];
        }
    }
    else
    {
        if( transA == 'T' && transB == 'T' )
        {
            // C := alpha A^T B^T + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*B[j+l*ldb];
        }
        else if( transA == 'T' )
        {
            // C := alpha A^T B^H + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*A[l+i*lda]*Conj(B[j+l*ldb]);
        }
        else if( transB == 'T' )
        {
            // C := alpha A^H B^T + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*B[j+l*ldb];
        }
        else
        {
            // C := alpha A^H B^H + C
            for( int j=0; j<n; ++j )
                for( int i=0; i<m; ++i )
                    for( int l=0; l<k; ++l )
                        C[i+j*ldc] += alpha*Conj(A[l+i*lda])*Conj(B[j+l*ldb]);
        }
    }
}

// TODO: templated Hemm
// TODO: templated Her2k
// TODO: templated Herk
// TODO: templated Symm
// TODO: templated Syr2k
// TODO: templated Syrk
// TODO: templated Trmm
// TODO: templated Trsm

} // namespace blas
} // namespace elem
