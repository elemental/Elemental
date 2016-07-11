/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_BLAS_HPP
#define EL_IMPORTS_BLAS_HPP

namespace El {

// NOTE: The EL_LAPACK macro is defined here since many of the BLAS overloads
//       (e.g., zsyr) are technically LAPACK routines
#if defined(EL_BUILT_BLIS_LAPACK) || defined(EL_BUILT_OPENBLAS)

# ifdef EL_HAVE_BLAS_SUFFIX
#  define EL_BLAS(name) EL_CONCAT(name,EL_BLAS_SUFFIX)
# else
#  define EL_BLAS(name) FC_GLOBAL(name,name)
# endif
# ifdef EL_HAVE_LAPACK_SUFFIX
#  define EL_LAPACK(name) EL_CONCAT(name,EL_LAPACK_SUFFIX)
# else
#  define EL_LAPACK(name) FC_GLOBAL(name,name)
# endif

#else

# if defined(EL_HAVE_BLAS_SUFFIX)
#  define EL_BLAS(name) EL_CONCAT(name,EL_BLAS_SUFFIX)
# else
#  define EL_BLAS(name) name
# endif

# if defined(EL_HAVE_LAPACK_SUFFIX)
#  define EL_LAPACK(name) EL_CONCAT(name,EL_LAPACK_SUFFIX)
# else
#  define EL_LAPACK(name) name
# endif

#endif

#ifdef EL_USE_64BIT_BLAS_INTS
typedef long long int BlasInt;
#else
typedef int BlasInt;
#endif

namespace blas {

// NOTE: templated routines are custom and not wrappers

// Level 1 BLAS 
// ============
template<typename T>
void Axpy
( BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx,
        T* y, BlasInt incy );

void Axpy
( BlasInt n,
  const float& alpha, 
  const float* x, BlasInt incx,
        float* y, BlasInt incy );
void Axpy
( BlasInt n,
  const double& alpha, 
  const double* x, BlasInt incx,
        double* y, BlasInt incy );
void Axpy
( BlasInt n,
  const scomplex& alpha, 
  const scomplex* x, BlasInt incx,
        scomplex* y, BlasInt incy );
void Axpy
( BlasInt n,
  const dcomplex& alpha, 
  const dcomplex* x, BlasInt incx,
        dcomplex* y, BlasInt incy );

template<typename T>
void Copy
( BlasInt n,
  const T* x, BlasInt incx,
        T* y, BlasInt incy );

void Copy
( BlasInt n,
  const float* x, BlasInt incx,
        float* y, BlasInt incy );
void Copy
( BlasInt n,
  const double* x, BlasInt incx,
        double* y, BlasInt incy );
void Copy
( BlasInt n,
  const scomplex* x, BlasInt incx,
        scomplex* y, BlasInt incy );
void Copy
( BlasInt n,
  const dcomplex* x, BlasInt incx,
        dcomplex* y, BlasInt incy );

template<typename T>
T Dot
( BlasInt n,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy );
double Dot
( BlasInt n,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy );

template<typename T>
T Dotc
( BlasInt n,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy );
double Dotc
( BlasInt n,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy );

template<typename T>
T Dotu
( BlasInt n,
  const T* x, BlasInt incx,
  const T* y, BlasInt incy );
double Dotu
( BlasInt n,
  const double* x, BlasInt incx,
  const double* y, BlasInt incy );

template<typename F>
Base<F> Nrm2( BlasInt n, const F* x, BlasInt incx );
double Nrm2( BlasInt n, const double  * x, BlasInt incx );
double Nrm2( BlasInt n, const dcomplex* x, BlasInt incx );

template<typename F>
BlasInt MaxInd( BlasInt n, const F* x, BlasInt incx );
BlasInt MaxInd( BlasInt n, const float* x, BlasInt incx );
BlasInt MaxInd( BlasInt n, const double* x, BlasInt incx );
BlasInt MaxInd( BlasInt n, const scomplex* x, BlasInt incx );
BlasInt MaxInd( BlasInt n, const dcomplex* x, BlasInt incx );

template<typename F>
void Rot
( BlasInt n,
  F* x, BlasInt incx,
  F* y, BlasInt incy,
  const Base<F>& c,
  const F& s );
// Useful for QR algorithms with real condensed matrices
template<typename Real>
void Rot
( BlasInt n,
  Complex<Real>* x, BlasInt incx,
  Complex<Real>* y, BlasInt incy,
  const Real& c,
  const Real& s );

void Rot
( BlasInt n,
  float* x, BlasInt incx, 
  float* y, BlasInt incy,
  const float& c,
  const float& s );
void Rot
( BlasInt n,
  double* x, BlasInt incx, 
  double* y, BlasInt incy,
  const double& c,
  const double& s );
void Rot
( BlasInt n,
  scomplex* x, BlasInt incx, 
  scomplex* y, BlasInt incy,
  const float& c,
  const scomplex& s );
void Rot
( BlasInt n,
  dcomplex* x, BlasInt incx, 
  dcomplex* y, BlasInt incy,
  const double& c,
  const dcomplex& s );

template<typename T> 
void Scal( BlasInt n, const T& alpha, T*  x, BlasInt incx );
void Scal( BlasInt n, const float&    alpha, float   * x, BlasInt incx );
void Scal( BlasInt n, const double&   alpha, double  * x, BlasInt incx );
void Scal( BlasInt n, const scomplex& alpha, scomplex* x, BlasInt incx );
void Scal( BlasInt n, const dcomplex& alpha, dcomplex* x, BlasInt incx );

template<typename T> 
void Scal( BlasInt n, const T& alpha, Complex<T>* x, BlasInt incx );

// NOTE: Nrm1 is not the official name but is consistent with Nrm2
template<typename F>
Base<F> Nrm1( BlasInt n, const F* x, BlasInt incx );
double Nrm1( BlasInt n, const double  * x, BlasInt incx );
double Nrm1( BlasInt n, const Complex<double>* x, BlasInt incx );

template<typename F>
Base<F> NrmInf( BlasInt n, const F* x, BlasInt incx );

void Swap( BlasInt n, float   * x, BlasInt incx, float   * y, BlasInt incy );
void Swap( BlasInt n, double  * x, BlasInt incx, double  * y, BlasInt incy );
void Swap( BlasInt n, scomplex* x, BlasInt incx, scomplex* y, BlasInt incy );
void Swap( BlasInt n, dcomplex* x, BlasInt incx, dcomplex* y, BlasInt incy );
template<typename T> 
void Swap( BlasInt n, T* x, BlasInt incx, T* y, BlasInt incy );
            
// Level 2 BLAS
// ============
template<typename T>
void Gemv
( char trans, BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy );

void Gemv
( char trans, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy );
void Gemv
( char trans, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy );
void Gemv
( char trans, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy );
void Gemv
( char trans, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy );

template<typename T>
void Ger
( BlasInt m, BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim );

void Ger
( BlasInt m, BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim );
void Ger
( BlasInt m, BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim );
void Ger
( BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim );
void Ger
( BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim );

template<typename T>
void Geru
( BlasInt m, BlasInt n,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim );

void Geru
( BlasInt m, BlasInt n,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim );
void Geru
( BlasInt m, BlasInt n,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim );
void Geru
( BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim );
void Geru
( BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim );

template<typename T>
void Hemv
( char uplo, BlasInt m,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy );

void Hemv
( char uplo, BlasInt m,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy );
void Hemv
( char uplo, BlasInt m,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy );
void Hemv
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy );
void Hemv
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy );

template<typename T>
void Her
( char uplo, BlasInt m,
  const Base<T>& alpha,
  const T* x, BlasInt incx,
        T* A, BlasInt ALDim );

void Her
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
        float* A, BlasInt ALDim );
void Her
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
        double* A, BlasInt ALDim );
void Her
( char uplo, BlasInt m,
  const float& alpha,
  const scomplex* x, BlasInt incx,
        scomplex* A, BlasInt ALDim );
void Her
( char uplo, BlasInt m,
  const double& alpha,
  const dcomplex* x, BlasInt incx,
        dcomplex* A, BlasInt ALDim );

template<typename T>
void Her2
( char uplo, BlasInt m,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim );

void Her2
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim );
void Her2
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim );
void Her2
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim );
void Her2
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim );

template<typename T>
void Symv
( char uplo, BlasInt m,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* x, BlasInt incx,
  const T& beta,
        T* y, BlasInt incy );

void Symv
( char uplo, BlasInt m,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* x, BlasInt incx,
  const float& beta,
        float* y, BlasInt incy );
void Symv
( char uplo, BlasInt m, 
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* x, BlasInt incx,
  const double& beta,
        double* y, BlasInt incy );
void Symv
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* x, BlasInt incx,
  const scomplex& beta,
        scomplex* y, BlasInt incy );
void Symv
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* x, BlasInt incx,
  const dcomplex& beta,
        dcomplex* y, BlasInt incy );

template<typename T>
void Syr
( char uplo, BlasInt m,
  const T& alpha,
  const T* x, BlasInt incx,
        T* A, BlasInt ALDim );

void Syr
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx,
        float* A, BlasInt ALDim );
void Syr
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx,
        double* A, BlasInt ALDim );
void Syr
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx,
        scomplex* A, BlasInt ALDim ); 
void Syr
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx,
        dcomplex* A, BlasInt ALDim );

template<typename T>
void Syr2
( char uplo, BlasInt m,
  const T& alpha,
  const T* x, BlasInt incx, 
  const T* y, BlasInt incy,
        T* A, BlasInt ALDim );

void Syr2
( char uplo, BlasInt m,
  const float& alpha,
  const float* x, BlasInt incx, 
  const float* y, BlasInt incy,
        float* A, BlasInt ALDim );
void Syr2
( char uplo, BlasInt m,
  const double& alpha,
  const double* x, BlasInt incx, 
  const double* y, BlasInt incy,
        double* A, BlasInt ALDim );
void Syr2
( char uplo, BlasInt m,
  const scomplex& alpha,
  const scomplex* x, BlasInt incx, 
  const scomplex* y, BlasInt incy,
        scomplex* A, BlasInt ALDim );
void Syr2
( char uplo, BlasInt m,
  const dcomplex& alpha,
  const dcomplex* x, BlasInt incx, 
  const dcomplex* y, BlasInt incy,
        dcomplex* A, BlasInt ALDim );

template<typename T>
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const T* A, BlasInt ALDim,
        T* x, BlasInt incx );

void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt ALDim,
        float* x, BlasInt incx );
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt ALDim,
        double* x, BlasInt incx );
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt ALDim,
        scomplex* x, BlasInt incx );
void Trmv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* x, BlasInt incx );

template<typename F>
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const F* A, BlasInt ALDim,
        F* x, BlasInt incx );

void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const float* A, BlasInt ALDim,
        float* x, BlasInt incx );
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const double* A, BlasInt ALDim,
        double* x, BlasInt incx );
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const scomplex* A, BlasInt ALDim,
        scomplex* x, BlasInt incx );
void Trsv
( char uplo, char trans, char diag, BlasInt m,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* x, BlasInt incx );

// Level 3 BLAS
// ============
template<typename T>
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim );

void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim );
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim );
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim );
void Gemm
( char transA, char transB, BlasInt m, BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim );

template<typename T>
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim );

void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim );
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim );
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim );
void Hemm
( char side, char uplo, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim );

template<typename T>
void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* B, BlasInt BLDim,
  const Base<T>& beta,
        T* C, BlasInt CLDim );

void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim );
void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim );
void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const float& beta,
        scomplex* C, BlasInt CLDim );
void Her2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const double& beta,
        dcomplex* C, BlasInt CLDim );

template<typename T>
void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const Base<T>& alpha,
  const T* A, BlasInt ALDim, 
  const Base<T>& beta,
        T* C, BlasInt CLDim );

void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float& beta,
        float* C, BlasInt CLDim );
void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double& beta,
        double* C, BlasInt CLDim );
void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const scomplex* A, BlasInt ALDim,
  const float& beta,
        scomplex* C, BlasInt CLDim );
void Herk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const dcomplex* A, BlasInt ALDim,
  const double& beta,
        dcomplex* C, BlasInt CLDim );

template<typename T>
void Symm
( char side, char uplo,
  BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim );

void Symm
( char side, char uplo,
  BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim );
void Symm
( char side, char uplo,
  BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim );
void Symm
( char side, char uplo,
  BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim );
void Symm
( char side, char uplo,
  BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim );

template<typename T>
void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const T& alpha,
  const T* A, BlasInt ALDim,
  const T* B, BlasInt BLDim,
  const T& beta,
        T* C, BlasInt CLDim );

void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim, 
  const float* B, BlasInt BLDim,
  const float& beta,
        float* C, BlasInt CLDim );
void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim, 
  const double* B, BlasInt BLDim,
  const double& beta,
        double* C, BlasInt CLDim );
void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim, 
  const scomplex* B, BlasInt BLDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim );
void Syr2k
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim, 
  const dcomplex* B, BlasInt BLDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim );

template<typename T>
void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k, 
  const T& alpha,
  const T* A, BlasInt ALDim, 
  const T& beta,
        T* C, BlasInt CLDim );

void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const float& alpha,
  const float* A, BlasInt ALDim,
  const float& beta,
        float* C, BlasInt CLDim );
void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const double& alpha,
  const double* A, BlasInt ALDim,
  const double& beta,
        double* C, BlasInt CLDim );
void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex& beta,
        scomplex* C, BlasInt CLDim );
void Syrk
( char uplo, char trans,
  BlasInt n, BlasInt k,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
  const dcomplex& beta,
        dcomplex* C, BlasInt CLDim );

template<typename T>
void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim,
        T* B, BlasInt BLDim );

void Trmm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
        float* B, BlasInt BLDim );
void Trmm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
        double* B, BlasInt BLDim );
void Trmm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
        scomplex* B, BlasInt BLDim );
void Trmm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* B, BlasInt BLDim );

template<typename F>
void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const F& alpha,
  const F* A, BlasInt ALDim,
        F* B, BlasInt BLDim );

void Trsm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
        float* B, BlasInt BLDim );
void Trsm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
        double* B, BlasInt BLDim );
void Trsm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
        scomplex* B, BlasInt BLDim );
void Trsm
( char side,  char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* B, BlasInt BLDim );

} // namespace blas
} // namespace El

#endif // ifndef EL_IMPORTS_BLAS_DECL_HPP
