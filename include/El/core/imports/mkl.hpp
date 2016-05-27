/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_MKL_HPP
#define EL_IMPORTS_MKL_HPP

#ifdef EL_HAVE_MKL
namespace El {
namespace mkl {

void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  float alpha, const char* matDescrA,
  const float* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const float* x, float beta, float* y );
void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  double alpha, const char* matDescrA,
  const double* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const double* x, double beta, double* y );
void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  Complex<float> alpha, const char* matDescrA,
  const Complex<float>* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const Complex<float>* x, Complex<float> beta, Complex<float>* y );
void csrmv
( Orientation orientation, BlasInt m, BlasInt k,
  Complex<double> alpha, const char* matDescrA,
  const Complex<double>* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const Complex<double>* x, Complex<double> beta, Complex<double>* y );

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, float* B, BlasInt ldb );
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, double* B, BlasInt ldb );
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, scomplex* B, BlasInt ldb );
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, dcomplex* B, BlasInt ldb );

// Filler routine not provided by MKL
template<typename T>
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha,
  const T* A, BlasInt ALDim,
        T* B, BlasInt BLDim )
{
    if( orientation == NORMAL )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                B[i+j*BLDim] = A[i+j*ALDim];
    }
    else if( orientation == TRANSPOSE )
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j+i*BLDim] = A[i+j*ALDim];
    }
    else
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j+i*BLDim] = Conj(A[i+j*ALDim]);
    }
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, BlasInt stridea,
                     float* B, BlasInt ldb, BlasInt strideb );
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, BlasInt stridea,
                      double* B, BlasInt ldb, BlasInt strideb );
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, BlasInt stridea,
                        scomplex* B, BlasInt ldb, BlasInt strideb );
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, BlasInt stridea,
                        dcomplex* B, BlasInt ldb, BlasInt strideb );

// Filler routine not provided by MKL
template<typename T>
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha,
  const T* A, BlasInt ALDim, BlasInt stridea,
        T* B, BlasInt BLDim, BlasInt strideb )
{
    if( orientation == NORMAL )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                B[i*strideb+j*BLDim] = A[i*stridea+j*ALDim];
    }
    else if( orientation == TRANSPOSE )
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j*strideb+i*BLDim] = A[i*stridea+j*ALDim];
    }
    else
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j*strideb+i*BLDim] = Conj(A[i*stridea+j*ALDim]);
    }
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, float* A, BlasInt lda, BlasInt ldb );
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, double* A, BlasInt lda, BlasInt ldb );
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, scomplex* A, BlasInt lda, BlasInt ldb );
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, dcomplex* A, BlasInt lda, BlasInt ldb );

// Filler routine not provided by MKL
template<typename T>
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha, T* A, BlasInt ALDim, BlasInt BLDim )
{
    // TODO
    LogicError("This routine not yet written");
}

void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        float alpha,
  const float* A, BlasInt ALDim,
  const float* B, BlasInt BLDim,
        float beta,
        float* C, BlasInt CLDim );
void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        double alpha,
  const double* A, BlasInt ALDim,
  const double* B, BlasInt BLDim,
        double beta,
        double* C, BlasInt CLDim );
void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        scomplex alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex* B, BlasInt BLDim,
        scomplex beta,
        scomplex* C, BlasInt CLDim );
void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        dcomplex alpha,
  const dcomplex* A, BlasInt ALDim,
  const dcomplex* B, BlasInt BLDim,
        dcomplex beta,
        dcomplex* C, BlasInt CLDim );

} // namespace mkl
} // namespace El
#endif // ifdef EL_HAVE_MKL

#endif // ifndef EL_IMPORTS_MKL_HPP
