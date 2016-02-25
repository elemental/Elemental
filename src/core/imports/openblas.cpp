/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#ifdef EL_HAVE_OPENBLAS

using El::BlasInt;
using El::scomplex;
using El::dcomplex;

extern "C" {

void EL_BLAS(somatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const float* alpha, const float* A, const BlasInt* lda,
                            float* B, const BlasInt* ldb );
void EL_BLAS(domatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const double* alpha, const double* A, const BlasInt* lda,
                             double* B, const BlasInt* ldb );
void EL_BLAS(comatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const scomplex* alpha, const scomplex* A, const BlasInt* lda,
                               scomplex* B, const BlasInt* ldb );
void EL_BLAS(zomatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const dcomplex* alpha, const dcomplex* A, const BlasInt* lda,
                               dcomplex* B, const BlasInt* ldb );

void EL_BLAS(simatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const float* alpha, float* A, const BlasInt* lda, const BlasInt* ldb );
void EL_BLAS(dimatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const double* alpha, double* A, const BlasInt* lda, const BlasInt* ldb );
void EL_BLAS(cimatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const scomplex* alpha, scomplex* A, const BlasInt* lda, const BlasInt* ldb );
void EL_BLAS(zimatcopy)
( const char* ordering, const char* trans,
  const BlasInt* rows, const BlasInt* cols, 
  const dcomplex* alpha, dcomplex* A, const BlasInt* lda, const BlasInt* ldb );

} // extern "C"

namespace El {
namespace openblas {

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda,
                     float* B, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(somatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, B, &ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda,
                      double* B, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(domatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, B, &ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda,
                        scomplex* B, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(comatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, B, &ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda,
                        dcomplex* B, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(zomatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, B, &ldb );
}

template<typename T>
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha, const T* A, BlasInt lda,
                 T* B, BlasInt ldb )
{
    if( orientation == NORMAL )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                B[i+j*ldb] = A[i+j*lda];
    }
    else if( orientation == TRANSPOSE )
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j+i*ldb] = A[i+j*lda];
    }
    else
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j+i*ldb] = Conj(A[i+j*lda]);
    }
}

template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Int alpha, const Int* A, BlasInt lda,
                   Int* B, BlasInt ldb );
#ifdef EL_HAVE_QD
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  DoubleDouble alpha, const DoubleDouble* A, BlasInt lda,
                            DoubleDouble* B, BlasInt ldb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  QuadDouble alpha, const QuadDouble* A, BlasInt lda,
                          QuadDouble* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_QUAD
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Quad alpha, const Quad* A, BlasInt lda,
                    Quad* B, BlasInt ldb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Complex<Quad> alpha, const Complex<Quad>* A, BlasInt lda,
                             Complex<Quad>* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_MPC
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigInt alpha, const BigInt* A, BlasInt lda,
                      BigInt* B, BlasInt ldb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigFloat alpha, const BigFloat* A, BlasInt lda,
                        BigFloat* B, BlasInt ldb );
#endif

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, float* A, BlasInt lda, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(simatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, &ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, double* A, BlasInt lda, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(dimatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, &ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, scomplex* A, BlasInt lda, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(cimatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, &ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, dcomplex* A, BlasInt lda, BlasInt ldb )
{
    if( m < 1 || n < 1 )
        return;
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    EL_BLAS(zimatcopy)( &ordering, &trans, &m, &n, &alpha, A, &lda, &ldb );
}

template<typename T>
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha, T* A, BlasInt lda, BlasInt ldb )
{
    LogicError("This version of imatcopy not yet supported");
}

template void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Int alpha, Int* A, BlasInt lda, BlasInt ldb );
#ifdef EL_HAVE_QD
template void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  DoubleDouble alpha, DoubleDouble* A, BlasInt lda, BlasInt ldb );
template void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  QuadDouble alpha, QuadDouble* A, BlasInt lda, BlasInt ldb );
#endif
#ifdef EL_HAVE_QUAD
template void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Quad alpha, Quad* A, BlasInt lda, BlasInt ldb );
template void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Complex<Quad> alpha, Complex<Quad>* A, BlasInt lda, BlasInt ldb );
#endif
#ifdef EL_HAVE_MPC
template void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigInt alpha, BigInt* A, BlasInt lda, BlasInt ldb );
template void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigFloat alpha, BigFloat* A, BlasInt lda, BlasInt ldb );
#endif

} // namespace openblas
} // namespace El

#endif // ifdef EL_HAVE_OPENBLAS
