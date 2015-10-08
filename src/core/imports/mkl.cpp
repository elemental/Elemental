/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#ifdef EL_HAVE_MKL

using El::BlasInt;
using El::scomplex;
using El::dcomplex;

extern "C" {

void mkl_scsrmv
( const char* transA, const BlasInt* m, const BlasInt* k, 
  const float* alpha, const char* matDescrA,
  const float* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const float* x, const float* beta, float* y );
void mkl_dcsrmv
( const char* transA, const BlasInt* m, const BlasInt* k, 
  const double* alpha, const char* matDescrA,
  const double* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const double* x, const double* beta, double* y );
void mkl_ccsrmv
( const char* transA, const BlasInt* m, const BlasInt* k, 
  const scomplex* alpha, const char* matDescrA,
  const scomplex* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const scomplex* x, const scomplex* beta, scomplex* y );
void mkl_zcsrmv
( const char* transA, const BlasInt* m, const BlasInt* k, 
  const dcomplex* alpha, const char* matDescrA,
  const dcomplex* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const dcomplex* x, const dcomplex* beta, dcomplex* y );

void mkl_somatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  float alpha, const float* A, size_t lda,
                     float* B, size_t ldb );
void mkl_domatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  double alpha, const double* A, size_t lda,
                      double* B, size_t ldb );
void mkl_comatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  scomplex alpha, const scomplex* A, size_t lda,
                        scomplex* B, size_t ldb );
void mkl_zomatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  dcomplex alpha, const dcomplex* A, size_t lda,
                        dcomplex* B, size_t ldb );

void mkl_simatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  float alpha, float* A, size_t lda, size_t ldb );
void mkl_dimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  double alpha, double* A, size_t lda, size_t ldb );
void mkl_cimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  scomplex alpha, scomplex* A, size_t lda, size_t ldb );
void mkl_zimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  dcomplex alpha, dcomplex* A, size_t lda, size_t ldb );

} // extern "C"

namespace El {
namespace mkl {

// NOTE: For the usual Elemental sparse format, we can simply set 
//           ptre = ptrb+1
//       and 
//           matDescrA[0] = 'G'
//           matDescrA[3] = 'C'
void csrmv
( Orientation orientation, BlasInt m, BlasInt k, 
  float alpha, const char* matDescrA, 
  const float* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const float* x, float beta, float* y )
{
    DEBUG_ONLY(CSE cse("mkl::csrmv"))
    char transA = OrientationToChar( orientation );
    if( transA == 'C' )
        transA = 'T';
    mkl_scsrmv
    ( &transA, &m, &k, &alpha, matDescrA, val, indx, pntrb, pntre,
      x, &beta, y );
}

void csrmv
( Orientation orientation, BlasInt m, BlasInt k, 
  double alpha, const char* matDescrA, 
  const double* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const double* x, double beta, double* y )
{
    DEBUG_ONLY(CSE cse("mkl::csrmv"))
    char transA = OrientationToChar( orientation );
    if( transA == 'C' )
        transA = 'T';
    mkl_dcsrmv
    ( &transA, &m, &k, &alpha, matDescrA, val, indx, pntrb, pntre,
      x, &beta, y );
}

void csrmv
( Orientation orientation, BlasInt m, BlasInt k, 
  Complex<float> alpha, const char* matDescrA, 
  const Complex<float>* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const Complex<float>* x, Complex<float> beta, Complex<float>* y )
{
    DEBUG_ONLY(CSE cse("mkl::csrmv"))
    char transA = OrientationToChar( orientation );
    mkl_ccsrmv
    ( &transA, &m, &k, &alpha, matDescrA, val, indx, pntrb, pntre,
      x, &beta, y );
}

void csrmv
( Orientation orientation, BlasInt m, BlasInt k, 
  Complex<double> alpha, const char* matDescrA, 
  const Complex<double>* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const Complex<double>* x, Complex<double> beta, Complex<double>* y )
{
    DEBUG_ONLY(CSE cse("mkl::csrmv"))
    char transA = OrientationToChar( orientation );
    mkl_zcsrmv
    ( &transA, &m, &k, &alpha, matDescrA, val, indx, pntrb, pntre,
      x, &beta, y );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Int alpha, const Int* A, BlasInt lda,
                   Int* B, BlasInt ldb )
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

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda,
                     float* B, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_somatcopy( ordering, trans, m, n, alpha, A, lda, B, ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda,
                      double* B, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_domatcopy( ordering, trans, m, n, alpha, A, lda, B, ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda,
                        scomplex* B, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_comatcopy( ordering, trans, m, n, alpha, A, lda, B, ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda,
                        dcomplex* B, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_zomatcopy( ordering, trans, m, n, alpha, A, lda, B, ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Int alpha, Int* A, BlasInt lda, BlasInt ldb )
{
    LogicError("Integer MKL imatcopy not yet supported");
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, float* A, BlasInt lda, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_simatcopy( ordering, trans, m, n, alpha, A, lda, ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, double* A, BlasInt lda, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_dimatcopy( ordering, trans, m, n, alpha, A, lda, ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, scomplex* A, BlasInt lda, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_cimatcopy( ordering, trans, m, n, alpha, A, lda, ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, dcomplex* A, BlasInt lda, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_zimatcopy( ordering, trans, m, n, alpha, A, lda, ldb );
}

} // namespace mkl
} // namespace El

#endif // ifdef EL_HAVE_MKL
