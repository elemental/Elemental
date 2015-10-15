/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#ifdef EL_HAVE_MKL
#include "mkl.h"

using El::BlasInt;
using El::scomplex;
using El::dcomplex;

namespace El {
namespace mkl {

// TODO: Ensure that BlasInt is compatible with MKL_INT in the following
//       routines

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
    ( &transA, &m, &k, 
      reinterpret_cast<const MKL_Complex8*>(&alpha),
      matDescrA,
      reinterpret_cast<const MKL_Complex8*>(val),
      indx, pntrb, pntre,
      reinterpret_cast<const MKL_Complex8*>(x),
      reinterpret_cast<const MKL_Complex8*>(&beta),
      reinterpret_cast<      MKL_Complex8*>(y) );
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
    ( &transA, &m, &k,
      reinterpret_cast<const MKL_Complex16*>(&alpha),
      matDescrA,
      reinterpret_cast<const MKL_Complex16*>(val),
      indx, pntrb, pntre,
      reinterpret_cast<const MKL_Complex16*>(x),
      reinterpret_cast<const MKL_Complex16*>(&beta),
      reinterpret_cast<      MKL_Complex16*>(y) );
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
    MKL_Complex8 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    mkl_comatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex8*>(A), lda,
      reinterpret_cast<      MKL_Complex8*>(B), ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda,
                        dcomplex* B, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex16 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    mkl_zomatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex16*>(A), lda,
      reinterpret_cast<      MKL_Complex16*>(B), ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Int alpha, const Int* A, BlasInt lda, BlasInt stridea,
                   Int* B, BlasInt ldb, BlasInt strideb )
{
    if( orientation == NORMAL )
    {
        for( BlasInt j=0; j<n; ++j )
            for( BlasInt i=0; i<m; ++i )
                B[i*strideb+j*ldb] = A[i*stridea+j*lda];
    }
    else if( orientation == TRANSPOSE )
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j*strideb+i*ldb] = A[i*stridea+j*lda];
    }
    else
    {
        for( BlasInt i=0; i<m; ++i )
            for( BlasInt j=0; j<n; ++j )
                B[j*strideb+i*ldb] = Conj(A[i*stridea+j*lda]);
    }
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, const float* A, BlasInt lda, BlasInt stridea,
                     float* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_somatcopy2
    ( ordering, trans, m, n, alpha, A, lda, stridea, B, ldb, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda, BlasInt stridea,
                      double* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    mkl_domatcopy2
    ( ordering, trans, m, n, alpha, A, lda, stridea, B, ldb, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt lda, BlasInt stridea,
                        scomplex* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex8 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    mkl_comatcopy2
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex8*>(A), lda, stridea,
      reinterpret_cast<      MKL_Complex8*>(B), ldb, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt lda, BlasInt stridea,
                        dcomplex* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex16 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    mkl_zomatcopy2
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex16*>(A), lda, stridea,
      reinterpret_cast<      MKL_Complex16*>(B), ldb, strideb );
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
    MKL_Complex8 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    mkl_cimatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<MKL_Complex8*>(A), lda, ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, dcomplex* A, BlasInt lda, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex16 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    mkl_zimatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<MKL_Complex16*>(A), lda, ldb );
}

} // namespace mkl
} // namespace El

#endif // ifdef EL_HAVE_MKL
