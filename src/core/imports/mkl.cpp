/*
   Copyright (c) 2009-2016, Jack Poulson
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

typedef struct _MKL_Complex8 {
  float real;
  float imag;
} MKL_Complex8;

typedef struct _MKL_Complex16 {
  double real;
  double imag;
} MKL_Complex16;

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
  const MKL_Complex8* alpha, const char* matDescrA,
  const MKL_Complex8* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const MKL_Complex8* x, const MKL_Complex8* beta, MKL_Complex8* y );
void mkl_zcsrmv
( const char* transA, const BlasInt* m, const BlasInt* k, 
  const MKL_Complex16* alpha, const char* matDescrA,
  const MKL_Complex16* val, const BlasInt* indx, 
  const BlasInt* pntrb, const BlasInt* pntre,
  const MKL_Complex16* x, const MKL_Complex16* beta, MKL_Complex16* y );

void MKL_Somatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  float alpha,
  const float* A, size_t lda,
        float* B, size_t ldb );
void MKL_Domatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  double alpha,
  const double* A, size_t lda,
        double* B, size_t ldb );
void MKL_Comatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex8 alpha,
  const MKL_Complex8* A, size_t lda,
        MKL_Complex8* B, size_t ldb );
void MKL_Zomatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex16 alpha,
  const MKL_Complex16* A, size_t lda,
        MKL_Complex16* B, size_t ldb );

void MKL_Somatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  float alpha,
  const float* A, size_t lda, size_t stridea,
        float* B, size_t ldb, size_t strideb );
void MKL_Domatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  double alpha,
  const double* A, size_t lda, size_t stridea,
        double* B, size_t ldb, size_t strideb );
void MKL_Comatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  MKL_Complex8 alpha,
  const MKL_Complex8* A, size_t lda, size_t stridea,
        MKL_Complex8* B, size_t ldb, size_t strideb );
void MKL_Zomatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  MKL_Complex16 alpha,
  const MKL_Complex16* A, size_t lda, size_t stridea,
        MKL_Complex16* B, size_t ldb, size_t strideb );

void MKL_Simatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  float alpha, float* A, size_t lda, size_t ldb );
void MKL_Dimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  double alpha, double* A, size_t lda, size_t ldb );
void MKL_Cimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex8 alpha, MKL_Complex8* A, size_t lda, size_t ldb );
void MKL_Zimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex16 alpha, MKL_Complex16* A, size_t lda, size_t ldb );

} // extern "C"

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
  float alpha, const float* A, BlasInt lda,
                     float* B, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Somatcopy( ordering, trans, m, n, alpha, A, lda, B, ldb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt lda,
                      double* B, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Domatcopy( ordering, trans, m, n, alpha, A, lda, B, ldb );
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
    MKL_Comatcopy
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
    MKL_Zomatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex16*>(A), lda,
      reinterpret_cast<      MKL_Complex16*>(B), ldb );
}

template<typename T>
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha,
  const T* A, BlasInt lda,
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
  DoubleDouble alpha,
  const DoubleDouble* A, BlasInt lda,
        DoubleDouble* B, BlasInt ldb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  QuadDouble alpha,
  const QuadDouble* A, BlasInt lda,
        QuadDouble* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_QUAD
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Quad alpha,
  const Quad* A, BlasInt lda,
        Quad* B, BlasInt ldb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Complex<Quad> alpha,
  const Complex<Quad>* A, BlasInt lda,
        Complex<Quad>* B, BlasInt ldb );
#endif
#ifdef EL_HAVE_MPC
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigInt alpha,
  const BigInt* A, BlasInt lda,
        BigInt* B, BlasInt ldb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigFloat alpha,
  const BigFloat* A, BlasInt lda,
        BigFloat* B, BlasInt ldb );
#endif

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha,
  const float* A, BlasInt lda, BlasInt stridea,
        float* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Somatcopy2
    ( ordering, trans, m, n, alpha, A, lda, stridea, B, ldb, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha,
  const double* A, BlasInt lda, BlasInt stridea,
        double* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Domatcopy2
    ( ordering, trans, m, n, alpha, A, lda, stridea, B, ldb, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha,
  const scomplex* A, BlasInt lda, BlasInt stridea,
        scomplex* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex8 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Comatcopy2
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex8*>(A), lda, stridea,
      reinterpret_cast<      MKL_Complex8*>(B), ldb, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha,
  const dcomplex* A, BlasInt lda, BlasInt stridea,
        dcomplex* B, BlasInt ldb, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex16 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Zomatcopy2
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex16*>(A), lda, stridea,
      reinterpret_cast<      MKL_Complex16*>(B), ldb, strideb );
}

template<typename T>
void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha,
  const T* A, BlasInt lda, BlasInt stridea,
        T* B, BlasInt ldb, BlasInt strideb )
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

template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Int alpha,
  const Int* A, BlasInt lda, BlasInt stridea,
        Int* B, BlasInt ldb, BlasInt strideb );
#ifdef EL_HAVE_QD
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  DoubleDouble alpha,
  const DoubleDouble* A, BlasInt lda, BlasInt stridea,
        DoubleDouble* B, BlasInt ldb, BlasInt strideb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  QuadDouble alpha,
  const QuadDouble* A, BlasInt lda, BlasInt stridea,
        QuadDouble* B, BlasInt ldb, BlasInt strideb );
#endif
#ifdef EL_HAVE_QUAD
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Quad alpha,
  const Quad* A, BlasInt lda, BlasInt stridea,
        Quad* B, BlasInt ldb, BlasInt strideb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  Complex<Quad> alpha,
  const Complex<Quad>* A, BlasInt lda, BlasInt stridea,
        Complex<Quad>* B, BlasInt ldb, BlasInt strideb );
#endif
#ifdef EL_HAVE_MPC
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigInt alpha,
  const BigInt* A, BlasInt lda, BlasInt stridea,
        BigInt* B, BlasInt ldb, BlasInt strideb );
template void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  BigFloat alpha,
  const BigFloat* A, BlasInt lda, BlasInt stridea,
        BigFloat* B, BlasInt ldb, BlasInt strideb );
#endif

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, float* A, BlasInt lda, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Simatcopy( ordering, trans, m, n, alpha, A, lda, ldb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, double* A, BlasInt lda, BlasInt ldb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Dimatcopy( ordering, trans, m, n, alpha, A, lda, ldb );
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
    MKL_Cimatcopy
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
    MKL_Zimatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<MKL_Complex16*>(A), lda, ldb );
}

template<typename T>
void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  T alpha, T* A, BlasInt lda, BlasInt ldb )
{
    LogicError("This routine not yet written");
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

} // namespace mkl
} // namespace El

#endif // ifdef EL_HAVE_MKL
