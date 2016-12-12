/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
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
  const float* A, size_t ALDim,
        float* B, size_t BLDim );
void MKL_Domatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  double alpha,
  const double* A, size_t ALDim,
        double* B, size_t BLDim );
void MKL_Comatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex8 alpha,
  const MKL_Complex8* A, size_t ALDim,
        MKL_Complex8* B, size_t BLDim );
void MKL_Zomatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex16 alpha,
  const MKL_Complex16* A, size_t ALDim,
        MKL_Complex16* B, size_t BLDim );

void MKL_Somatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  float alpha,
  const float* A, size_t ALDim, size_t stridea,
        float* B, size_t BLDim, size_t strideb );
void MKL_Domatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  double alpha,
  const double* A, size_t ALDim, size_t stridea,
        double* B, size_t BLDim, size_t strideb );
void MKL_Comatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  MKL_Complex8 alpha,
  const MKL_Complex8* A, size_t ALDim, size_t stridea,
        MKL_Complex8* B, size_t BLDim, size_t strideb );
void MKL_Zomatcopy2
( char ordering, char trans, size_t rows, size_t cols,
  MKL_Complex16 alpha,
  const MKL_Complex16* A, size_t ALDim, size_t stridea,
        MKL_Complex16* B, size_t BLDim, size_t strideb );

void MKL_Simatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  float alpha, float* A, size_t ALDim, size_t BLDim );
void MKL_Dimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  double alpha, double* A, size_t ALDim, size_t BLDim );
void MKL_Cimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex8 alpha, MKL_Complex8* A, size_t ALDim, size_t BLDim );
void MKL_Zimatcopy
( char ordering, char trans, size_t rows, size_t cols, 
  MKL_Complex16 alpha, MKL_Complex16* A, size_t ALDim, size_t BLDim );

#ifdef EL_HAVE_MKL_GEMMT
// This routine corresponds to Elemental's Trrk (Triangular rank-k);
// Intel chose this name because it updates a single triangle using a GEMM
void EL_BLAS(sgemmt)
( const char* uplo,
  const char* transA,
  const char* transB,
  const BlasInt* n,
  const BlasInt* k,
  const float* alpha,
  const float* A, const BlasInt* ALDim,
  const float* B, const BlasInt* BLDim,
  const float* beta,
        float* C, const BlasInt* CLDim );
void EL_BLAS(dgemmt)
( const char* uplo,
  const char* transA,
  const char* transB,
  const BlasInt* n,
  const BlasInt* k,
  const double* alpha,
  const double* A, const BlasInt* ALDim,
  const double* B, const BlasInt* BLDim,
  const double* beta,
        double* C, const BlasInt* CLDim );
void EL_BLAS(cgemmt)
( const char* uplo,
  const char* transA,
  const char* transB,
  const BlasInt* n,
  const BlasInt* k,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
  const scomplex* B, const BlasInt* BLDim,
  const scomplex* beta,
        scomplex* C, const BlasInt* CLDim );
void EL_BLAS(zgemmt)
( const char* uplo,
  const char* transA,
  const char* transB,
  const BlasInt* n,
  const BlasInt* k,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
  const dcomplex* B, const BlasInt* BLDim,
  const dcomplex* beta,
        dcomplex* C, const BlasInt* CLDim );
#endif

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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
  float alpha, const float* A, BlasInt ALDim,
                     float* B, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Somatcopy( ordering, trans, m, n, alpha, A, ALDim, B, BLDim );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, const double* A, BlasInt ALDim,
                      double* B, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Domatcopy( ordering, trans, m, n, alpha, A, ALDim, B, BLDim );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, const scomplex* A, BlasInt ALDim,
                        scomplex* B, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex8 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Comatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex8*>(A), ALDim,
      reinterpret_cast<      MKL_Complex8*>(B), BLDim );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, const dcomplex* A, BlasInt ALDim,
                        dcomplex* B, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex16 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Zomatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex16*>(A), ALDim,
      reinterpret_cast<      MKL_Complex16*>(B), BLDim );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha,
  const float* A, BlasInt ALDim, BlasInt stridea,
        float* B, BlasInt BLDim, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Somatcopy2
    ( ordering, trans, m, n, alpha, A, ALDim, stridea, B, BLDim, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha,
  const double* A, BlasInt ALDim, BlasInt stridea,
        double* B, BlasInt BLDim, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Domatcopy2
    ( ordering, trans, m, n, alpha, A, ALDim, stridea, B, BLDim, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha,
  const scomplex* A, BlasInt ALDim, BlasInt stridea,
        scomplex* B, BlasInt BLDim, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex8 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Comatcopy2
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex8*>(A), ALDim, stridea,
      reinterpret_cast<      MKL_Complex8*>(B), BLDim, strideb );
}

void omatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha,
  const dcomplex* A, BlasInt ALDim, BlasInt stridea,
        dcomplex* B, BlasInt BLDim, BlasInt strideb )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex16 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Zomatcopy2
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<const MKL_Complex16*>(A), ALDim, stridea,
      reinterpret_cast<      MKL_Complex16*>(B), BLDim, strideb );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  float alpha, float* A, BlasInt ALDim, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Simatcopy( ordering, trans, m, n, alpha, A, ALDim, BLDim );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  double alpha, double* A, BlasInt ALDim, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Dimatcopy( ordering, trans, m, n, alpha, A, ALDim, BLDim );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  scomplex alpha, scomplex* A, BlasInt ALDim, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex8 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Cimatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<MKL_Complex8*>(A), ALDim, BLDim );
}

void imatcopy
( Orientation orientation, BlasInt m, BlasInt n,
  dcomplex alpha, dcomplex* A, BlasInt ALDim, BlasInt BLDim )
{
    char ordering = 'C';
    char trans = OrientationToChar( orientation );
    MKL_Complex16 alphaMKL;
    alphaMKL.real = alpha.real();
    alphaMKL.imag = alpha.imag();
    MKL_Zimatcopy
    ( ordering, trans, m, n, alphaMKL,
      reinterpret_cast<MKL_Complex16*>(A), ALDim, BLDim );
}

#ifdef EL_HAVE_MKL_GEMMT
void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        float alpha,
  const float* A, BlasInt ALDim,
  const float* B, BlasInt BLDim,
        float beta,
        float* C, BlasInt CLDim )
{
    EL_BLAS(sgemmt)
    ( &uplo, &transA, &transB,
      &n, &k,
      &alpha, A, &ALDim, B, &BLDim,
      &beta,  C, &CLDim );
}
void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        double alpha,
  const double* A, BlasInt ALDim,
  const double* B, BlasInt BLDim,
        double beta,
        double* C, BlasInt CLDim )
{
    EL_BLAS(dgemmt)
    ( &uplo, &transA, &transB,
      &n, &k,
      &alpha, A, &ALDim, B, &BLDim,
      &beta,  C, &CLDim );
}
void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        scomplex alpha,
  const scomplex* A, BlasInt ALDim,
  const scomplex* B, BlasInt BLDim,
        scomplex beta,
        scomplex* C, BlasInt CLDim )
{
    EL_BLAS(cgemmt)
    ( &uplo, &transA, &transB,
      &n, &k,
      &alpha, A, &ALDim, B, &BLDim,
      &beta,  C, &CLDim );
}
void Trrk
( char uplo, char transA, char transB,
  BlasInt n, BlasInt k,
        dcomplex alpha,
  const dcomplex* A, BlasInt ALDim,
  const dcomplex* B, BlasInt BLDim,
        dcomplex beta,
        dcomplex* C, BlasInt CLDim )
{
    EL_BLAS(zgemmt)
    ( &uplo, &transA, &transB,
      &n, &k,
      &alpha, A, &ALDim, B, &BLDim,
      &beta,  C, &CLDim );
}
#endif

} // namespace mkl
} // namespace El

#endif // ifdef EL_HAVE_MKL
