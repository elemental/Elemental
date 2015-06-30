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

} // namespace mkl
} // namespace El

#endif // ifdef EL_HAVE_MKL
