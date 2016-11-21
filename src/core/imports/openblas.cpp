/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#ifdef EL_HAVE_OPENBLAS

using El::BlasInt;
using El::scomplex;
using El::dcomplex;

extern "C" {

// The following routines seem to not be defined by the shared library on Mac OS X.
// Further, there appear to be open bugs with OpenBLAS's imatcopy:
// https://github.com/xianyi/OpenBLAS/issues/899
//
// For these reasons, these functions are being temporarily disabled.
//
/*
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
*/

} // extern "C"

namespace El {
namespace openblas {

/*
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
*/

} // namespace openblas
} // namespace El

#endif // ifdef EL_HAVE_OPENBLAS
