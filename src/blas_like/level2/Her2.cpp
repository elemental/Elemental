/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>

namespace El {

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
    EL_DEBUG_CSE
    Syr2( uplo, T(alpha), x, y, A, true );
}

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& x, const AbstractDistMatrix<T>& y,
                 AbstractDistMatrix<T>& A )
{
    EL_DEBUG_CSE
    Syr2( uplo, T(alpha), x, y, A, true );
}

#define PROTO(T) \
  template void Her2 \
  ( UpperOrLower uplo, T alpha, \
    const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A ); \
  template void Her2 \
  ( UpperOrLower uplo, T alpha, \
    const AbstractDistMatrix<T>& x, const AbstractDistMatrix<T>& y, \
          AbstractDistMatrix<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
