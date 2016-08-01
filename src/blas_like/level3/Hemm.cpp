/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level3.hpp>

namespace El {

template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
    DEBUG_CSE
    Symm( side, uplo, alpha, A, B, beta, C, true );
}

template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C )
{
    DEBUG_CSE
    Symm( side, uplo, alpha, A, B, beta, C, true );
}

#define PROTO(T) \
  template void Hemm \
  ( LeftOrRight side, UpperOrLower uplo, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void Hemm \
  ( LeftOrRight side, UpperOrLower uplo, \
    T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, \
    T beta,        AbstractDistMatrix<T>& C );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
