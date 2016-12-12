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
void Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
    EL_DEBUG_CSE
    Symv( uplo, alpha, A, x, beta, y, true );
}

template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& A,
           const AbstractDistMatrix<T>& x,
  T beta,        AbstractDistMatrix<T>& y,
  const SymvCtrl<T>& ctrl )
{
    EL_DEBUG_CSE
    Symv( uplo, alpha, A, x, beta, y, true, ctrl );
}

#define PROTO(T) \
  template void Hemv \
  ( UpperOrLower uplo, T alpha, \
    const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y ); \
  template void Hemv \
  ( UpperOrLower uplo, T alpha, \
    const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& x, \
    T beta, AbstractDistMatrix<T>& y, \
    const SymvCtrl<T>& ctrl );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
