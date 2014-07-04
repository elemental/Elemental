/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Hemv"))
    Symv( uplo, alpha, A, x, beta, y, true );
}

template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Hemv"))
    Symv( uplo, alpha, A, x, beta, y, true );
}

#define PROTO(T) \
  template void Hemv \
  ( UpperOrLower uplo, T alpha, \
    const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y ); \
  template void Hemv \
  ( UpperOrLower uplo, T alpha, \
    const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y );

// blas::Hemv not yet supported
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
