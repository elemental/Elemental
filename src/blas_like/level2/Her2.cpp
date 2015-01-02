/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her2"))
    Syr2( uplo, T(alpha), x, y, A, true );
}

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const AbstractDistMatrix<T>& x, const AbstractDistMatrix<T>& y,
                 AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her2"))
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

// blas::Her2 not yet supported
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
