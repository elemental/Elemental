/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her2"))
    Syr2( uplo, alpha, x, y, A, true );
}

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x,
           const DistMatrix<T>& y,
                 DistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her2"))
    Syr2( uplo, alpha, x, y, A, true );
}

#define PROTO(T) \
  template void Her2 \
  ( UpperOrLower uplo, T alpha, \
    const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A ); \
  template void Her2 \
  ( UpperOrLower uplo, T alpha, \
    const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A );

// blas::Her2 not yet supported
//PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
