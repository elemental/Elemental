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
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Hemm"))
    Symm( side, uplo, alpha, A, B, beta, C, true );
}

template<typename T>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Hemm"))
    Symm( side, uplo, alpha, A, B, beta, C, true );
}

#define PROTO(T) \
  template void Hemm \
  ( LeftOrRight side, UpperOrLower uplo, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void Hemm \
  ( LeftOrRight side, UpperOrLower uplo, \
    T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, \
    T beta,        DistMatrix<T>& C );

// blas::Hemm not yet supported for Int
//PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
