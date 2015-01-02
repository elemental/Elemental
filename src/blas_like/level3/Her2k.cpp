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
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha,      const Matrix<T>& A, const Matrix<T>& B, 
  Base<T> beta,       Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Her2k"))
    Syr2k( uplo, orientation, alpha, A, B, T(beta), C, true );
}

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Her2k"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syr2k( uplo, orientation, alpha, A, B, T(0), C, true );
}

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha,      const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  Base<T> beta,       AbstractDistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Her2k"))
    Syr2k( uplo, orientation, alpha, A, B, T(beta), C, true );
}

template<typename T>
void Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
                 AbstractDistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Her2k"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syr2k( uplo, orientation, alpha, A, B, T(0), C, true );
}

#define PROTO(T) \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha,      const Matrix<T>& A, const Matrix<T>& B, \
    Base<T> beta,       Matrix<T>& C ); \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
                   Matrix<T>& C ); \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, \
             const AbstractDistMatrix<T>& B, \
                   AbstractDistMatrix<T>& C ); \
  template void Her2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha,      const AbstractDistMatrix<T>& A, \
                  const AbstractDistMatrix<T>& B, \
    Base<T> beta,       AbstractDistMatrix<T>& C );

// blas::Her2k not yet supported for Int
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
