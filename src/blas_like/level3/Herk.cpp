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
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const Matrix<T>& A, Base<T> beta, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    Syrk( uplo, orientation, T(alpha), A, T(beta), C, true );
}

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const Matrix<T>& A, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, T(alpha), A, T(0), C, true );
}

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const AbstractDistMatrix<T>& A, 
  Base<T> beta,        AbstractDistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    Syrk( uplo, orientation, T(alpha), A, T(beta), C, true );
}

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, T(alpha), A, T(0), C, true );
}

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const SparseMatrix<T>& A,
  Base<T> beta,        SparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    Syrk( uplo, orientation, T(alpha), A, T(beta), C, true );
}

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const SparseMatrix<T>& A,
                       SparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    Syrk( uplo, orientation, T(alpha), A, C, true );
}

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const DistSparseMatrix<T>& A,
  Base<T> beta,        DistSparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    Syrk( uplo, orientation, T(alpha), A, T(beta), C, true );
}

template<typename T>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const DistSparseMatrix<T>& A,
                       DistSparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Herk"))
    Syrk( uplo, orientation, T(alpha), A, C, true );
}

#define PROTO(T) \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const Matrix<T>& A, \
    Base<T> beta,        Matrix<T>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const Matrix<T>& A, Matrix<T>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const AbstractDistMatrix<T>& A, \
    Base<T> beta,        AbstractDistMatrix<T>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const SparseMatrix<T>& A, \
    Base<T> beta,        SparseMatrix<T>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const SparseMatrix<T>& A, \
                         SparseMatrix<T>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const DistSparseMatrix<T>& A, \
    Base<T> beta,        DistSparseMatrix<T>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<T> alpha, const DistSparseMatrix<T>& A, \
                         DistSparseMatrix<T>& C );

// blas::Herk not yet supported for Int
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
