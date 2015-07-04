/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const Matrix<Ring>& A, Base<Ring> beta, Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    Syrk( uplo, orientation, Ring(alpha), A, Ring(beta), C, true );
}

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const Matrix<Ring>& A, Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, Ring(alpha), A, Ring(0), C, true );
}

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const AbstractDistMatrix<Ring>& A, 
  Base<Ring> beta,        AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    Syrk( uplo, orientation, Ring(alpha), A, Ring(beta), C, true );
}

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const AbstractDistMatrix<Ring>& A, 
                          AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, Ring(alpha), A, Ring(0), C, true );
}

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const SparseMatrix<Ring>& A,
  Base<Ring> beta,        SparseMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    Syrk( uplo, orientation, Ring(alpha), A, Ring(beta), C, true );
}

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const SparseMatrix<Ring>& A,
                          SparseMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    Syrk( uplo, orientation, Ring(alpha), A, C, true );
}

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const DistSparseMatrix<Ring>& A,
  Base<Ring> beta,        DistSparseMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    Syrk( uplo, orientation, Ring(alpha), A, Ring(beta), C, true );
}

template<typename Ring>
void Herk
( UpperOrLower uplo, Orientation orientation,
  Base<Ring> alpha, const DistSparseMatrix<Ring>& A,
                          DistSparseMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Herk"))
    Syrk( uplo, orientation, Ring(alpha), A, C, true );
}

#define PROTO(Ring) \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const Matrix<Ring>& A, \
    Base<Ring> beta,        Matrix<Ring>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const Matrix<Ring>& A, Matrix<Ring>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const AbstractDistMatrix<Ring>& A, \
                            AbstractDistMatrix<Ring>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const AbstractDistMatrix<Ring>& A, \
    Base<Ring> beta,        AbstractDistMatrix<Ring>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const SparseMatrix<Ring>& A, \
    Base<Ring> beta,        SparseMatrix<Ring>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const SparseMatrix<Ring>& A, \
                            SparseMatrix<Ring>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const DistSparseMatrix<Ring>& A, \
    Base<Ring> beta,        DistSparseMatrix<Ring>& C ); \
  template void Herk \
  ( UpperOrLower uplo, Orientation orientation, \
    Base<Ring> alpha, const DistSparseMatrix<Ring>& A, \
                            DistSparseMatrix<Ring>& C );

// blas::Herk not yet supported for Int
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
