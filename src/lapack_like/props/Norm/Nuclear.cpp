/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
Base<Field> NuclearNorm( const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    return SchattenNorm( A, Base<Field>(1) );
}

template<typename Field>
Base<Field> HermitianNuclearNorm( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    return HermitianSchattenNorm( uplo, A, Base<Field>(1) );
}

template<typename Field>
Base<Field> SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    return SymmetricSchattenNorm( uplo, A, Base<Field>(1) );
}

template<typename Field>
Base<Field> NuclearNorm( const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    return SchattenNorm( A, Base<Field>(1) );
}

template<typename Field>
Base<Field> HermitianNuclearNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    return HermitianSchattenNorm( uplo, A, Base<Field>(1) );
}

template<typename Field>
Base<Field> SymmetricNuclearNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    return SymmetricSchattenNorm( uplo, A, Base<Field>(1) );
}

#define PROTO(Field) \
  template Base<Field> NuclearNorm( const Matrix<Field>& A ); \
  template Base<Field> NuclearNorm( const AbstractDistMatrix<Field>& A ); \
  template Base<Field> HermitianNuclearNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> HermitianNuclearNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template Base<Field> SymmetricNuclearNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> SymmetricNuclearNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
