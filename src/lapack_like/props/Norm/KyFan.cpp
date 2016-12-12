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
Base<Field> KyFanNorm( const Matrix<Field>& A, Int k )
{
    EL_DEBUG_CSE
    return KyFanSchattenNorm( A, k, Base<Field>(1) );
}

template<typename Field>
Base<Field>
HermitianKyFanNorm( UpperOrLower uplo, const Matrix<Field>& A, Int k )
{
    EL_DEBUG_CSE
    return HermitianKyFanSchattenNorm( uplo, A, k, Base<Field>(1) );
}

template<typename Field>
Base<Field>
SymmetricKyFanNorm( UpperOrLower uplo, const Matrix<Field>& A, Int k )
{
    EL_DEBUG_CSE
    return SymmetricKyFanSchattenNorm( uplo, A, k, Base<Field>(1) );
}

template<typename Field>
Base<Field> KyFanNorm( const AbstractDistMatrix<Field>& A, Int k )
{
    EL_DEBUG_CSE
    return KyFanSchattenNorm( A, k, Base<Field>(1) );
}

template<typename Field>
Base<Field> HermitianKyFanNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Int k )
{
    EL_DEBUG_CSE
    return HermitianKyFanSchattenNorm( uplo, A, k, Base<Field>(1) );
}

template<typename Field>
Base<Field> SymmetricKyFanNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Int k )
{
    EL_DEBUG_CSE
    return SymmetricKyFanSchattenNorm( uplo, A, k, Base<Field>(1) );
}

#define PROTO(Field) \
  template Base<Field> KyFanNorm( const Matrix<Field>& A, Int k ); \
  template Base<Field> KyFanNorm( const AbstractDistMatrix<Field>& A, Int k ); \
  template Base<Field> HermitianKyFanNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A, Int k ); \
  template Base<Field> HermitianKyFanNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Int k ); \
  template Base<Field> SymmetricKyFanNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A, Int k ); \
  template Base<Field> SymmetricKyFanNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Int k );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
