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
Base<Field> LogBarrier( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename Field>
Base<Field> LogBarrier( UpperOrLower uplo, Matrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet =
      SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

template<typename Field>
Base<Field> LogBarrier( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename Field>
Base<Field> LogBarrier
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet =
      SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

#define PROTO(Field) \
  template Base<Field> LogBarrier \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> LogBarrier \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template Base<Field> LogBarrier \
  ( UpperOrLower uplo, Matrix<Field>& A, bool canOverwrite ); \
  template Base<Field> LogBarrier \
  ( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool canOverwrite );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
