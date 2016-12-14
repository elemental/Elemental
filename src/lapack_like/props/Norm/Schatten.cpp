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
Base<Field> SchattenNorm( const Matrix<Field>& A, Base<Field> p )
{
    EL_DEBUG_CSE
    const Int minDim = Min(A.Height(),A.Width());
    return KyFanSchattenNorm( A, minDim, p );
}

template<typename Field>
Base<Field> HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<Field>& A, Base<Field> p )
{
    EL_DEBUG_CSE
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename Field>
Base<Field> SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<Field>& A, Base<Field> p )
{
    EL_DEBUG_CSE
    const Int minDim = A.Height();
    return SymmetricKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename Field>
Base<Field> SchattenNorm( const AbstractDistMatrix<Field>& A, Base<Field> p )
{
    EL_DEBUG_CSE
    const Int minDim = Min(A.Height(),A.Width());
    return KyFanSchattenNorm( A, minDim, p );
}

template<typename Field>
Base<Field> HermitianSchattenNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Base<Field> p )
{
    EL_DEBUG_CSE
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename Field>
Base<Field> SymmetricSchattenNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Base<Field> p )
{
    EL_DEBUG_CSE
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

#define PROTO(Field) \
  template Base<Field> SchattenNorm( const Matrix<Field>& A, Base<Field> p ); \
  template Base<Field> \
  SchattenNorm( const AbstractDistMatrix<Field>& A, Base<Field> p ); \
  template Base<Field> HermitianSchattenNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A, Base<Field> p ); \
  template Base<Field> HermitianSchattenNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Base<Field> p ); \
  template Base<Field> SymmetricSchattenNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A, Base<Field> p ); \
  template Base<Field> SymmetricSchattenNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Base<Field> p );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
