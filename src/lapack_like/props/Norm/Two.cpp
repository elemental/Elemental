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
Base<Field> TwoNorm( const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    Matrix<Base<Field>> s;
    SVD( A, s );
    return InfinityNorm( s );
}

template<typename Field>
Base<Field> HermitianTwoNorm( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    Matrix<Base<Field>> s;
    HermitianSVD( uplo, A, s );
    return InfinityNorm( s );
}

template<typename Field>
Base<Field> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    Matrix<Field> B( A );
    Matrix<Base<Field>> s;
    MakeSymmetric( uplo, B );
    SVDCtrl<Base<Field>> ctrl;
    ctrl.overwrite = true;
    SVD( B, s, ctrl );
    return MaxNorm( s );
}

template<typename Field>
Base<Field> TwoNorm( const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    DistMatrix<Base<Field>,VR,STAR> s( A.Grid() );
    SVD( A, s );
    return InfinityNorm( s );
}

template<typename Field>
Base<Field>
HermitianTwoNorm( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    DistMatrix<Base<Field>,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, A, s );
    return InfinityNorm( s );
}

template<typename Field>
Base<Field>
SymmetricTwoNorm( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    DistMatrix<Field> B( A );
    DistMatrix<Base<Field>,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVDCtrl<Base<Field>> ctrl;
    ctrl.overwrite = true;
    SVD( B, s, ctrl );
    return MaxNorm( s );
}

#define PROTO(Field) \
  template Base<Field> TwoNorm( const Matrix<Field>& A ); \
  template Base<Field> TwoNorm( const AbstractDistMatrix<Field>& A ); \
  template Base<Field> \
  HermitianTwoNorm( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> HermitianTwoNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template Base<Field> \
  SymmetricTwoNorm( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> SymmetricTwoNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
