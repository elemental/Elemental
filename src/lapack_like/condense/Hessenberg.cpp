/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Hessenberg/LowerBlocked.hpp"
#include "./Hessenberg/UpperBlocked.hpp"
#include "./Hessenberg/ApplyQ.hpp"
#include "./Hessenberg/FormQ.hpp"

namespace El {

template<typename F>
void Hessenberg
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<F>& householderScalars )
{
    EL_DEBUG_CSE
    if( uplo == UPPER )
        hessenberg::UpperBlocked( A, householderScalars );
    else
        hessenberg::LowerBlocked( A, householderScalars );
}

template<typename F> 
void Hessenberg
( UpperOrLower uplo,
  AbstractDistMatrix<F>& A,
  AbstractDistMatrix<F>& householderScalars )
{
    EL_DEBUG_CSE
    if( uplo == UPPER )
        hessenberg::UpperBlocked( A, householderScalars );
    else
        hessenberg::LowerBlocked( A, householderScalars );
}

namespace hessenberg {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A )
{
    EL_DEBUG_CSE
    Matrix<F> householderScalars;
    Hessenberg( uplo, A, householderScalars );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F> 
void ExplicitCondensed( UpperOrLower uplo, AbstractDistMatrix<F>& A )
{
    EL_DEBUG_CSE
    DistMatrix<F,STAR,STAR> householderScalars(A.Grid());
    Hessenberg( uplo, A, householderScalars );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

} // namespace hessenberg

#define PROTO(F) \
  template void Hessenberg \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<F>& householderScalars ); \
  template void Hessenberg \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalars ); \
  template void hessenberg::ExplicitCondensed \
  ( UpperOrLower uplo, Matrix<F>& A ); \
  template void hessenberg::ExplicitCondensed \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A ); \
  template void hessenberg::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
          Matrix<F>& H ); \
  template void hessenberg::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
          AbstractDistMatrix<F>& B ); \
  template void hessenberg::FormQ \
  ( UpperOrLower uplo, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
          Matrix<F>& Q ); \
  template void hessenberg::FormQ \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
          AbstractDistMatrix<F>& Q );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
