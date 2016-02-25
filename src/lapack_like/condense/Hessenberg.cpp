/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Hessenberg/L.hpp"
#include "./Hessenberg/U.hpp"
#include "./Hessenberg/ApplyQ.hpp"

namespace El {

template<typename F>
void Hessenberg( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CSE cse("Hessenberg"))
    if( uplo == UPPER )
        hessenberg::U( A, t );
    else
        hessenberg::L( A, t );
}

template<typename F> 
void Hessenberg
( UpperOrLower uplo, ElementalMatrix<F>& A, ElementalMatrix<F>& t )
{
    DEBUG_ONLY(CSE cse("Hessenberg"))
    if( uplo == UPPER )
        hessenberg::U( A, t );
    else
        hessenberg::L( A, t );
}

namespace hessenberg {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("hessenberg::ExplicitCondensed"))
    Matrix<F> t;
    Hessenberg( uplo, A, t );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F> 
void ExplicitCondensed( UpperOrLower uplo, ElementalMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("hessenberg::ExplicitCondensed"))
    DistMatrix<F,STAR,STAR> t(A.Grid());
    Hessenberg( uplo, A, t );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

} // namespace hessenberg

#define PROTO(F) \
  template void Hessenberg( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t ); \
  template void Hessenberg \
  ( UpperOrLower uplo, ElementalMatrix<F>& A, ElementalMatrix<F>& t ); \
  template void hessenberg::ExplicitCondensed \
  ( UpperOrLower uplo, Matrix<F>& A ); \
  template void hessenberg::ExplicitCondensed \
  ( UpperOrLower uplo, ElementalMatrix<F>& A ); \
  template void hessenberg::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& H ); \
  template void hessenberg::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const ElementalMatrix<F>& A, const ElementalMatrix<F>& t, \
          ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
