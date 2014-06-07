/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./Hessenberg/L.hpp"
#include "./Hessenberg/U.hpp"
#include "./Hessenberg/ApplyQ.hpp"

namespace El {

template<typename F>
void Hessenberg( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    if( uplo == UPPER )
        hessenberg::U( A, t );
    else
        hessenberg::L( A, t );
}

template<typename F> 
void Hessenberg
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    if( uplo == UPPER )
        hessenberg::U( A, t );
    else
        hessenberg::L( A, t );
}

template<typename F>
void Hessenberg( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    Matrix<F> t;
    Hessenberg( uplo, A, t );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F> 
void Hessenberg( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Hessenberg"))
    DistMatrix<F,STAR,STAR> t(A.Grid());
    Hessenberg( uplo, A, t );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

#define PROTO(F) \
  template void Hessenberg( UpperOrLower uplo, Matrix<F>& A ); \
  template void Hessenberg( UpperOrLower uplo, DistMatrix<F>& A ); \
  template void Hessenberg( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t ); \
  template void Hessenberg \
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t ); \
  template void hessenberg::ApplyQ \
  ( UpperOrLower uplo, LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& H ); \
  template void hessenberg::ApplyQ \
  ( UpperOrLower uplo, LeftOrRight side, Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& H ); \
  template void hessenberg::ApplyQ \
  ( UpperOrLower uplo, LeftOrRight side, Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& H );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
