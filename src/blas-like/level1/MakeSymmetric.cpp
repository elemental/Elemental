/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T>
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    Matrix<T> d;
    A.GetDiagonal( d );
    if( conjugate )
        MakeReal( d );

    if( uplo == LOWER )
        MakeTriangular( LOWER, A );
    else
        MakeTriangular( UPPER, A );
    SetDiagonal( A, T(0) );
    Matrix<T> ATrans;
    Transpose( A, ATrans, conjugate );
    Axpy( T(1), ATrans, A );

    A.SetDiagonal( d );
}

template<typename T,Dist U,Dist V>
void MakeSymmetric
( UpperOrLower uplo, DistMatrix<T,U,V>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    auto d = A.GetDiagonal();
    if( conjugate )
        MakeReal( d );

    if( uplo == LOWER )
        MakeTriangular( LOWER, A );
    else
        MakeTriangular( UPPER, A );
    SetDiagonal( A, T(0) );

    const Grid& g = A.Grid();
    DistMatrix<T,U,V> ATrans(g);
    Transpose( A, ATrans, conjugate );
    Axpy( T(1), ATrans, A );

    A.SetDiagonal( d );
}

template<typename F>
void MakeSymmetric
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<DistMatrix<F,CDIST,RDIST>&>(A); \
        MakeSymmetric( uplo, ACast, conjugate );
    #include "El/core/GuardAndPayload.h"
}

#define DIST_PROTO(F,U,V) \
  template void MakeSymmetric \
  ( UpperOrLower uplo, DistMatrix<F,U,V>& A, bool conjugate );

#define PROTO(F) \
  template void MakeSymmetric \
  ( UpperOrLower uplo, Matrix<F>& A, bool conjugate ); \
  template void MakeSymmetric \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool conjugate ); \
  DIST_PROTO(F,CIRC,CIRC) \
  DIST_PROTO(F,MC,  MR  ) \
  DIST_PROTO(F,MC,  STAR) \
  DIST_PROTO(F,MD,  STAR) \
  DIST_PROTO(F,MR,  MC  ) \
  DIST_PROTO(F,MR,  STAR) \
  DIST_PROTO(F,STAR,MC  ) \
  DIST_PROTO(F,STAR,MD  ) \
  DIST_PROTO(F,STAR,MR  ) \
  DIST_PROTO(F,STAR,STAR) \
  DIST_PROTO(F,STAR,VC  ) \
  DIST_PROTO(F,STAR,VR  ) \
  DIST_PROTO(F,VC,  STAR) \
  DIST_PROTO(F,VR,  STAR)

PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
