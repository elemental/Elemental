/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T,typename S>
void AxpyTriangle
( UpperOrLower uplo, S alphaS, const Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("AxpyTriangle");
        if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
            X.Height() != Y.Height() )
            LogicError("Nonconformal AxpyTriangle");
    )
    const T alpha = T(alphaS);
    if( uplo == UPPER )
    {
        for( Int j=0; j<X.Width(); ++j )
            blas::Axpy( j+1, alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
    else
    {
        const Int n = X.Height();
        for( Int j=0; j<X.Width(); ++j )
            blas::Axpy( n-j, alpha, X.LockedBuffer(j,j), 1, Y.Buffer(j,j), 1 );
    }
}

template<typename T,typename S,Dist U,Dist V>
void AxpyTriangle
( UpperOrLower uplo, S alphaS, 
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("AxpyTriangle");
        if( X.Grid() != Y.Grid() )
            LogicError
            ("X and Y must be distributed over the same grid");
        if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
            X.Height() != Y.Height() )
            LogicError("Nonconformal AxpyTriangle");
    )
    const T alpha = T(alphaS);
    if( X.ColAlign() == Y.ColAlign() && X.RowAlign() == Y.RowAlign() )
    {
        const Int localHeight = X.LocalHeight();
        const Int localWidth = X.LocalWidth();
        const T* XBuffer = X.LockedBuffer();
        T* YBuffer = Y.Buffer();
        const Int XLDim = X.LDim();
        const Int YLDim = Y.LDim();
        if( uplo == UPPER )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = X.GlobalCol(jLoc);
                const Int localHeightAbove = X.LocalRowOffset(j+1);
                blas::Axpy
                ( localHeightAbove, alpha, 
                  &XBuffer[jLoc*XLDim], 1, &YBuffer[jLoc*YLDim], 1 );
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = X.GlobalCol(jLoc);
                const Int localHeightAbove = X.LocalRowOffset(j);
                const Int localHeightBelow = localHeight - localHeightAbove;
                blas::Axpy
                ( localHeightBelow, alpha, 
                  &XBuffer[localHeightAbove+jLoc*XLDim], 1,
                  &YBuffer[localHeightAbove+jLoc*YLDim], 1 );
            }
        }
    }
    else
    {
        DistMatrix<T,U,V> XCopy( X.Grid() );
        XCopy.AlignWith( Y );
        XCopy = X;
        AxpyTriangle( uplo, alpha, XCopy, Y );
    }
}

template<typename T,typename S>
void AxpyTriangle
( UpperOrLower uplo, S alpha, 
  const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("AxpyTriangle"))
    if( A.DistData().colDist != B.DistData().colDist || 
        A.DistData().rowDist != B.DistData().rowDist )
        RuntimeError("Distributions of A and B must match");
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A); \
        auto& BCast = dynamic_cast<      DistMatrix<T,CDIST,RDIST>&>(B); \
        AxpyTriangle( uplo, alpha, ACast, BCast );
    #include "El/core/GuardAndPayload.h"
}

#define DIST_PROTO(T,S,U,V) \
  template void AxpyTriangle \
  ( UpperOrLower uplo, S alpha, \
    const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B ); 

#define PROTO_TYPES(T,S) \
  template void AxpyTriangle \
  ( UpperOrLower uplo, S alpha, const Matrix<T>& A, Matrix<T>& B ); \
  template void AxpyTriangle \
  ( UpperOrLower uplo, S alpha, \
    const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  DIST_PROTO(T,S,CIRC,CIRC) \
  DIST_PROTO(T,S,MC,  MR  ) \
  DIST_PROTO(T,S,MC,  STAR) \
  DIST_PROTO(T,S,MD,  STAR) \
  DIST_PROTO(T,S,MR,  MC  ) \
  DIST_PROTO(T,S,MR,  STAR) \
  DIST_PROTO(T,S,STAR,MC  ) \
  DIST_PROTO(T,S,STAR,MD  ) \
  DIST_PROTO(T,S,STAR,MR  ) \
  DIST_PROTO(T,S,STAR,STAR) \
  DIST_PROTO(T,S,STAR,VC  ) \
  DIST_PROTO(T,S,STAR,VR  ) \
  DIST_PROTO(T,S,VC  ,STAR) \
  DIST_PROTO(T,S,VR  ,STAR)

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T)

#define PROTO_CPX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_TYPES(T,T)

PROTO_INT(Int)
PROTO_REAL(float)
PROTO_REAL(double)
PROTO_CPX(Complex<float>)
PROTO_CPX(Complex<double>)

} // namespace El
