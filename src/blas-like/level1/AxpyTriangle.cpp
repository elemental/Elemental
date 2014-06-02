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
void AxpyTriangle
( UpperOrLower uplo, T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("AxpyTriangle");
        if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
            X.Height() != Y.Height() )
            LogicError("Nonconformal AxpyTriangle");
    )
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

template<typename T,Dist U,Dist V>
void AxpyTriangle
( UpperOrLower uplo, T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
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

template<typename T>
void AxpyTriangle
( UpperOrLower uplo, T alpha, 
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

#define DIST_PROTO(T,U,V) \
  template void AxpyTriangle \
  ( UpperOrLower uplo, T       alpha, \
    const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B ); 

#define PROTO(T) \
  template void AxpyTriangle \
  ( UpperOrLower uplo, T alpha, const Matrix<T>& A, Matrix<T>& B ); \
  template void AxpyTriangle \
  ( UpperOrLower uplo, T alpha, \
    const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  DIST_PROTO(T,CIRC,CIRC); \
  DIST_PROTO(T,MC,  MR  ); \
  DIST_PROTO(T,MC,  STAR); \
  DIST_PROTO(T,MD,  STAR); \
  DIST_PROTO(T,MR,  MC  ); \
  DIST_PROTO(T,MR,  STAR); \
  DIST_PROTO(T,STAR,MC  ); \
  DIST_PROTO(T,STAR,MD  ); \
  DIST_PROTO(T,STAR,MR  ); \
  DIST_PROTO(T,STAR,STAR); \
  DIST_PROTO(T,STAR,VC  ); \
  DIST_PROTO(T,STAR,VR  ); \
  DIST_PROTO(T,VC  ,STAR); \
  DIST_PROTO(T,VR  ,STAR);

PROTO(Int);
PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
