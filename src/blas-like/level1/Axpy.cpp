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
void Axpy( S alphaS, const Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("Axpy"))
    const T alpha = T(alphaS);
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( X.Height()==1 || X.Width()==1 )
    {
        const Int XLength = ( X.Width()==1 ? X.Height() : X.Width() );
        const Int XStride = ( X.Width()==1 ? 1          : X.LDim() );
        const Int YStride = ( Y.Width()==1 ? 1          : Y.LDim() );
        DEBUG_ONLY(
            const Int YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
            if( XLength != YLength )
                LogicError("Nonconformal Axpy");
        )
        blas::Axpy
        ( XLength, alpha, X.LockedBuffer(), XStride, Y.Buffer(), YStride );
    }
    else
    {
        DEBUG_ONLY(
            if( X.Height() != Y.Height() || X.Width() != Y.Width() )
                LogicError("Nonconformal Axpy");
        )
        if( X.Width() <= X.Height() )
        {
            for( Int j=0; j<X.Width(); ++j )
            {
                blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( Int i=0; i<X.Height(); ++i )
            {
                blas::Axpy
                ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                    Y.Buffer(i,0),       Y.LDim() );
            }
        }
    }
}

template<typename T,typename S,Dist U,Dist V,Dist W,Dist Z>
void Axpy( S alphaS, const DistMatrix<T,U,V>& X, DistMatrix<T,W,Z>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("Axpy");
        if( X.Grid() != Y.Grid() )
            LogicError("X and Y must be distributed over the same grid");
    )
    const T alpha = T(alphaS);
    if( U == V && W == Z && 
        X.ColAlign() == Y.ColAlign() && X.RowAlign() == Y.RowAlign() )
    {
        Axpy( alpha, X.LockedMatrix(), Y.Matrix() );
    }
    else
    {
        DistMatrix<T,W,Z> XCopy( X.Grid() );
        XCopy.AlignWith( Y );
        XCopy = X;
        Axpy( alpha, XCopy.LockedMatrix(), Y.Matrix() );
    }
}

template<typename T,typename S>
void Axpy( S alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Axpy"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(B); \
        Axpy( alpha, ACast, BCast );
    #include "El/core/NestedGuardAndPayload.h"
}

#define DIST_PROTO_INNER(T,S,U,V,W,Z) \
  template void Axpy \
  ( S alpha, const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

#define DIST_PROTO(T,S,U,V) \
  DIST_PROTO_INNER(T,S,U,V,CIRC,CIRC) \
  DIST_PROTO_INNER(T,S,U,V,MC,  MR  ) \
  DIST_PROTO_INNER(T,S,U,V,MC,  STAR) \
  DIST_PROTO_INNER(T,S,U,V,MD,  STAR) \
  DIST_PROTO_INNER(T,S,U,V,MR,  MC  ) \
  DIST_PROTO_INNER(T,S,U,V,MR,  STAR) \
  DIST_PROTO_INNER(T,S,U,V,STAR,MD  ) \
  DIST_PROTO_INNER(T,S,U,V,STAR,MR  ) \
  DIST_PROTO_INNER(T,S,U,V,STAR,STAR) \
  DIST_PROTO_INNER(T,S,U,V,STAR,VC  ) \
  DIST_PROTO_INNER(T,S,U,V,STAR,VR  ) \
  DIST_PROTO_INNER(T,S,U,V,VC,  STAR) \
  DIST_PROTO_INNER(T,S,U,V,VR,  STAR)

#define PROTO_TYPES(T,S) \
  template void Axpy( S alpha, const Matrix<T>& A, Matrix<T>& B ); \
  template void Axpy \
  ( S alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
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
