/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int ldim = X.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const T delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(delta);
            else
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(delta);
            else
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= delta;
        }
    }
}

template<typename TDiag,typename T,Dist U,Dist V>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<TDiag>& d, DistMatrix<T,U,V>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    if( side == LEFT )
    {
        DistMatrix<TDiag,U,GatheredDist<V>()> d_U_VGath( X.Grid() );
        d_U_VGath.AlignWith( X );
        Copy( d, d_U_VGath, READ_PROXY );
        DiagonalScale
        ( LEFT, orientation, d_U_VGath.LockedMatrix(), X.Matrix() );
    }
    else
    {
        DistMatrix<TDiag,V,GatheredDist<U>()> d_V_UGath( X.Grid() );
        d_V_UGath.AlignWith( X );
        Copy( d, d_V_UGath, READ_PROXY );
        DiagonalScale
        ( RIGHT, orientation, d_V_UGath.LockedMatrix(), X.Matrix() );
    }
}

template<typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    #define GUARD(CDIST,RDIST) \
        X.DistData().colDist == CDIST && X.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(X); \
        DiagonalScale( side, orientation, d, XCast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    #define GUARD(CDIST,RDIST) \
        X.DistData().colDist == CDIST && X.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<Complex<T>,CDIST,RDIST>&>(X); \
        DiagonalScale( side, orientation, d, XCast );
    #include "El/macros/GuardAndPayload.h"
}

#define DIST_PROTO(T,U,V) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<T,U,V>& X );

#define DIST_PROTO_REAL(T,U,V) \
  DIST_PROTO(T,U,V) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<Complex<T>,U,V>& X );

#define PROTO(T) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X ); \
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

#define PROTO_REAL(T) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<Complex<T>>& X ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& X ); \
  DIST_PROTO_REAL(T,CIRC,CIRC); \
  DIST_PROTO_REAL(T,MC,  MR  ); \
  DIST_PROTO_REAL(T,MC,  STAR); \
  DIST_PROTO_REAL(T,MD,  STAR); \
  DIST_PROTO_REAL(T,MR,  MC  ); \
  DIST_PROTO_REAL(T,MR,  STAR); \
  DIST_PROTO_REAL(T,STAR,MC  ); \
  DIST_PROTO_REAL(T,STAR,MD  ); \
  DIST_PROTO_REAL(T,STAR,MR  ); \
  DIST_PROTO_REAL(T,STAR,STAR); \
  DIST_PROTO_REAL(T,STAR,VC  ); \
  DIST_PROTO_REAL(T,STAR,VR  ); \
  DIST_PROTO_REAL(T,VC  ,STAR); \
  DIST_PROTO_REAL(T,VR  ,STAR);

#include "El/macros/Instantiate.h"

} // namespace El
