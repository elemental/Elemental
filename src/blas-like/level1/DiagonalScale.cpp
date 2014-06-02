/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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

template<typename TDiag,typename T,Dist U,Dist V,Dist W,Dist Z>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    if( side == LEFT )
    {
        const Dist ZGath = GatheredDist<Z>();
        DistMatrix<TDiag,W,ZGath> d_W_ZGath( X.Grid() );
        if( U == W && V == STAR && d.ColAlign() == X.ColAlign() )
        {
            d_W_ZGath = LockedView( d );
        }
        else
        {
            d_W_ZGath.AlignWith( X );
            d_W_ZGath = d;
        }
        DiagonalScale
        ( LEFT, orientation, d_W_ZGath.LockedMatrix(), X.Matrix() );
    }
    else
    {
        const Dist WGath = GatheredDist<W>();
        DistMatrix<TDiag,Z,WGath> d_Z_WGath( X.Grid() );
        if( U == Z && V == STAR && d.ColAlign() == X.RowAlign() )
        {
            d_Z_WGath = LockedView( d );
        }
        else
        {
            d_Z_WGath.AlignWith( X );
            d_Z_WGath = d;
        }
        DiagonalScale
        ( RIGHT, orientation, d_Z_WGath.LockedMatrix(), X.Matrix() );
    }
}

template<typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    #define GUARD(CDIST,RDIST) \
        d.DistData().colDist == CDIST && d.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        X.DistData().colDist == CDIST && X.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& dCast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(d);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(X); \
        DiagonalScale( side, orientation, dCast, XCast );
    #include "El/core/NestedGuardAndPayload.h"
}

template<typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    #define GUARD(CDIST,RDIST) \
        d.DistData().colDist == CDIST && d.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        X.DistData().colDist == CDIST && X.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& dCast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(d);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<Complex<T>,CDIST,RDIST>&>(X); \
        DiagonalScale( side, orientation, dCast, XCast );
    #include "El/core/NestedGuardAndPayload.h"
}

#define DIST_PROTO_INNER_BASE(T,U,V,W,Z) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X );

#define DIST_PROTO_INNER_REAL(T,U,V,W,Z) \
  DIST_PROTO_INNER_BASE(T,U,V,W,Z) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMatrix<T,U,V>& d, DistMatrix<Complex<T>,W,Z>& X );

#define DIST_PROTO_BASE(T,U,V) \
  DIST_PROTO_INNER_BASE(T,U,V,CIRC,CIRC); \
  DIST_PROTO_INNER_BASE(T,U,V,MC,  MR  ); \
  DIST_PROTO_INNER_BASE(T,U,V,MC,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,MD,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,MR,  MC  ); \
  DIST_PROTO_INNER_BASE(T,U,V,MR,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,MD  ); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,MR  ); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,VC  ); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,VR  ); \
  DIST_PROTO_INNER_BASE(T,U,V,VC,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,VR,  STAR);

#define DIST_PROTO_REAL(T,U,V) \
  DIST_PROTO_INNER_REAL(T,U,V,CIRC,CIRC); \
  DIST_PROTO_INNER_REAL(T,U,V,MC,  MR  ); \
  DIST_PROTO_INNER_REAL(T,U,V,MC,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,MD,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,MR,  MC  ); \
  DIST_PROTO_INNER_REAL(T,U,V,MR,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,MD  ); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,MR  ); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,VC  ); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,VR  ); \
  DIST_PROTO_INNER_REAL(T,U,V,VC,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,VR,  STAR);

#define PROTO_BASE(T) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X ); \
  DIST_PROTO_BASE(T,CIRC,CIRC); \
  DIST_PROTO_BASE(T,MC,  MR  ); \
  DIST_PROTO_BASE(T,MC,  STAR); \
  DIST_PROTO_BASE(T,MD,  STAR); \
  DIST_PROTO_BASE(T,MR,  MC  ); \
  DIST_PROTO_BASE(T,MR,  STAR); \
  DIST_PROTO_BASE(T,STAR,MC  ); \
  DIST_PROTO_BASE(T,STAR,MD  ); \
  DIST_PROTO_BASE(T,STAR,MR  ); \
  DIST_PROTO_BASE(T,STAR,STAR); \
  DIST_PROTO_BASE(T,STAR,VC  ); \
  DIST_PROTO_BASE(T,STAR,VR  ); \
  DIST_PROTO_BASE(T,VC  ,STAR); \
  DIST_PROTO_BASE(T,VR  ,STAR);

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

PROTO_BASE(Int);
PROTO_REAL(float);
PROTO_REAL(double);
PROTO_BASE(Complex<float>);
PROTO_BASE(Complex<double>);

} // namespace El
