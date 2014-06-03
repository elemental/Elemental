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
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  const Matrix<TDiag>& d, Matrix<T>& A, Int offset )
{
    DEBUG_ONLY(
        CallStackEntry cse("DiagonalScaleTrapezoid");
        if( side==LEFT && (d.Height()!=A.Height() || d.Width()!=1) )
            LogicError("d should have been a vector of the height of A");
        if( side==RIGHT && (d.Height()!=A.Width() || d.Width()!=1) )
            LogicError("d should have been a vector of the width of A");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    const Int ldim = A.LDim();
    T* ABuf = A.Buffer();
    const bool conjugate = ( orientation==ADJOINT );

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    if( uplo == LOWER && side == LEFT )
    {
        // Scale from the left up to the diagonal
        for( Int i=iOff; i<m; ++i )
        {
            const Int k = i-iOff;
            const Int j = k+jOff;
            const TDiag alpha = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            blas::Scal( Min(j+1,n), alpha, &ABuf[i], ldim );
        }
    }
    else if( uplo == UPPER && side == LEFT )
    {
        // Scale from the diagonal to the right
        for( Int i=0; i<iOff+diagLength; ++i )
        {
            const Int k = i-iOff;
            const Int j = k+jOff;
            const Int jLeft = Max(j,0);
            const TDiag alpha = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            blas::Scal( n-jLeft, alpha, &ABuf[i+jLeft*ldim], ldim );
        }
    }
    else if( uplo == LOWER && side == RIGHT )
    {
        // Scale from the diagonal downwards
        for( Int j=0; j<jOff+diagLength; ++j )
        {
            const Int k = j-jOff;
            const Int i = k+iOff;
            const Int iTop = Max(i,0);
            const TDiag alpha = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            blas::Scal( m-iTop, alpha, &ABuf[iTop+j*ldim], 1 );
        }
    }
    else /* uplo == UPPER && side == RIGHT */
    {
        // Scale downward to the diagonal
        for( Int j=jOff; j<n; ++j )
        {
            const Int k = j-jOff;
            const Int i = k+iOff;
            const TDiag alpha = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            blas::Scal( Min(i+1,m), alpha, &ABuf[j*ldim], 1 );
        }
    }
}

template<typename TDiag,typename T,Dist U,Dist V,Dist W,Dist Z>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& A, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScaleTrapezoid"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    const bool conjugate = ( orientation==ADJOINT );

    const Int diagLength = A.DiagonalLength(offset);
    const Int ldim = A.LDim();
    T* ABuf = A.Buffer();

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    if( side == LEFT )
    {
        const Dist ZGath = GatheredDist<Z>();
        DistMatrix<TDiag,W,ZGath> d_W_ZGath( A.Grid() );
        if( U == W && V == STAR && d.ColAlign() == A.ColAlign() )
        {
            d_W_ZGath = LockedView( d );
        }
        else
        {
            d_W_ZGath.AlignWith( A );
            d_W_ZGath = d;
        }

        if( uplo == LOWER )
        {
            // Scale from the left up to the diagonal
            for( Int iLoc=0; iLoc<mLoc; ++iLoc )            
            {
                const Int i = A.GlobalRow(iLoc);
                if( i >= iOff )
                {
                    const Int k = i-iOff;
                    const Int j = k+jOff;
                    const Int width = Min(j+1,n);
                    const Int localWidth = A.LocalColOffset(width);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_W_ZGath.GetLocal(iLoc,0))
                                    : d_W_ZGath.GetLocal(iLoc,0) );
                    blas::Scal( localWidth, alpha, &ABuf[iLoc], ldim );
                }
            }
        }
        else
        {
            // Scale from the diagonal to the right
            for( Int iLoc=0; iLoc<mLoc; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                if( i < iOff+diagLength )
                {
                    const Int k = i-iOff;
                    const Int j = k+jOff;
                    const Int jLeft = Max(j,0);
                    const Int jLeftLoc = A.LocalColOffset(jLeft);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_W_ZGath.GetLocal(iLoc,0))
                                    : d_W_ZGath.GetLocal(iLoc,0) );
                    blas::Scal
                    ( nLoc-jLeftLoc, alpha, &ABuf[iLoc+jLeftLoc*ldim], ldim );
                }
            }
        }    
    }
    else
    {
        const Dist WGath = GatheredDist<W>();
        DistMatrix<TDiag,Z,WGath> d_Z_WGath( A.Grid() );
        if( U == Z && V == STAR && d.ColAlign() == A.RowAlign() )
        {
            d_Z_WGath = LockedView( d );
        }
        else
        {
            d_Z_WGath.AlignWith( A );
            d_Z_WGath = d;
        }

        if( uplo == LOWER )
        {
            // Scale from the diagonal downwards
            for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                if( j < jOff+diagLength )
                {
                    const Int k = j-jOff;
                    const Int i = k+iOff;
                    const Int iTop = Max(i,0);
                    const Int iTopLoc = A.LocalRowOffset(iTop);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_Z_WGath.GetLocal(jLoc,0))
                                    : d_Z_WGath.GetLocal(jLoc,0) );
                    blas::Scal
                    ( mLoc-iTopLoc, alpha, &ABuf[iTopLoc+jLoc*ldim], 1 );
                }
            }
        }
        else 
        {
            // Scale downward to the diagonal
            for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                if( j >= jOff )
                {
                    const Int k = j-jOff;
                    const Int i = k+iOff;
                    const Int height = Min(i+1,m);
                    const Int localHeight = A.LocalRowOffset(height);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d_Z_WGath.GetLocal(jLoc,0))
                                    : d_Z_WGath.GetLocal(jLoc,0) );
                    blas::Scal( localHeight, alpha, &ABuf[jLoc*ldim], 1 );
                }
            }
        }
    }
}

template<typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X, Int offset )
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
        DiagonalScaleTrapezoid( side, uplo, orientation, dCast, XCast, offset );
    #include "El/core/NestedGuardAndPayload.h"
}

template<typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& X, 
  Int offset )
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
        DiagonalScaleTrapezoid \
        ( side, uplo, orientation, dCast, XCast, offset );
    #include "El/core/NestedGuardAndPayload.h"
}

#define DIST_PROTO_INNER_BASE(T,U,V,W,Z) \
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X, Int offset );

#define DIST_PROTO_INNER_REAL(T,U,V,W,Z) \
  DIST_PROTO_INNER_BASE(T,U,V,W,Z) \
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const DistMatrix<T,U,V>& d, DistMatrix<Complex<T>,W,Z>& X, Int offset );

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
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X, Int offset ); \
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X, Int offset ); \
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
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X, Int offset ); \
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const Matrix<T>& d, Matrix<Complex<T>>& X, Int offset ); \
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X, Int offset ); \
  template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& X, \
    Int offset ); \
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
