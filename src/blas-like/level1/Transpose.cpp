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
void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( n, m );
    if( conjugate )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B.Set(j,i,Conj(A.Get(i,j)));
    }
    else
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B.Set(j,i,A.Get(i,j));
    }
}

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Transpose
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    if( U == Z && V == W && 
        (A.ColAlign()==B.RowAlign() || !B.RowConstrained()) &&
        (A.RowAlign()==B.ColAlign() || !B.ColConstrained()) )
    {
        B.Align( A.RowAlign(), A.ColAlign() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        DistMatrix<T,Z,W> C( B.Grid() );
        C.AlignRowsWith( B );
        C.AlignColsWith( B );
        C = A;
        B.Resize( A.Width(), A.Height() );
        Transpose( C.LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Transpose
( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    // TODO: Member function to quickly query this information?!?
    //       Perhaps A.ColsAlignedWithRows( B )?
    //               A.ColsAlignedWithCols( B )?
    const bool transposeColAligned = 
        U               == Z              &&
        A.BlockHeight() == B.BlockWidth() &&
        A.ColAlign()    == B.RowAlign()   &&
        A.ColCut()      == B.RowCut();
    const bool transposeRowAligned =
        V              == W               &&
        A.BlockWidth() == B.BlockHeight() &&
        A.RowAlign()   == B.ColAlign()    &&
        A.RowCut()     == B.ColCut();
    if( (transposeColAligned || !B.RowConstrained()) &&
        (transposeRowAligned || !B.ColConstrained()) )
    {
        B.Align
        ( A.BlockWidth(), A.BlockHeight(), 
          A.RowAlign(), A.ColAlign(), A.RowCut(), A.ColCut() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        BlockDistMatrix<T,Z,W> C( B.Grid() );
        C.AlignRowsWith( B );
        C.AlignColsWith( B );
        C = A;
        B.Resize( A.Width(), A.Height() );
        Transpose( C.LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(B); \
        Transpose( ACast, BCast, conjugate );
    #include "El/core/NestedGuardAndPayload.h"
}

template<typename T>
void Transpose
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B, 
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const BlockDistMatrix<T,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<BlockDistMatrix<T,CDIST,RDIST>&>(B); \
        Transpose( ACast, BCast, conjugate );
    #include "El/core/NestedGuardAndPayload.h"
}

#define DIST_PROTO_INNER(T,U,V,W,Z) \
  template void Transpose \
  ( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B, bool conjugate ); \
  template void Transpose \
  ( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B, \
    bool conjugate );

#define DIST_PROTO(T,U,V) \
  DIST_PROTO_INNER(T,U,V,CIRC,CIRC); \
  DIST_PROTO_INNER(T,U,V,MC,  MR  ); \
  DIST_PROTO_INNER(T,U,V,MC,  STAR); \
  DIST_PROTO_INNER(T,U,V,MD,  STAR); \
  DIST_PROTO_INNER(T,U,V,MR,  MC  ); \
  DIST_PROTO_INNER(T,U,V,MR,  STAR); \
  DIST_PROTO_INNER(T,U,V,STAR,MD  ); \
  DIST_PROTO_INNER(T,U,V,STAR,MR  ); \
  DIST_PROTO_INNER(T,U,V,STAR,STAR); \
  DIST_PROTO_INNER(T,U,V,STAR,VC  ); \
  DIST_PROTO_INNER(T,U,V,STAR,VR  ); \
  DIST_PROTO_INNER(T,U,V,VC,  STAR); \
  DIST_PROTO_INNER(T,U,V,VR,  STAR); \

#define PROTO(T) \
  template void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<T>& B, bool conjugate ); \
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
