/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Copy( const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

template<typename Real>
void Copy( const Matrix<Real>& A, Matrix<Complex<Real>>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B.Set( i, j, A.Get(i,j) );
}

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
void Copy( const DistMatrix<Real,U,V>& A, DistMatrix<Complex<Real>,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))

    if( U == W && V == Z )
    {
        if( !B.ColConstrained() )
            B.AlignCols( A.ColAlign() );
        if( !B.RowConstrained() )
            B.AlignRows( A.RowAlign() );
        if( A.ColAlign() == B.ColAlign() && A.RowAlign() == B.RowAlign() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }

    DistMatrix<Real,W,Z> BReal(A.Grid());
    BReal.AlignWith( B );
    BReal = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BReal.LockedMatrix(), B.Matrix() );
}

template<typename T,Dist U,Dist V,Dist W,Dist Z>
void Copy( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
void Copy
( const BlockDistMatrix<Real,U,V>& A, BlockDistMatrix<Complex<Real>,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))

    if( U == W && V == Z )
    {
        if( !B.ColConstrained() )
            B.AlignColsWith( A.DistData() );
        if( !B.RowConstrained() )
            B.AlignRowsWith( A.DistData() );
        if( A.ColAlign() == B.ColAlign() && 
            A.RowAlign() == B.RowAlign() && 
            A.ColCut() == B.ColCut() &&
            A.RowCut() == B.RowCut() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }

    BlockDistMatrix<Real,W,Z> BReal(A.Grid());
    BReal.AlignWith( B );
    BReal = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BReal.LockedMatrix(), B.Matrix() );
}

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const DistMatrix<T,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(B); \
        Copy( ACast, BCast );
    #include "El/core/NestedGuardAndPayload.h"
}

template<typename T>
void Copy( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const BlockDistMatrix<T,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<BlockDistMatrix<T,CDIST,RDIST>&>(B); \
        Copy( ACast, BCast );
    #include "El/core/NestedGuardAndPayload.h"
}

template<typename Real>
void Copy
( const AbstractDistMatrix<Real>& A, AbstractDistMatrix<Complex<Real>>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const DistMatrix<Real,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<Complex<Real>,CDIST,RDIST>&>(B); \
        Copy( ACast, BCast );
    #include "El/core/NestedGuardAndPayload.h"
}

template<typename Real>
void Copy
( const AbstractBlockDistMatrix<Real>& A, 
        AbstractBlockDistMatrix<Complex<Real>>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) \
        A.DistData().colDist == CDIST && A.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        B.DistData().colDist == CDIST && B.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<const BlockDistMatrix<Real,CDIST,RDIST>&>(A);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& BCast = \
            dynamic_cast<BlockDistMatrix<Complex<Real>,CDIST,RDIST>&>(B); \
        Copy( ACast, BCast );
    #include "El/core/NestedGuardAndPayload.h"
}

#define DIST_PROTO_INNER(T,U,V,W,Z) \
  template void Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B ); \
  template void Copy \
  ( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B );

#define DIST_PROTO_INNER_REAL(T,U,V,W,Z) \
  DIST_PROTO_INNER(T,U,V,W,Z) \
  template void Copy \
  ( const DistMatrix<T,U,V>& A, DistMatrix<Complex<T>,W,Z>& B ); \
  template void Copy \
  ( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<Complex<T>,W,Z>& B );

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
  DIST_PROTO_INNER(T,U,V,VR,  STAR);

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

#define PROTO(T) \
  template void Copy( const Matrix<T>& A, Matrix<T>& B ); \
  template void Copy \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void Copy \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B ); \
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
  template void Copy( const Matrix<T>& A, Matrix<T>& B ); \
  template void Copy( const Matrix<T>& A, Matrix<Complex<T>>& B ); \
  template void Copy \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void Copy \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B ); \
  template void Copy \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Complex<T>>& B ); \
  template void Copy \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<Complex<T>>& B ); \
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
