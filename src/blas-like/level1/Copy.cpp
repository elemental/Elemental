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
    auto convert = []( const Real alpha ) { return Complex<Real>(alpha); };
    EntrywiseMap( A, B, std::function<Complex<Real>(Real)>(convert) );
}

template<typename T,Dist U,Dist V>
void Copy( const AbstractDistMatrix<T>& A, DistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A; 
}

template<typename Real,Dist U,Dist V>
void Copy
( const AbstractDistMatrix<Real>& A, DistMatrix<Complex<Real>,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))

    DistMatrix<Real,U,V> BReal(A.Grid());
    if( A.Grid() == B.Grid() && A.ColDist() == U && A.RowDist() == V )
    {
        if( !B.RootConstrained() )
            B.SetRoot( A.Root() );
        if( !B.ColConstrained() )
            B.AlignCols( A.ColAlign() );
        if( !B.RowConstrained() )
            B.AlignRows( A.RowAlign() );
        if( A.Root() == B.Root() && 
            A.ColAlign() == B.ColAlign() && A.RowAlign() == B.RowAlign() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }
    BReal.AlignWith( B );
    BReal = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BReal.LockedMatrix(), B.Matrix() );
}

template<typename T,Dist U,Dist V>
void Copy( const AbstractBlockDistMatrix<T>& A, BlockDistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    // TODO: Add support for shallow copies
    B = A;
}

template<typename Real,Dist U,Dist V>
void Copy
( const AbstractBlockDistMatrix<Real>& A, 
  BlockDistMatrix<Complex<Real>,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))

    if( A.Grid() == B.Grid() && A.ColDist() == U && A.RowDist() == V )
    {
        if( !B.RootConstrained() )
            B.SetRoot( A.Root() );
        if( !B.ColConstrained() )
            B.AlignColsWith( A.DistData() );
        if( !B.RowConstrained() )
            B.AlignRowsWith( A.DistData() );
        if( A.Root() == B.Root() && 
            A.ColAlign() == B.ColAlign() && 
            A.RowAlign() == B.RowAlign() && 
            A.ColCut() == B.ColCut() &&
            A.RowCut() == B.RowCut() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }

    BlockDistMatrix<Real,U,V> BReal(A.Grid());
    BReal.AlignWith( B );
    BReal = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BReal.LockedMatrix(), B.Matrix() );
}

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(B); \
        Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
void Copy
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<BlockDistMatrix<T,CDIST,RDIST>&>(B); \
        Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename Real>
void Copy
( const AbstractDistMatrix<Real>& A, AbstractDistMatrix<Complex<Real>>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<Complex<Real>,CDIST,RDIST>&>(B); \
        Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename Real>
void Copy
( const AbstractBlockDistMatrix<Real>& A, 
        AbstractBlockDistMatrix<Complex<Real>>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& BCast = \
            dynamic_cast<BlockDistMatrix<Complex<Real>,CDIST,RDIST>&>(B); \
        Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

#define DIST_PROTO(T,U,V) \
  template void Copy( const AbstractDistMatrix<T>& A, DistMatrix<T,U,V>& B ); \
  template void Copy \
  ( const AbstractBlockDistMatrix<T>& A, BlockDistMatrix<T,U,V>& B );

#define DIST_PROTO_REAL(T,U,V) \
  DIST_PROTO(T,U,V) \
  template void Copy \
  ( const AbstractDistMatrix<T>& A, DistMatrix<Complex<T>,U,V>& B ); \
  template void Copy \
  ( const AbstractBlockDistMatrix<T>& A, BlockDistMatrix<Complex<T>,U,V>& B );

#define PROTO_BASE(T) \
  template void Copy( const Matrix<T>& A, Matrix<T>& B ); \
  template void Copy \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void Copy \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

#define PROTO(T) \
  PROTO_BASE(T) \
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
  PROTO_BASE(T) \
  template void Copy( const Matrix<T>& A, Matrix<Complex<T>>& B ); \
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
