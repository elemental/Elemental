/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_INDEXDEPENDENTMAP_HPP
#define EL_BLAS_INDEXDEPENDENTMAP_HPP

namespace El {

template<typename T>
void IndexDependentMap( Matrix<T>& A, function<T(Int,Int,const T&)> func )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A(i,j) = func(i,j,A(i,j));
}

template<typename T>
void IndexDependentMap
( AbstractDistMatrix<T>& A, function<T(Int,Int,const T&)> func )
{
    EL_DEBUG_CSE
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    auto& ALoc = A.Matrix();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            ALoc(iLoc,jLoc) = func(i,j,ALoc(iLoc,jLoc));
        }
    }
}

template<typename S,typename T>
void IndexDependentMap
( const Matrix<S>& A, Matrix<T>& B, function<T(Int,Int,const S&)> func )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B(i,j) = func(i,j,A(i,j));
}

template<typename S,typename T,Dist U,Dist V,DistWrap wrap>
void IndexDependentMap
( const DistMatrix<S,U,V,wrap>& A,
        DistMatrix<T,U,V,wrap>& B,
  function<T(Int,Int,const S&)> func )
{
    EL_DEBUG_CSE
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    B.AlignWith( A.DistData() );
    B.Resize( A.Height(), A.Width() );
    auto& ALoc = A.LockedMatrix();
    auto& BLoc = B.Matrix();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            BLoc(iLoc,jLoc) = func(i,j,ALoc(iLoc,jLoc));
        }
    }
}

template<typename S,typename T,Dist U,Dist V>
void IndexDependentMap
( const AbstractDistMatrix<S>& A,
        DistMatrix<T,U,V>& B,
  function<T(Int,Int,const S&)> func )
{
    EL_DEBUG_CSE
    if( A.Wrap() == ELEMENT && A.DistData() == B.DistData() )
    {
        auto& ACast = static_cast<const DistMatrix<T,U,V>&>(A);
        IndexDependentMap( ACast, B, func );
    }
    else
    {
        ElementalProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.rowConstrain = true;
        ctrl.root = B.Root();
        ctrl.colAlign = B.ColAlign();
        ctrl.rowAlign = B.RowAlign();

        DistMatrixReadProxy<S,S,U,V> ALikeBProx( A, ctrl );
        auto& ALikeB = ALikeBProx.GetLocked();
        IndexDependentMap( ALikeB, B, func );
    }
}

template<typename S,typename T,Dist U,Dist V>
void IndexDependentMap
( const AbstractDistMatrix<S>& A,
        DistMatrix<T,U,V,BLOCK>& B,
  function<T(Int,Int,const S&)> func )
{
    EL_DEBUG_CSE
    if( A.Wrap() == BLOCK && A.DistData() == B.DistData() )
    {
        auto& ACast = static_cast<const DistMatrix<T,U,V,BLOCK>&>(A);
        IndexDependentMap( ACast, B, func );
    }
    else
    {
        ProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.rowConstrain = true;
        ctrl.root = B.Root();
        ctrl.colAlign = B.ColAlign();
        ctrl.rowAlign = B.RowAlign();
        ctrl.blockHeight = B.BlockHeight();
        ctrl.blockWidth = B.BlockWidth();
        ctrl.rowCut = B.RowCut();
        ctrl.colCut = B.ColCut();

        DistMatrixReadProxy<S,S,U,V,BLOCK> ALikeBProx( A, ctrl );
        auto& ALikeB = ALikeBProx.GetLocked();
        IndexDependentMap( ALikeB, B, func );
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO_DIST(T,U,V) \
  EL_EXTERN template void IndexDependentMap \
  ( const DistMatrix<T,U,V>& A, \
          DistMatrix<T,U,V>& B, \
          function<T(Int,Int,const T&)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( const DistMatrix<T,U,V,BLOCK>& A, \
          DistMatrix<T,U,V,BLOCK>& B, \
          function<T(Int,Int,const T&)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( const AbstractDistMatrix<T>& A, \
          DistMatrix<T,U,V>& B, \
          function<T(Int,Int,const T&)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( const AbstractDistMatrix<T>& A, \
          DistMatrix<T,U,V,BLOCK>& B, \
          function<T(Int,Int,const T&)> func );

#define PROTO(T) \
  EL_EXTERN template void IndexDependentMap \
  ( Matrix<T>& A, \
    function<T(Int,Int,const T&)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( AbstractDistMatrix<T>& A, \
    function<T(Int,Int,const T&)> func ); \
  EL_EXTERN template void IndexDependentMap \
  ( const Matrix<T>& A, \
          Matrix<T>& B, \
          function<T(Int,Int,const T&)> func ); \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR)

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>
#undef PROTO_DIST

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_INDEXDEPENDENTMAP_HPP
