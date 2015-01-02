/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void IndexDependentMap( Matrix<T>& A, std::function<T(Int,Int,T)> func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func(i,j,A.Get(i,j)) );
}

template<typename T>
void IndexDependentMap
( AbstractDistMatrix<T>& A, std::function<T(Int,Int,T)> func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentMap"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, func(i,j,A.GetLocal(iLoc,jLoc)) );
        }
    }
}

template<typename T>
void IndexDependentMap
( AbstractBlockDistMatrix<T>& A, std::function<T(Int,Int,T)> func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentMap"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, func(i,j,A.GetLocal(iLoc,jLoc)) );
        }
    }
}

template<typename S,typename T>
void IndexDependentMap
( const Matrix<S>& A, Matrix<T>& B, std::function<T(Int,Int,S)> func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B.Set( i, j, func(i,j,A.Get(i,j)) );
}

template<typename S,typename T>
void IndexDependentMap
( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B, 
  std::function<T(Int,Int,S)> func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentMap"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    B.AlignWith( A.DistData() );
    B.Resize( A.Height(), A.Width() );
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            B.SetLocal( iLoc, jLoc, func(i,j,A.GetLocal(iLoc,jLoc)) );
        }
    }
}

template<typename S,typename T>
void IndexDependentMap
( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B, 
  std::function<T(Int,Int,S)> func )
{
    DEBUG_ONLY(CallStackEntry cse("IndexDependentMap"))
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    B.AlignWith( A.DistData() );
    B.Resize( A.Height(), A.Width() );
    for( Int jLoc=0; jLoc<nLoc; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLoc; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            B.SetLocal( iLoc, jLoc, func(i,j,A.GetLocal(iLoc,jLoc)) );
        }
    }
}

#define PROTO(T) \
  template void IndexDependentMap \
  ( Matrix<T>& A, std::function<T(Int,Int,T)> func ); \
  template void IndexDependentMap \
  ( AbstractDistMatrix<T>& A, std::function<T(Int,Int,T)> func ); \
  template void IndexDependentMap \
  ( AbstractBlockDistMatrix<T>& A, std::function<T(Int,Int,T)> func );

#define PROTO_TYPES(S,T) \
  template void IndexDependentMap \
  ( const Matrix<S>& A, Matrix<T>& B, std::function<T(Int,Int,S)> func ); \
  template void IndexDependentMap \
  ( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B, \
    std::function<T(Int,Int,S)> func ); \
  template void IndexDependentMap \
  ( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B, \
    std::function<T(Int,Int,S)> func );

#define PROTO_INT(T) \
  PROTO(T) \
  PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO(T) \
  PROTO_TYPES(Int,T) \
  PROTO_TYPES(T,T)

#define PROTO_COMPLEX(T) \
  PROTO(T) \
  PROTO_TYPES(Int,T) \
  PROTO_TYPES(Base<T>,T) 

#include "El/macros/Instantiate.h"

} // namespace El
